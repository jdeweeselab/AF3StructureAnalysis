#!/usr/bin/env python3
import argparse
import csv
import json
import re
import shlex
import time
import traceback
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["pdf.fonttype"] = 42

HERE = Path(__file__).resolve().parent
THREE_TO_ONE = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V","SEP":"s","TPO":"t","PTR":"y"}
SS_COLORS = {"H":"#d62728","S":"#1f77b4","L":"#c7c7c7"}


def log(msg, fh):
    line=f"[{time.strftime('%H:%M:%S')}] {msg}"; print(line,flush=True); fh.write(line+"\n"); fh.flush()


def read_cif_backbone(path):
    cols, rows, in_loop, reading = [], [], False, False
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            s=line.strip()
            if s=="loop_":
                if reading: break
                in_loop,cols=True,[]; continue
            if in_loop and s.startswith("_atom_site."): cols.append(s.split()[0][len("_atom_site."):]); continue
            if in_loop and cols and not reading:
                if s.startswith("_"): in_loop,cols=False,[]; continue
                reading=True
            if reading:
                if not s or s.startswith("#") or s.startswith("_") or s=="loop_": break
                rows.append(shlex.split(s) if ('"' in s or "'" in s) else s.split())
    if not cols: raise ValueError(f"no _atom_site loop in {path}")
    ix={c:i for i,c in enumerate(cols)}
    def col(*names):
        for n in names:
            if n in ix:return ix[n]
        raise KeyError(names)
    i_atom=col("label_atom_id","auth_atom_id"); i_res=col("label_comp_id","auth_comp_id"); i_ch=col("auth_asym_id","label_asym_id")
    i_seq=col("label_seq_id","auth_seq_id"); i_x=col("Cartn_x"); i_y=col("Cartn_y"); i_z=col("Cartn_z"); i_b=col("B_iso_or_equiv")
    i_model=ix.get("pdbx_PDB_model_num"); chains={}; first_model=None
    for r in rows:
        if i_model is not None:
            if first_model is None:first_model=r[i_model]
            elif r[i_model]!=first_model:break
        if r[i_seq] in (".","?"):continue
        resname=r[i_res]; atom=r[i_atom].strip('"')
        if atom not in ("N","CA","C","O"):continue
        ch=chains.setdefault(r[i_ch],{}); res=ch.setdefault(int(r[i_seq]),{"name":resname})
        res[atom]=(float(r[i_x]),float(r[i_y]),float(r[i_z]))
        if atom=="CA":res["plddt"]=float(r[i_b])
    out={}
    for c,residues in chains.items():
        nums=sorted(residues)
        if len(nums)<20 or sum("CA" in residues[n] for n in nums)<0.9*len(nums):continue
        arr={a:np.full((len(nums),3),np.nan) for a in ("N","CA","C","O")}
        for k,n in enumerate(nums):
            for a in arr:
                if a in residues[n]:arr[a][k]=residues[n][a]
        out[c]=dict(resnum=np.array(nums),resname=[residues[n]["name"] for n in nums],plddt=np.array([residues[n].get("plddt",np.nan) for n in nums]),**arr)
    return out


def assign_ss(bb):
    N,C,O=bb["N"],bb["C"],bb["O"]; n=len(N)
    H=np.full_like(N,np.nan); co=C[:-1]-O[:-1]
    with np.errstate(invalid="ignore",divide="ignore"):co/=np.linalg.norm(co,axis=1,keepdims=True)
    H[1:]=N[1:]+co
    brk=np.linalg.norm(N[1:]-C[:-1],axis=1)>2.5; H[1:][brk]=np.nan
    pro=np.array([r=="PRO" for r in bb["resname"]]); H[pro]=np.nan
    def d(a,b):return np.linalg.norm(a[:,None,:]-b[None,:,:],axis=2)
    with np.errstate(invalid="ignore",divide="ignore"):
        E=0.084*332*(1/d(O,N)+1/d(C,H)-1/d(O,H)-1/d(C,N))
    hb=np.nan_to_num(E,nan=0.0)<-0.5; idx=np.arange(n)
    hb[np.abs(idx[:,None]-idx[None,:])<2]=False
    def hbond(i,j):
        ok=(i>=0)&(i<n)&(j>=0)&(j<n); out=np.zeros(np.broadcast(i,j).shape,dtype=bool)
        out[ok]=hb[np.broadcast_to(i,out.shape)[ok],np.broadcast_to(j,out.shape)[ok]]; return out
    helix=np.zeros(n,dtype=bool)
    for t in (4,3,5):
        turn=np.zeros(n,dtype=bool); turn[:n-t]=hb[idx[:n-t],idx[t:]]
        start=np.zeros(n,dtype=bool); start[1:]=turn[1:]&turn[:-1]
        for i in np.where(start)[0]:helix[i:i+t]=True
    I,J=np.meshgrid(idx,idx,indexing="ij"); far=np.abs(I-J)>2
    par=(hbond(I-1,J)&hbond(J,I+1))|(hbond(J-1,I)&hbond(I,J+1))
    anti=(hbond(I,J)&hbond(J,I))|(hbond(I-1,J+1)&hbond(J-1,I+1))
    strand=((par|anti)&far).any(axis=1)
    ss=np.full(n,"L",dtype="<U1"); ss[strand]="S"; ss[helix]="H"; return ss


def find_cifs(folder):
    out=[]
    for p in sorted(folder.rglob("*.cif")):
        m=re.search(r"(.*)_model_(\d+)(?:_\d+)?\.cif$",p.name)
        if m:out.append((f"{m.group(1)}_model{m.group(2)}",p));continue
        if re.search(r"seed-\d+_sample-\d+",p.parent.name):out.append((f"{p.parent.parent.name}/{p.parent.name}",p))
    seen=set();uniq=[]
    for name,p in out:
        if name not in seen:seen.add(name);uniq.append((name,p))
    return uniq


def run_condition(cond,data_dir,out_root,args,fh):
    name=cond["name"]; folder=data_dir/cond["folder"]
    out=out_root/name/"secondary_structure"
    if not folder.exists():log(f"SKIP {name}: folder not found ({folder})",fh);return None
    if (out/"summary.json").exists() and not args.force:return json.loads((out/"summary.json").read_text())
    cifs=find_cifs(folder)
    if not cifs:log(f"SKIP {name}: no model CIF files found in {folder}",fh);return None
    log(f"{name}: {len(cifs)} model CIF files",fh)
    if len(cifs)!=15:log(f"  WARNING {name}: expected 15 models, found {len(cifs)}",fh)
    out.mkdir(parents=True,exist_ok=True); lo,hi=(int(x) for x in cond["ctd"].split("-"))
    per_chain={};model_rows=[]
    for mname,path in cifs:
        chains=read_cif_backbone(path)
        for c,bb in chains.items():
            ss=assign_ss(bb); pc=per_chain.setdefault(c,dict(resnum=bb["resnum"],resname=bb["resname"],ss=[],plddt=[],models=[]))
            if len(bb["resnum"])!=len(pc["resnum"]) or (bb["resnum"]!=pc["resnum"]).any():raise ValueError(f"{mname} chain {c}: residue numbering differs")
            pc["ss"].append(ss);pc["plddt"].append(bb["plddt"]);pc["models"].append(mname)
            ctd=(bb["resnum"]>=lo)&(bb["resnum"]<=hi)
            model_rows.append({"model":mname,"chain":c,"CTD_pct_helix":round(100*float((ss[ctd]=="H").mean()),1),"CTD_pct_strand":round(100*float((ss[ctd]=="S").mean()),1),"CTD_mean_pLDDT":round(float(np.nanmean(bb["plddt"][ctd])),1),"core_mean_pLDDT":round(float(np.nanmean(bb["plddt"][~ctd])),1)})
    tables={};summ={"condition":name,"protein":cond["protein"],"ptm":cond["ptm"],"dna":cond["dna"],"n_models":len(cifs),"ctd":f"{lo}-{hi}"};pooled={"H":[],"L":[],"S":[]};cons_all=[]
    for c,pc in per_chain.items():
        ss=np.array(pc["ss"]);pl=np.array(pc["plddt"]);nm=ss.shape[0]
        counts={k:(ss==k).sum(0) for k in ("H","L","S")};stack=np.stack([counts["H"],counts["L"],counts["S"]]);agree=stack.max(0);cons=np.array(["H","L","S"])[stack.argmax(0)];tie=(stack==agree).sum(0)>1
        rows=[]
        for k,rn in enumerate(pc["resnum"]):
            rows.append(dict(residue=int(rn),aa=THREE_TO_ONE.get(pc["resname"][k],"X"),**{f"m{m+1}":ss[m,k] for m in range(nm)},n_H=int(counts["H"][k]),n_L=int(counts["L"][k]),n_S=int(counts["S"][k]),consensus=cons[k],agreement=int(agree[k]),tie="tie" if tie[k] else "",pct_helix=round(100*counts["H"][k]/nm,1),pLDDT_mean=round(float(np.nanmean(pl[:,k])),1),pLDDT_sd=round(float(np.nanstd(pl[:,k],ddof=1)),1) if nm>1 else 0.0,in_CTD="yes" if lo<=rn<=hi else ""))
        tables[c]=(rows,[f"m{m+1}" for m in range(nm)],pc["models"])
        with open(out/f"{name}_chain{c}_SS.csv","w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]);w.writeheader();w.writerows(rows)
        ctd=(pc["resnum"]>=lo)&(pc["resnum"]<=hi)
        for k in ("H","L","S"):pooled[k].append(pl[:,ctd][ss[:,ctd]==k])
        cons_all.append((c,pc["resnum"][ctd],cons[ctd],agree[ctd]/nm,np.nanmean(pl[:,ctd],0)))
    allpl=np.concatenate([np.array(per_chain[c]["plddt")[:,(per_chain[c]["resnum"]>=lo)&(per_chain[c]["resnum"]<=hi)].ravel() for c in per_chain])
    for k,lab in (("H","helix"),("L","loop"),("S","strand")):
        v=np.concatenate(pooled[k]) if pooled[k] else np.array([]);summ[f"CTD_pct_{lab}_calls"]="" if not v.size else round(100*v.size/max(allpl.size,1),1);summ[f"CTD_pLDDT_{lab}_mean"]="" if not v.size else round(float(np.nanmean(v)),1);summ[f"CTD_pLDDT_{lab}_sd"]="" if not v.size else round(float(np.nanstd(v)),1)
    cons_cat=np.concatenate([x[2] for x in cons_all]);agr_cat=np.concatenate([x[3] for x in cons_all])
    summ["CTD_consensus_helix_pct"]=round(100*float((cons_cat=="H").mean()),1);summ["CTD_strong_consensus_pct"]=round(100*float((agr_cat>=0.8).mean()),1);summ["CTD_no_majority_pct"]=round(100*float((agr_cat<=0.5).mean()),1)
    helix_ctd=np.concatenate(pooled["H"]) if pooled["H"] else np.array([]);summ["CTD_helix_calls_pLDDT_ge70_pct"]="" if not helix_ctd.size else round(100*float((helix_ctd>=70).mean()),1)
    (out/"summary.json").write_text(json.dumps(summ));log(f"  CTD: {summ['CTD_consensus_helix_pct']}% consensus helix",fh);return summ


def main():
    ap=argparse.ArgumentParser();ap.add_argument("--conditions",default=str(HERE/"conditions.csv"));ap.add_argument("--data",default=str(HERE/"data"));ap.add_argument("--results",default=str(HERE/"results"));ap.add_argument("--only");ap.add_argument("--force");ap.add_argument("--format",default="svg,pdf,png");a=ap.parse_args()
    conds=[r for r in csv.DictReader(open(a.conditions,newline="",encoding="utf-8-sig")) if r.get("name") and not r["name"].startswith("#")]
    out_root=Path(a.results);out_root.mkdir(parents=True,exist_ok=True);only=set(a.only.split(",")) if a.only else None;summaries=[]
    with open(out_root/"ss_run_log.txt","a",encoding="utf-8") as fh:
        for c in conds:
            if only and c["name"] not in only:continue
            try:
                s=run_condition(c,Path(a.data),out_root,a,fh)
                if s:summaries.append(s)
            except Exception:log(f"ERROR {c['name']}:\n{traceback.format_exc()}",fh)
    if summaries:
        fields=list(dict.fromkeys(k for s in summaries for k in s))
        with open(out_root/"all_conditions_SS_summary.csv","w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(summaries)


if __name__=="__main__":main()
