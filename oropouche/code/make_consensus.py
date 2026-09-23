#!/usr/bin/env python3
import csv, os, re
from collections import Counter
from Bio import SeqIO
import parasail
from parasail import blosum62
D=os.path.expanduser("~/oropouche/data")
era={r["acc"]:r["era"] for r in csv.DictReader(open(f"{D}/manifest.tsv"),delimiter="\t")}
def cigar_map(ref,qry,cig):
    ops=re.findall(r'(\d+)([=XID])',cig)
    def build(swap):
        m,rp,qp={},0,0
        for n,op in ops:
            n=int(n)
            if op in '=X':
                for i in range(n): m[rp+i]=qp+i
                rp+=n; qp+=n
            elif op=='I':
                if swap: qp+=n
                else: rp+=n
            else:
                if swap: rp+=n
                else: qp+=n
        return m
    ma,mb=build(False),build(True)
    ca=sum(1 for p,q in ma.items() if q<len(qry) and p<len(ref) and ref[p]==qry[q])
    cb=sum(1 for p,q in mb.items() if q<len(qry) and p<len(ref) and ref[p]==qry[q])
    return ma if ca>=cb else mb
def align_map(ref,qry):
    res=parasail.nw_trace_striped_16(ref,qry,10,1,blosum62)
    c=res.cigar.decode
    c=c.decode() if isinstance(c,bytes) else c
    return cigar_map(ref,qry,c)
def outbreak_consensus(ref_acc, seqs):
    ref=seqs[ref_acc]
    prof={}
    for a,s in seqs.items():
        if era.get(a)!="outbreak": continue
        m=align_map(ref,s)
        for rp,qp in m.items():
            if qp is not None and qp<len(s): prof.setdefault(rp,Counter())[s[qp]]+=1
    cons=list(ref)
    for rp,c in prof.items():
        if sum(c.values())>=10: cons[rp]=c.most_common(1)[0][0]
    return "".join(cons)
L={r.id:str(r.seq) for r in SeqIO.parse(f"{D}/prot_L.fasta","fasta")}
M={r.id:str(r.seq) for r in SeqIO.parse(f"{D}/prot_Mpoly.fasta","fasta")}
N={r.id:str(r.seq) for r in SeqIO.parse(f"{D}/prot_N.fasta","fasta")}
consL=outbreak_consensus("KP026179.1",L)
consM=outbreak_consensus("KC759123.1",M)
consN=outbreak_consensus("HM470107.1",N)
targets={
 "N_full": consN,
 "Gn": consM[:350],
 "Gc_head": consM[481:703],
 "L_endo": consL[:185],
}
h3=str(list(SeqIO.parse(f"{D}/pdb_6H3X.fasta","fasta"))[0].seq)
targets["Gc_head_6H3Xseq"]=h3
with open(f"{D}/structure_targets.fasta","w") as fh:
    for k,v in targets.items():
        fh.write(f">{k}\n")
        for i in range(0,len(v),70): fh.write(v[i:i+70]+"\n")
for k,v in targets.items(): print(k,len(v),"aa")
