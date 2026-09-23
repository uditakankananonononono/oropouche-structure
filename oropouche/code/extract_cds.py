#!/usr/bin/env python3
"""Extract CDS translations from GenBank records for byte-locked accessions."""
import time, os, sys, csv, json
from urllib.request import urlopen
from urllib.parse import urlencode
from io import StringIO
from Bio import SeqIO

BASE="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
OUT=os.path.expanduser("~/oropouche/data")
def get(url):
    last=None
    for a in range(5):
        try: return urlopen(url, timeout=90).read()
        except Exception as e: last=e; time.sleep(2*(a+1))
    raise last

PRODUCTS = {
 "L":  ["polymerase","rdrp","l protein","rna-directed rna polymerase"],
 "N":  ["nucleoprotein","nucleocapsid"],
 "NSs":["nonstructural s","nss","non-structural s"],
 "M":  ["polyprotein","envelopment polyprotein","glycoprotein precursor","m polyprotein"],
}
def classify(seg, product):
    p=(product or "").lower()
    if seg=="L":
        return "L" if any(k in p for k in PRODUCTS["L"]) else None
    if seg=="S":
        if any(k in p for k in PRODUCTS["NSs"]): return "NSs"
        if any(k in p for k in PRODUCTS["N"]): return "N"
        return None
    if seg=="M":
        return "Mpoly" if any(k in p for k in PRODUCTS["M"]) else None
    return None

rows=list(csv.DictReader(open(f"{OUT}/manifest.tsv"),delimiter="\t"))
by_seg={}
for r in rows: by_seg.setdefault(r["segment"],[]).append(r["acc"])
stats={}
seqs={}
for seg, accs in by_seg.items():
    stats[seg]={"records":len(accs),"with_target_cds":0,"no_gb":0,"no_match":0}
    for i in range(0,len(accs),100):
        chunk=accs[i:i+100]
        url=BASE+"efetch.fcgi?"+urlencode({"db":"nuccore","id":",".join(chunk),"rettype":"gb","retmode":"text"})
        txt=get(url).decode()
        for rec in SeqIO.parse(StringIO(txt),"genbank"):
            acc=rec.id
            found=False
            for f in rec.features:
                if f.type!="CDS": continue
                prod=f.qualifiers.get("product",[""])[0]
                kind=classify(seg,prod)
                if not kind: continue
                if kind=="Mpoly" and len(f.extract(rec.seq))<3900: continue  # want full polyprotein
                transl=f.qualifiers.get("translation",[None])[0]
                if not transl:
                    try: transl=str(f.extract(rec.seq).translate(to_stop=False)).rstrip("*")
                    except Exception: continue
                key=(kind,acc)
                if key not in seqs or len(transl)>len(seqs[key][1]):
                    seqs[key]=(prod,transl)
                found=True
            if found: stats[seg]["with_target_cds"]+=1
            else: stats[seg]["no_match"]+=1
        time.sleep(0.35)
    print(seg, stats[seg], flush=True)

for kind in ["L","Mpoly","N","NSs"]:
    n=0
    with open(f"{OUT}/prot_{kind}.fasta","w") as fh:
        for (k,acc),(prod,tr) in sorted(seqs.items()):
            if k!=kind: continue
            fh.write(f">{acc}\n")
            for i in range(0,len(tr),70): fh.write(tr[i:i+70]+"\n")
            n+=1
    print(kind, n, "proteins")
