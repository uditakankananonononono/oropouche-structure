#!/usr/bin/env python3
"""Fetch + byte-lock OROV L/M/S segment sequences from NCBI nuccore."""
import json, time, hashlib, os, re, sys
from urllib.request import urlopen
from urllib.parse import urlencode

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
OUT = os.path.expanduser("~/oropouche/data")
os.makedirs(OUT, exist_ok=True)

def get(url):
    last = None
    for a in range(5):
        try:
            return urlopen(url, timeout=60).read()
        except Exception as e:
            last = e; time.sleep(2*(a+1))
    raise last

def esearch(term):
    url = BASE+"esearch.fcgi?"+urlencode({"db":"nuccore","term":term,"retmax":5000,"retmode":"json"})
    return json.loads(get(url))["esearchresult"]["idlist"]

def esummary(ids):
    out = []
    for i in range(0, len(ids), 200):
        url = BASE+"esummary.fcgi?"+urlencode({"db":"nuccore","id":",".join(ids[i:i+200]),"retmode":"json","version":"2.0"})
        j = json.loads(get(url))
        for uid in j["result"]["uids"]:
            out.append(j["result"][uid])
        time.sleep(0.35)
    return out

def efetch_fasta(ids):
    recs = []
    for i in range(0, len(ids), 150):
        url = BASE+"efetch.fcgi?"+urlencode({"db":"nuccore","id":",".join(ids[i:i+150]),"rettype":"fasta","retmode":"text"})
        txt = get(url).decode()
        hdr, seq = None, []
        for line in txt.splitlines():
            if line.startswith(">"):
                if hdr: recs.append((hdr, "".join(seq)))
                hdr, seq = line[1:], []
            else: seq.append(line.strip())
        if hdr: recs.append((hdr, "".join(seq)))
        time.sleep(0.35)
    return recs

QUERIES = {
 "L": 'Oropouche virus[Organism] AND 6000:8000[SLEN]',
 "M": 'Oropouche virus[Organism] AND 3800:5200[SLEN]',
 "S": 'Oropouche virus[Organism] AND 650:1300[SLEN]',
}
def seg_ok(seg, title):
    t = title.lower()
    if seg=="L": return ("segment l" in t or "large segment" in t or "l segment" in t or "rna polymerase" in t or "polymerase (l)" in t) and "segment m" not in t
    if seg=="M": return ("segment m" in t or "medium segment" in t or "m segment" in t or "glycoprotein" in t or "polyprotein" in t)
    if seg=="S": return ("segment s" in t or "small segment" in t or "s segment" in t or "nucleocapsid" in t or "nucleoprotein" in t)
    return False

def year_of(cdate):
    if not cdate: return None
    m = re.findall(r'(19|20)\d{2}', cdate)
    return int(re.search(r'((?:19|20)\d{2})', cdate).group(1)) if re.search(r'((?:19|20)\d{2})', cdate) else None

manifest_rows = []
for seg, q in QUERIES.items():
    ids = esearch(q)
    print(f"{seg}: esearch {len(ids)} ids", flush=True)
    summ = esummary(ids)
    keep, n_nocomplete, n_wrongseg = [], 0, 0
    for r in summ:
        title = r.get("title","")
        subtypes = (r.get("subtype","") or "").split("|")
        subnames = (r.get("subname","") or "").split("|")
        meta = dict(zip(subtypes, subnames))
        complete = ("complete" in title.lower())
        if not seg_ok(seg, title): n_wrongseg += 1; continue
        if not complete: n_nocomplete += 1; continue
        keep.append({"uid": r["uid"], "acc": r.get("accessionversion", r.get("caption","")),
                     "title": title, "slen": r.get("slen",0),
                     "cdate": meta.get("collection_date",""), "country": meta.get("country",""),
                     "host": meta.get("host","")})
    print(f"{seg}: kept {len(keep)} complete+segment-matched (dropped {n_nocomplete} not-complete, {n_wrongseg} wrong-seg)", flush=True)
    recs = efetch_fasta([k["uid"] for k in keep])
    by_uid = {k["uid"]: k for k in keep}
    fasta_path = os.path.join(OUT, f"orov_{seg}_complete.fasta")
    with open(fasta_path, "w") as fh:
        for hdr, seq in recs:
            acc = hdr.split()[0]
            k = next((v for v in keep if v["acc"]==acc), None)
            if k is None:
                k = next((v for v in keep if v["acc"].split(".")[0]==acc.split(".")[0]), None)
            if k is None: continue
            fh.write(f">{acc}\n")
            for i in range(0, len(seq), 70): fh.write(seq[i:i+70]+"\n")
            yr = year_of(k["cdate"])
            era = "outbreak" if yr and yr>=2022 else ("historical" if yr and yr<2020 else "buffer/unknown")
            rec_text = f">{acc}\n{seq}\n"
            manifest_rows.append({"acc":acc,"segment":seg,"slen":len(seq),"year":yr or "",
                "era":era,"cdate":k["cdate"],"country":k["country"],"host":k["host"],
                "sha256":hashlib.sha256(rec_text.encode()).hexdigest(),"title":k["title"]})
    print(f"{seg}: wrote {len(recs)} fasta records -> {fasta_path}", flush=True)

import csv
mpath = os.path.join(OUT, "manifest.tsv")
with open(mpath,"w",newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(manifest_rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(manifest_rows)
from collections import Counter
c = Counter((r["segment"], r["era"]) for r in manifest_rows)
print("== MANIFEST COUNTS ==")
for k in sorted(c): print(k, c[k])
print("manifest:", mpath)
