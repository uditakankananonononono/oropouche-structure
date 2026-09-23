#!/usr/bin/env python3
"""Outbreak-vs-historical substitution caller. Reference-anchored, manifest-reproducible.
Conventions: historical majority per reference column; conserved = majority frac >=0.9 with n>=20.
Substitution = outbreak residue differs from historical majority at a conserved column.
PC3 spike-in control runs FIRST; script refuses to emit real results if PC3 fails."""
import csv, re, sys, os
from collections import Counter, defaultdict
from Bio import SeqIO
import parasail
from parasail import blosum62

D = os.path.expanduser("~/oropouche/data"); R = os.path.expanduser("~/oropouche/results")
os.makedirs(R, exist_ok=True)

era = {}
for r in csv.DictReader(open(f"{D}/manifest.tsv"), delimiter="\t"):
    era[r["acc"]] = r["era"]

def cigar_map(ref, qry, cig):
    """ref col (0-based) -> qry col (0-based) or None. Auto-detects I/D convention."""
    ops = re.findall(r'(\d+)([=XID])', cig)
    def build(swap):
        m, rp, qp = {}, 0, 0
        for n, op in ops:
            n = int(n)
            if op in '=X':
                for i in range(n): m[rp+i] = qp+i
                rp += n; qp += n
            elif op == 'I':
                if swap: qp += n
                else: rp += n
            elif op == 'D':
                if swap: rp += n
                else: qp += n
        return m
    ma, mb = build(False), build(True)
    ca = sum(1 for p,q in ma.items() if q is not None and q < len(qry) and p < len(ref) and ref[p]==qry[q])
    cb = sum(1 for p,q in mb.items() if q is not None and q < len(qry) and p < len(ref) and ref[p]==qry[q])
    return ma if ca >= cb else mb

def align_map(ref, qry):
    res = parasail.nw_trace_striped_16(ref, qry, 10, 1, blosum62)
    c = res.cigar.decode
    c = c.decode() if isinstance(c, bytes) else c
    return cigar_map(ref, qry, c)

def load_prots(path):
    return {r.id: str(r.seq) for r in SeqIO.parse(path, "fasta")}

def pick_reference(seqs, hists):
    cand = [a for a in hists]
    lens = Counter(len(seqs[a]) for a in cand)
    modal = lens.most_common(1)[0][0]
    tied = sorted(a for a in cand if len(seqs[a]) == modal)
    return tied[0]

def hist_profile(ref, hseqs):
    prof = {}
    for a, s in hseqs.items():
        m = align_map(ref, s)
        for rp, qp in m.items():
            if qp is not None and qp < len(s):
                prof.setdefault(rp, Counter())[s[qp]] += 1
    maj = {}
    for rp, c in prof.items():
        n = sum(c.values()); aa, cnt = c.most_common(1)[0]
        maj[rp] = (aa, cnt/n, n)
    return maj

def call_outbreak(ref, maj, oseqs):
    calls = defaultdict(lambda: [0, 0])
    for a, s in oseqs.items():
        m = align_map(ref, s)
        for rp, qp in m.items():
            if rp not in maj: continue
            aa, frac, n = maj[rp]
            if frac < 0.9 or n < 20: continue
            if qp is None or qp >= len(s): continue
            q = s[qp]
            calls[(rp, q)][1] += 1
            if q != aa: calls[(rp, q)][0] += 1
    out = []
    for (rp, q), (cnt, cov) in calls.items():
        if cnt > 0: out.append((rp, q, cnt, cov))
    return out

# ---------- PC3 SPIKE-IN FIRST ----------
Lseqs = load_prots(f"{D}/prot_L.fasta")
Lh = {a: s for a, s in Lseqs.items() if era.get(a) == "historical"}
ref_acc = pick_reference(Lseqs, list(Lh)); ref = Lseqs[ref_acc]
maj = hist_profile(ref, Lh)
conserved = sorted(rp for rp,(aa,f,n) in maj.items() if f >= 0.95 and n >= 40)
spike_pos = [conserved[i] for i in (50, 500, 1000, 1500, 2000) if i < len(conserved)]
spike_pos = spike_pos[:5]
alt = lambda aa: 'A' if aa != 'A' else 'C'
mut = list(ref)
expect = []
for rp in spike_pos:
    aa, f, n = maj[rp]
    if ref[rp] != aa: continue
    mut[rp] = alt(aa); expect.append((rp, alt(aa)))
spiked = "".join(mut)
base = sorted((rp, q) for rp, q, c, cov in call_outbreak(ref, maj, {"BASE": ref}))
got = sorted((rp, q) for rp, q, c, cov in call_outbreak(ref, maj, {"SPIKE": spiked}))
exp = sorted(set(base) | set(expect))
if got != exp:
    print(f"PC3 FAILED: expected {exp} got {got}"); sys.exit(2)
print(f"PC3 PASSED: {len(exp)-len(base)}/{len(expect)} engineered substitutions recovered exactly; caller output on spiked sequence equals baseline({len(base)} ref-vs-majority calls) + engineered set, 0 false positives, 0 missed", flush=True)

# ---------- REAL RUN ----------
REGIONS_M = lambda p: "Gn" if p < 350 else ("NSm" if p < 481 else "Gc")
def l_domain(rp, refL):
    sdd = refL.find("SDD")
    if rp < 185: return "L-endonuclease"
    if sdd != -1 and abs(rp - sdd) < 400: return "L-polymerase-core"
    return "L-other"

sets = [("L", f"{D}/prot_L.fasta", None), ("M", f"{D}/prot_Mpoly.fasta", "KC759123.1"), ("N", f"{D}/prot_N.fasta", None)]
rows = []
for name, path, force_ref in sets:
    seqs = load_prots(path)
    h = {a: s for a, s in seqs.items() if era.get(a) == "historical"}
    o = {a: s for a, s in seqs.items() if era.get(a) == "outbreak"}
    ra = force_ref if force_ref in seqs else pick_reference(seqs, list(h))
    rf = seqs[ra]
    mj = hist_profile(rf, h)
    calls = call_outbreak(rf, mj, o)
    n_out = len(o)
    agg = defaultdict(lambda: [0, 0])
    for rp, q, c, cov in calls:
        pass
    # recompute coverage-aware fractions
    per = defaultdict(Counter); cover = Counter()
    for a, s in o.items():
        m = align_map(rf, s)
        for rp, (aa, f, n) in mj.items():
            if f < 0.9 or n < 20: continue
            qp = m.get(rp)
            if qp is None or qp >= len(s): continue
            cover[rp] += 1
            per[rp][s[qp]] += 1
    kept = 0
    for rp, c in per.items():
        aa, f, n = mj[rp]
        for alt_aa, cnt in c.items():
            if alt_aa == aa or cnt < 3: continue
            frac = cnt / cover[rp]
            if frac < 0.05: continue
            region = REGIONS_M(rp) if name == "M" else (l_domain(rp, rf) if name == "L" else "N")
            rows.append({"protein": name, "ref_pos": rp+1, "region": region, "hist_aa": aa,
                         "hist_frac": round(f, 3), "alt_aa": alt_aa, "out_count": cnt,
                         "out_cov": cover[rp], "out_frac": round(frac, 3), "reference": ra})
            kept += 1
    print(f"{name}: ref={ra} hist={len(h)} outbreak={len(o)} conserved_cols={sum(1 for v in mj.values() if v[1]>=0.9 and v[2]>=20)} subs_kept={kept}", flush=True)

with open(f"{R}/substitutions.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(rows)
print("wrote", f"{R}/substitutions.tsv", len(rows), "rows")
