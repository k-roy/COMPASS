#!/usr/bin/env python3
"""Generate a COMPASS-format introns TSV (chrom start end strand type Name intron_type)
from a gencode GTF, matching the yeast convention: intron BODY in 1-based inclusive
coords = [prev_exon.end+1, next_exon.start-1]. Dedups identical introns across transcripts."""
import re, sys
from collections import defaultdict

gtf, out = sys.argv[1], sys.argv[2]
tx_exons = defaultdict(list)
tx_info = {}
for line in open(gtf):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] != "exon":
        continue
    m = re.search(r'transcript_id "([^"]+)"', f[8])
    if not m:
        continue
    tid = m.group(1)
    g = re.search(r'gene_name "([^"]+)"', f[8])
    tx_exons[tid].append((int(f[3]), int(f[4])))
    tx_info[tid] = (f[0], f[6], g.group(1) if g else tid)

introns = {}
for tid, exons in tx_exons.items():
    chrom, strand, gname = tx_info[tid]
    ex = sorted(exons)
    for i in range(len(ex) - 1):
        istart, iend = ex[i][1] + 1, ex[i + 1][0] - 1
        if iend < istart:
            continue
        introns.setdefault((chrom, istart, iend, strand), gname + "_intron")

with open(out, "w") as fh:
    fh.write("chrom\tstart\tend\tstrand\ttype\tName\tintron_type\n")
    for (chrom, s, e, strand), name in sorted(introns.items()):
        fh.write(f"{chrom}\t{s}\t{e}\t{strand}\tintron\t{name}\tspliceosomal_intron\n")
print(f"transcripts={len(tx_exons)} unique_introns={len(introns)}", file=sys.stderr)
