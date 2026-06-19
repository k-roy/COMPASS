#!/usr/bin/env python3
"""Dependency-light replacement for jvarkit samfixcigar (COMPASS needs SAM-1.4 =/X
ops so its edit-distance scoring can count mismatches from the CIGAR). Rewrites every
M (op 0) run to runs of = (op 7, match) / X (op 8, mismatch) by comparing read bases
to the reference; I/D/N/S/H are preserved. Avoids the jvarkit JDK-17/gradle build.

Usage: samfixcigar.py in.bam out.bam reference.fa
Equivalent to: java -jar samfixcigar --reference reference.fa --samoutputformat BAM in.bam

Robustness: the reference name is whitespace-trimmed (BBMap keeps the full FASTA
header e.g. 'chr1 1' while the .fai/other aligners use 'chr1'), and each read is
processed in a try/except so one problematic read writes through unchanged instead
of truncating the whole BAM.
"""
import sys
import pysam

inbam, outbam, ref = sys.argv[1], sys.argv[2], sys.argv[3]
fa = pysam.FastaFile(ref)
src = pysam.AlignmentFile(inbam, "rb")
dst = pysam.AlignmentFile(outbam, "wb", template=src)

n_fixed = n_passthru = 0
for r in src:
    if r.is_unmapped or not r.cigartuples:
        dst.write(r)
        continue
    try:
        seq = r.query_sequence or ""
        chrom = (r.reference_name or "").split()[0]   # BBMap keeps full header 'chr1 1' -> 'chr1'
        new = []
        qpos, rpos = 0, r.reference_start
        for op, ln in r.cigartuples:
            if op in (0, 7, 8):                       # M/=/X -> recompute =/X vs reference
                refseq = fa.fetch(chrom, rpos, rpos + ln).upper()
                run_op, run_len = None, 0
                for i in range(ln):
                    rb = seq[qpos + i].upper() if qpos + i < len(seq) else "N"
                    cb = refseq[i] if i < len(refseq) else "N"
                    o = 7 if (rb == cb and rb != "N") else 8
                    if o == run_op:
                        run_len += 1
                    else:
                        if run_op is not None:
                            new.append((run_op, run_len))
                        run_op, run_len = o, 1
                if run_op is not None:
                    new.append((run_op, run_len))
                qpos += ln
                rpos += ln
            elif op == 1:                             # I (consumes query)
                new.append((op, ln)); qpos += ln
            elif op in (2, 3):                        # D, N (consume ref)
                new.append((op, ln)); rpos += ln
            elif op == 4:                             # S (consumes query)
                new.append((op, ln)); qpos += ln
            else:                                     # H (5), P (6)
                new.append((op, ln))
        r.cigartuples = new
        dst.write(r)
        n_fixed += 1
    except Exception as exc:
        sys.stderr.write(f"samfixcigar: read {r.query_name} passthrough ({exc})\n")
        dst.write(r)
        n_passthru += 1

dst.close(); src.close(); fa.close()
sys.stderr.write(f"samfixcigar: fixed={n_fixed} passthrough={n_passthru}\n")
