#!/usr/bin/env python3
"""
make_fuzz_data.py -- runtime files and seed corpora for the samtools fuzzers
(view, consensus, calmd, fixmate, markdup, sort, flagstat, stats,
depth, mpileup, fastq, index, coverage, bedcov, regions, merge, cat,
collate, split, addreplacerg, reheader, ampliconclip, ampliconstats,
dict, import, reference, reset, depad, phase, targetcut, cram_size,
samples, checksum, head, quickcheck, rmdup, tview, faidx and faidx_idx).

Usage: ./make_fuzz_data.py [OUTDIR]          (default: fuzz_data)

Pure Python 3 with no dependencies, and deterministic.  If samtools is on
PATH, BAM and CRAM (no_ref) copies of the SAM seeds are produced as well.

Output:
  OUTDIR/runtime/        copy into the fuzzer's working directory
  OUTDIR/sam/, fasta/    the seed payloads without control bytes, to read
  OUTDIR/seeds/<fuzzer>/ initial corpora, control bytes prepended
"""

import gzip
import os
import random
import re
import shutil
import struct
import subprocess
import sys
import zlib

OUT = sys.argv[1] if len(sys.argv) > 1 else "fuzz_data"
rng = random.Random(12345)

BASES = "ACGT"
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def rand_seq(n):
    return "".join(rng.choice(BASES) for _ in range(n))


def rand_qual(n, lo=2, hi=40):
    return "".join(chr(33 + rng.randint(lo, hi)) for _ in range(n))


def write(path, data):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if isinstance(data, str):
        data = data.encode()
    with open(path, "wb") as f:
        f.write(data)


# ---------------------------------------------------------------------------
# FASTA / FASTQ, .fai, BGZF and .gzi

def fasta(seqs, width=60, eol="\n", desc=None):
    out = []
    for name, s in seqs:
        d = (" " + desc[name]) if desc and name in desc else ""
        out.append(">" + name + d + eol)
        if width:
            for i in range(0, len(s), width):
                out.append(s[i:i + width] + eol)
        else:
            out.append(s + eol)
    return "".join(out)


def fastq(recs, width=0):
    out = []
    for name, s, q in recs:
        out.append("@" + name + "\n")
        w = width or max(len(s), 1)
        out += [s[i:i + w] + "\n" for i in range(0, len(s), w)]
        out.append("+\n")
        out += [q[i:i + w] + "\n" for i in range(0, len(q), w)]
    return "".join(out)


def make_fai(data):
    """Index a well-formed FASTA or FASTQ (bytes), as fai_build would."""
    lines, off = [], 0
    for raw in data.splitlines(keepends=True):
        lines.append((off, raw.rstrip(b"\r\n"), len(raw)))
        off += len(raw)
    out, i = [], 0
    while i < len(lines):
        content = lines[i][1]
        if content.startswith(b">"):
            name = content[1:].split()[0].decode()
            i += 1
            seq_off, blen, wlen = lines[i][0], len(lines[i][1]), lines[i][2]
            length = 0
            while i < len(lines) and not lines[i][1].startswith(b">"):
                length += len(lines[i][1])
                i += 1
            out.append(f"{name}\t{length}\t{seq_off}\t{blen}\t{wlen}\n")
        elif content.startswith(b"@"):
            name = content[1:].split()[0].decode()
            i += 1
            seq_off, blen, wlen = lines[i][0], len(lines[i][1]), lines[i][2]
            length = 0
            while not lines[i][1].startswith(b"+"):
                length += len(lines[i][1])
                i += 1
            i += 1
            qual_off, q = lines[i][0], 0
            while q < length:
                q += len(lines[i][1])
                i += 1
            out.append(f"{name}\t{length}\t{seq_off}\t{blen}\t{wlen}\t{qual_off}\n")
        else:
            i += 1
    return "".join(out).encode()


BGZF_EOF = bytes.fromhex(
    "1f8b08040000000000ff0600424302001b0003000000000000000000")


def bgzf(data, block=65280):
    """BGZF-compress data; returns (compressed, gzi).  Small block sizes give
    multi-block files, so the .gzi has entries and random access is real."""
    out, offsets = bytearray(), []
    for i in range(0, len(data), block):
        chunk = data[i:i + block]
        offsets.append((len(out), i))
        c = zlib.compressobj(6, zlib.DEFLATED, -15)
        comp = c.compress(chunk) + c.flush()
        out += struct.pack("<BBBBIBBHBBHH", 31, 139, 8, 4, 0, 0, 255,
                           6, 66, 67, 2, len(comp) + 25)
        out += comp
        out += struct.pack("<II", zlib.crc32(chunk) & 0xffffffff, len(chunk))
    out += BGZF_EOF
    assert gzip.decompress(bytes(out)) == data
    entries = offsets[1:]  # the first block (0,0) is implicit
    gzi = struct.pack("<Q", len(entries)) + b"".join(
        struct.pack("<QQ", c, u) for c, u in entries)
    return bytes(out), gzi


# ---------------------------------------------------------------------------
# Reference

def make_reference():
    ref = {}
    s = list(rand_seq(900))
    s[300:310] = "A" * 10                       # homopolymer
    s[500:520] = "CA" * 10                      # dinucleotide repeat
    ref["chr1"] = "".join(s)
    s = list(rand_seq(600))
    s[50:56] = "N" * 6                          # N run
    s[200:240] = [c.lower() for c in s[200:240]]  # soft-masked
    ref["chr2"] = "".join(s)
    ref["chr3:1-10"] = rand_seq(200)            # name that parses as a region
    ref["chrM"] = rand_seq(120)
    return ref


REF = make_reference()
TIDS = {name: i for i, name in enumerate(REF)}


# ---------------------------------------------------------------------------
# SAM records

def parse_cigar(cigar):
    return [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)]


def qlen(cigar):
    return sum(n for n, op in parse_cigar(cigar) if op in "MIS=X")


def refend(pos0, cigar):
    return pos0 + sum(n for n, op in parse_cigar(cigar) if op in "MDN=X")


def simulate(refseq, pos0, cigar, err=0.0, force_mm=()):
    """Read sequence for an alignment at 0-based pos0 with this CIGAR."""
    seq, p, q = [], pos0, 0
    for n, op in parse_cigar(cigar):
        if op in "M=X":
            for _ in range(n):
                r = refseq[p].upper() if p < len(refseq) else "N"
                b = r if r in BASES else rng.choice(BASES)
                if op == "X" or (op == "M" and (rng.random() < err or q in force_mm)):
                    b = rng.choice([x for x in BASES if x != b])
                seq.append(b)
                p += 1
                q += 1
        elif op in "IS":
            seq.append(rand_seq(n))
            q += n
        elif op in "DN":
            p += n
    return "".join(seq)


def md_nm(refseq, pos0, cigar, seq):
    """MD and NM as calmd computes them (N never matches)."""
    md, run, nm, p, q = [], 0, 0, pos0, 0
    for n, op in parse_cigar(cigar):
        if op in "M=X":
            for _ in range(n):
                r = refseq[p].upper() if p < len(refseq) else "N"
                if seq[q].upper() == r and r != "N":
                    run += 1
                else:
                    md.append(f"{run}{r}")
                    run = 0
                    nm += 1
                p += 1
                q += 1
        elif op in "IS":
            q += n
            if op == "I":
                nm += n
        elif op == "D":
            md.append(f"{run}^{refseq[p:p + n].upper()}")
            run = 0
            nm += n
            p += n
        elif op == "N":
            p += n
    md.append(str(run))
    return "".join(md), nm


class Rec:
    def __init__(self, qname, flag, rname, pos, mapq, cigar, seq, qual,
                 tags=(), rnext="*", pnext=0, tlen=0):
        self.qname, self.flag, self.rname, self.pos = qname, flag, rname, pos
        self.mapq, self.cigar, self.seq, self.qual = mapq, cigar, seq, qual
        self.tags, self.rnext, self.pnext, self.tlen = list(tags), rnext, pnext, tlen

    def key(self):
        return (TIDS.get(self.rname, 1 << 30), self.pos)

    def line(self):
        f = [self.qname, self.flag, self.rname, self.pos, self.mapq, self.cigar,
             self.rnext, self.pnext, self.tlen, self.seq, self.qual] + self.tags
        return "\t".join(map(str, f))


def aligned(qname, flag, rname, pos, cigar, mapq=60, err=0.0, md=True,
            nm=True, rg=None, tags=(), force_mm=(), seq=None):
    refseq, pos0 = REF[rname], pos - 1
    if seq is None:
        seq = simulate(refseq, pos0, cigar, err, force_mm)
    t = []
    if seq != "*" and (md or nm):
        m, n = md_nm(refseq, pos0, cigar, seq)
        if md:
            t.append(f"MD:Z:{m}")
        if nm:
            t.append(f"NM:i:{n}")
    if rg:
        t.append(f"RG:Z:{rg}")
    qual = "*" if seq == "*" else rand_qual(len(seq))
    return Rec(qname, flag, rname, pos, mapq, cigar, seq, qual, t + list(tags))


def header(so="coordinate", go=None):
    h = ["@HD\tVN:1.6\tSO:" + so + (f"\tGO:{go}" if go else "")]
    # No M5 tags, so CRAM decoding never tries to fetch a reference
    h += [f"@SQ\tSN:{n}\tLN:{len(s)}" for n, s in REF.items()]
    h += ["@RG\tID:rg1\tSM:s1\tLB:L\tPL:ILLUMINA",
          "@RG\tID:rg2\tSM:s2\tLB:L2\tPL:ONT",
          "@PG\tID:bwa\tPN:bwa\tVN:0.7.17"]
    return "\n".join(h) + "\n"


def sam(hdr, recs):
    for r in recs:
        if r.seq != "*" and r.cigar != "*":
            assert qlen(r.cigar) == len(r.seq), (r.qname, r.cigar)
        if r.qual != "*":
            assert len(r.qual) == len(r.seq), r.qname
    return hdr + "".join(r.line() + "\n" for r in recs)


def by_pos(recs):
    return sorted(recs, key=Rec.key)


# ---- Single-ended reads exercising CIGAR, MD/NM and flag variety ----------

def mixed_reads():
    rg = lambda i: ("rg1", "rg2", None)[i % 3]
    spec = [
        ("r01", 0,    "chr1", 10,  "40M",        {}),
        ("r02", 0,    "chr1", 50,  "40M",        {"err": 0.08}),
        ("r03", 16,   "chr1", 100, "5S30M5S",    {}),
        ("r04", 0,    "chr1", 120, "3H30M",      {}),
        ("r05", 16,   "chr1", 150, "15M2I15M",   {}),
        ("r06", 0,    "chr1", 180, "15M3D15M",   {"force_mm": {15}}),
        ("r07", 0,    "chr1", 200, "10M100N20M", {}),
        ("r08", 16,   "chr1", 290, "10=1X10=",   {}),   # homopolymer
        ("r09", 0,    "chr1", 330, "10M1P1I10M", {}),
        ("r10", 0,    "chr1", 60,  "1I20M3D5M",  {}),
        ("r11", 0,    "chr1", 495, "4S30M",      {}),   # CA repeat
        ("r12", 2048, "chr1", 600, "5H25M",      {"tags": ["SA:Z:chr1,10,+,25M5S,60,0;"]}),
        ("r13", 256,  "chr1", 400, "30M",        {"seq": "*", "md": False, "nm": False}),
        ("r14", 0,    "chr1", 700, "30M",        {"mapq": 0}),
        ("r15", 0,    "chr1", 750, "30M",        {"mapq": 255, "md": False, "nm": False}),
        ("r17", 0,    "chr1", 850, "30M",        {"nm": False}),
        ("r18", 0,    "chr1", 820, "30M",        {"md": False}),
        ("r19", 0,    "chr1", 880, "40M",        {}),   # runs off the end
        ("r20", 16,   "chr2", 40,  "40M",        {}),   # over the N run
        ("r21", 0,    "chr2", 195, "50M",        {"err": 0.05}),  # lowercase
        ("r22", 0,    "chr2", 570, "50M",        {}),   # runs off the end
        ("r23", 0,    "chr3:1-10", 5, "30M",     {}),
        ("r24", 1024, "chrM", 1,   "30M",        {}),
        ("r25", 512,  "chrM", 50,  "30M",        {}),
        ("r26", 0,    "chr1", 30,  "25M",        {"tags": [
            "XA:A:x", "XI:i:-5", "XF:f:1.5", "XZ:Z:hello world", "BX:Z:ACGT-1",
            "XB:B:c,1,-2,3", "XS:B:S,1000,2000", "XH:H:1AE301",
            "XC:i:4294967295"]}),
    ]
    recs = []
    for i, (q, flag, rname, pos, cigar, kw) in enumerate(spec):
        recs.append(aligned(q, flag, rname, pos, cigar, rg=rg(i), **kw))

    # Wrong MD/NM, which calmd should notice and fix
    r = aligned("r16", 0, "chr1", 780, "30M", md=False, nm=False)
    r.tags = ["MD:Z:10A5^C3", "NM:i:99"]
    recs.append(r)

    # N bases in the read
    r = aligned("r27", 0, "chrM", 80, "20M")
    r.seq = r.seq[:5] + "NNN" + r.seq[8:]
    r.tags = [t for t in r.tags if t[:2] not in ("MD", "NM")]
    m, n = md_nm(REF["chrM"], 79, "20M", r.seq)
    r.tags += [f"MD:Z:{m}", f"NM:i:{n}"]
    recs.append(r)

    # Placed unmapped, and unplaced unmapped (which sort last)
    recs.append(Rec("r28", 4, "chr1", 700, 0, "*", rand_seq(30), rand_qual(30)))
    recs.append(Rec("r29", 4, "*", 0, 0, "*", rand_seq(30), rand_qual(30)))
    recs.append(Rec("r30", 4, "*", 0, 0, "*", rand_seq(25), "*"))
    return by_pos(recs)


# ---- Pairs ---------------------------------------------------------------

def make_pair(name, r1, r2, orient="FR", flags=0, mate_unmapped=False):
    """r1, r2 = (rname, pos, cigar).  Returns two records with consistent
    mate information, TLEN and MC tags."""
    s1, s2 = {"FR": (0, 1), "RF": (1, 0), "FF": (0, 0)}[orient]
    a = aligned(name, 0, *r1, rg="rg1")
    if mate_unmapped:
        b = Rec(name, 0, r1[0], r1[1], 0, "*", rand_seq(30), rand_qual(30),
                ["RG:Z:rg1"])
    else:
        b = aligned(name, 0, *r2, rg="rg1")
    a.flag = 0x1 | 0x40 | (0x10 if s1 else 0) | flags
    b.flag = 0x1 | 0x80 | (0x10 if s2 else 0) | flags
    if mate_unmapped:
        a.flag |= 0x8
        b.flag |= 0x4
        b.flag &= ~0x10
    else:
        a.flag |= 0x20 if s2 else 0
        b.flag |= 0x20 if s1 else 0
        if orient == "FR" and r1[0] == r2[0]:
            a.flag |= 0x2
            b.flag |= 0x2
    for x, y in ((a, b), (b, a)):
        x.rnext = "=" if x.rname == y.rname else y.rname
        x.pnext = y.pos
        if y.cigar != "*":
            x.tags.append(f"MC:Z:{y.cigar}")
    if a.rname == b.rname and not mate_unmapped:
        s = min(a.pos, b.pos) - 1
        e = max(refend(a.pos - 1, a.cigar), refend(b.pos - 1, b.cigar))
        left = a if a.pos <= b.pos else b
        left.tlen = e - s
        (b if left is a else a).tlen = -(e - s)
    return [a, b]


def pair_reads():
    recs = []
    recs += make_pair("p01", ("chr1", 100, "30M"), ("chr1", 250, "30M"))
    recs += make_pair("p02", ("chr1", 300, "40M"), ("chr1", 310, "40M"))  # overlap
    recs += make_pair("p03", ("chr1", 380, "30M"), ("chr1", 420, "30M"), "RF")
    recs += make_pair("p04", ("chr2", 100, "30M"), ("chr2", 200, "30M"), "FF")
    recs += make_pair("p05", ("chr1", 500, "30M"), ("chr2", 300, "30M"))
    recs += make_pair("p06", ("chr1", 600, "30M"), None, mate_unmapped=True)
    recs += make_pair("p08", ("chr1", 700, "20M2D20M"), ("chr1", 720, "5S35M"))
    recs += make_pair("p09", ("chrM", 10, "30M"), ("chrM", 60, "30M"), flags=0x400)
    recs += make_pair("p10", ("chr2", 500, "30M"), ("chr2", 560, "50M"))
    recs.append(Rec("p07", 77, "*", 0, 0, "*", rand_seq(30), rand_qual(30)))
    recs.append(Rec("p07", 141, "*", 0, 0, "*", rand_seq(30), rand_qual(30)))
    return recs


def collated_reads():
    """Name-grouped pairs with stale or missing mate information, for
    fixmate to fill in."""
    recs = pair_reads()
    for i, r in enumerate(recs):
        r.tags = [t for t in r.tags if not t.startswith("MC:")]
        if i % 3 == 0:
            r.rnext, r.pnext, r.tlen = "*", 0, 0
            r.flag &= ~(0x2 | 0x20)
        elif i % 3 == 1:
            r.pnext, r.tlen = r.pnext + 7, 12345      # wrong
    # Secondary and supplementary in their groups, a singleton and an
    # unpaired read, and a name with two read1s
    sec = aligned("p01", 0x1 | 0x40 | 0x100, "chr2", 400, "30M", mapq=3)
    sup = aligned("p02", 0x1 | 0x80 | 0x800 | 0x10, "chr1", 800, "10H30M",
                  tags=["SA:Z:chr1,310,-,40M,60,0;"])
    recs += [sec, sup]
    recs.append(aligned("s01", 0x1 | 0x40, "chr1", 140, "30M"))
    recs.append(aligned("u01", 0, "chr2", 450, "30M"))
    recs.append(aligned("d01", 0x1 | 0x40, "chr1", 160, "30M"))
    recs.append(aligned("d01", 0x1 | 0x40 | 0x10, "chr1", 260, "30M"))
    order = {}
    for r in recs:
        order.setdefault(r.qname, len(order))
    return sorted(recs, key=lambda r: (order[r.qname], r.flag & 0x900, r.flag & 0xc0))


# ---- Base modifications --------------------------------------------------

def mm_tag(orig, mods=(("C", "m", ""),), frac=0.5):
    """MM/ML for a sequence in original read orientation."""
    mm, ml = [], []
    for base, code, mode in mods:
        idx = [i for i, b in enumerate(orig) if b == base]
        chosen = [k for k in range(len(idx)) if rng.random() < frac] or [0]
        if not idx:
            chosen = []
        deltas, prev = [], -1
        for k in chosen:
            deltas.append(k - prev - 1)
            prev = k
        mm.append(f"{base}+{code}{mode}" + "".join(f",{d}" for d in deltas) + ";")
        ml += [rng.randint(0, 255) for _ in chosen]
    tags = ["MM:Z:" + "".join(mm)]
    if ml:
        tags.append("ML:B:C," + ",".join(map(str, ml)))
    return tags


def basemod_reads():
    recs = []

    # Pair with mods; hard-clipped secondary/supplementary copying the
    # primary's tags without MN (the case fixmate -M repairs), and one with MN
    a, b = make_pair("bm01", ("chr1", 100, "40M"), ("chr1", 200, "40M"))
    a.tags += mm_tag(a.seq)
    b.tags += mm_tag(revcomp(b.seq))
    sec = Rec("bm01", 0x1 | 0x40 | 0x100, "chr2", 300, 0, "10H30M",
              a.seq[10:], a.qual[10:], [t for t in a.tags if t[:2] in ("MM", "ML")])
    sec2 = Rec("bm01", 0x1 | 0x40 | 0x100, "chr2", 350, 0, "30M10H",
               a.seq[:30], a.qual[:30],
               [t for t in a.tags if t[:2] in ("MM", "ML")] + ["MN:i:40"])
    sup = Rec("bm01", 0x1 | 0x80 | 0x800 | 0x10, "chr1", 600, 0, "30M10H",
              b.seq[:30], b.qual[:30], [t for t in b.tags if t[:2] in ("MM", "ML")])
    for r in (a, b):
        r.tags = [t for t in r.tags if not t.startswith("MC:")]
    recs += [a, sec, sec2, b, sup]

    # Multiple modification types, implicit/explicit modes
    r = aligned("bm02", 0, "chr1", 400, "40M")
    r.tags += mm_tag(r.seq, (("C", "m", "?"), ("C", "h", "."), ("A", "a", "")))
    recs.append(r)

    # Reverse strand single read
    r = aligned("bm03", 16, "chr2", 100, "40M")
    r.tags += mm_tag(revcomp(r.seq))
    recs.append(r)

    # Invalid combinations: MM without ML, ML too short, MM skipping past
    # the available Cs
    r = aligned("bm04", 0, "chr1", 450, "30M")
    r.tags.append(mm_tag(r.seq)[0])
    recs.append(r)
    r = aligned("bm05", 0, "chr1", 500, "30M")
    t = mm_tag(r.seq)
    r.tags += [t[0], "ML:B:C,7"]
    recs.append(r)
    r = aligned("bm06", 0, "chr1", 550, "30M")
    r.tags += ["MM:Z:C+m,0,99;", "ML:B:C,10,20"]
    recs.append(r)

    # Unmapped read with mods
    s = rand_seq(30)
    recs.append(Rec("bm07", 4, "*", 0, 0, "*", s, rand_qual(30), mm_tag(s)))
    return recs


# ---- Duplicates, for markdup ---------------------------------------------

def illumina(tile, x, y, umi=None):
    """machine:run:flowcell:lane:tile:x:y[:UMI]"""
    return f"A00123:8:HFLTDDSXX:1:{tile}:{x}:{y}" + (f":{umi}" if umi else "")


def add_ms(recs):
    """ms (mate score) tags as fixmate -m would add: the sum of the mate's
    base qualities >= 15."""
    groups = {}
    for r in recs:
        groups.setdefault(r.qname, []).append(r)
    for rs in groups.values():
        prim = {r.flag & 0xc0: r for r in rs if not r.flag & 0x900}
        for r in rs:
            mate = prim.get((r.flag & 0xc0) ^ 0xc0)
            if mate is None or mate.qual == "*":
                continue
            ms = sum(ord(c) - 33 for c in mate.qual if ord(c) - 33 >= 15)
            r.tags.append(f"ms:i:{ms}")


def dup_reads():
    t1 = (("chr1", 100, "30M"), ("chr1", 250, "30M"))
    umi = "ACGTAC+TTGACA"
    recs = []
    # Duplicate pairs: two optical (same tile, close x/y), one optical
    # without a UMI, one on another tile, one with a soft clip but the same
    # unclipped start, and one with a different UMI
    recs += make_pair(illumina(1101, 1000, 2000, umi), *t1)
    recs += make_pair(illumina(1101, 1010, 2030, umi), *t1)
    recs += make_pair(illumina(1101, 1040, 1990), *t1)
    recs += make_pair(illumina(2204, 5000, 9000, umi), *t1)
    recs += make_pair(illumina(1101, 1500, 2500, umi),
                      ("chr1", 103, "3S27M"), ("chr1", 250, "30M"))
    recs += make_pair(illumina(1101, 1020, 2010, "GGGGGG+TTGACA"), *t1)

    # Same positions, other read group (only a dup without --use-read-groups)
    for r in make_pair(illumina(1102, 300, 400, umi), *t1):
        r.tags = ["RG:Z:rg2" if t == "RG:Z:rg1" else t for t in r.tags]
        recs.append(r)

    # Same positions, BX barcodes (one matching, one not)
    for i, bx in enumerate(("AAAA-1", "AAAA-1", "CCCC-1")):
        for r in make_pair(illumina(1105, 100 + i, 100), ("chr1", 330, "30M"),
                           ("chr1", 450, "30M")):
            r.tags.append(f"BX:Z:{bx}")
            recs.append(r)

    # A QC-failed duplicate, and one already marked with old dt/do tags
    recs += make_pair(illumina(1106, 10, 10), *t1, flags=0x200)
    for r in make_pair(illumina(1106, 20, 20), *t1, flags=0x400):
        r.tags += ["dt:Z:SQ", "do:Z:" + illumina(1101, 1000, 2000, umi)]
        recs.append(r)

    # Supplementary and secondary of a duplicate
    name = illumina(1101, 1010, 2030, umi)
    recs.append(aligned(name, 0x1 | 0x40 | 0x800, "chr2", 100, "10H20M",
                        tags=["SA:Z:chr1,100,+,20S10M,60,0;"]))
    recs.append(aligned(name, 0x1 | 0x80 | 0x100, "chr2", 150, "30M", mapq=0))

    # Single-ended reverse-strand duplicates, differing only in end clipping
    for i, cig in enumerate(("30M", "27M3S", "30M")):
        recs.append(aligned(illumina(1103, 100 + i * 5, 100), 16, "chr1",
                            500, cig, rg="rg1"))

    # Pairs with an unmapped mate, and pairs across chromosomes
    for i in range(2):
        recs += make_pair(illumina(1104, 50 + i, 60), ("chr1", 600, "30M"),
                          None, mate_unmapped=True)
        recs += make_pair(illumina(1107, 70 + i, 80), ("chr1", 700, "30M"),
                          ("chr2", 300, "30M"))

    # Names markdup can't parse for optical duplicates
    recs += make_pair("odd_name", ("chr2", 400, "30M"), ("chr2", 500, "30M"))
    recs += make_pair("x:y:z", ("chr2", 400, "30M"), ("chr2", 500, "30M"))

    add_ms(recs)
    return by_pos(recs)


# ---- Barcode / UMI tags, for fastq -i and the stats barcode sections ----

def barcode_reads():
    recs = []
    for i in range(4):
        pair = make_pair(f"bc{i}", ("chr1", 100 + 40 * i, "30M"),
                         ("chr1", 250 + 40 * i, "30M"))
        dual = i % 2 == 1
        bc = rand_seq(8) + ("-" + rand_seq(8) if dual else "")
        qt = rand_qual(8) + (" " + rand_qual(8) if dual else "")
        if i == 3:
            bc, qt = "ACNNTG-GG", rand_qual(5)   # lengths don't match
        tags = [f"BC:Z:{bc}", f"QT:Z:{qt}",
                f"RX:Z:{rand_seq(6)}-{rand_seq(6)}", f"QX:Z:{rand_qual(6)} {rand_qual(6)}",
                f"OX:Z:{rand_seq(8)}", f"BZ:Z:{rand_qual(8)}",
                f"CR:Z:{rand_seq(16)}", f"CY:Z:{rand_qual(16)}"]
        for r in pair:
            r.tags = [t for t in r.tags if not t.startswith("MC:")] + tags
        recs += pair
    return recs


# ---- Amplicons, for ampliconclip / ampliconstats -------------------------

# name: (ref, left primer start, end, right primer start, end), 0-based
AMPLICONS = {
    "amp1": ("chr1", 90, 110, 230, 250),
    "amp2": ("chr1", 210, 230, 380, 400),
    "amp3": ("chr2", 100, 120, 260, 280),
}


def primers_bed():
    lines = []
    for name, (ref, ls, le, rs, re_) in AMPLICONS.items():
        lines.append(f"{ref}\t{ls}\t{le}\t{name}_LEFT\t60\t+")
        lines.append(f"{ref}\t{rs}\t{re_}\t{name}_RIGHT\t60\t-")
    lines.append("chr1\t205\t225\tamp2_LEFT_alt\t60\t+")  # alternate primer
    return "\n".join(lines) + "\n"


def amplicon_reads():
    recs = []
    for name, (ref, ls, le, rs, re_) in AMPLICONS.items():
        for i in range(3):
            left = (ref, ls + 1 + i, "40M" if i != 2 else "3S37M")
            right = (ref, re_ - 39 - i, "40M")
            recs += make_pair(f"{name}_t{i}", left, right)
    # A template spanning two amplicons (wrong primer pairing)
    recs += make_pair("chimera", ("chr1", 91, "40M"), ("chr1", 361, "40M"))
    return by_pos(recs)


def conflict_sam():
    """A merge partner whose header disagrees with the main datasets:
    a different chr1 length, a clashing @RG ID, a colliding @PG ID."""
    h = ["@HD\tVN:1.6\tSO:coordinate",
         "@SQ\tSN:chr1\tLN:950", "@SQ\tSN:chrX\tLN:300",
         "@RG\tID:rg1\tSM:other\tLB:L9",
         "@PG\tID:bwa\tPN:bwa\tVN:0.7.18",
         "@CO\tconflicting header"]
    recs = [aligned("c01", 0, "chr1", 120, "30M", rg="rg1"),
            aligned("c02", 16, "chr1", 300, "30M", rg="rg1")]
    return "\n".join(h) + "\n" + "".join(r.line() + "\n" for r in recs)


def reheader_variants():
    base = header()
    lines = base.rstrip("\n").split("\n")
    return {
        "addco":   base + "@CO\tnew comment\n",
        "rename":  base.replace("SN:chr2", "SN:chrTwo"),
        "fewer":   "\n".join(l for l in lines if "SN:chrM" not in l) + "\n",
        "badline": base + "@XX\tnot a valid header line\n@SQ\tLN:10\n",
    }



# ---- Group 4 inputs: FASTQ for import/checksum, alt file, checksum report

def fastq_sets():
    """Named (READ1, READ2) FASTQ text pairs in assorted header styles."""
    def rec(name, seq=None, n=30):
        seq = seq or rand_seq(n)
        return f"@{name}\n{seq}\n+\n{rand_qual(len(seq))}\n"
    sets = {}
    sets["casava"] = (
        "".join(rec(f"M1:8:FC:1:1101:{1000+i}:{2000+i} 1:{'NY'[i == 2]}:0:ACGTAC+TTGACA")
                for i in range(4)),
        "".join(rec(f"M1:8:FC:1:1101:{1000+i}:{2000+i} 2:{'NY'[i == 2]}:0:ACGTAC+TTGACA")
                for i in range(4)))
    sets["slash"] = ("".join(rec(f"q{i}/1") for i in range(4)),
                     "".join(rec(f"q{i}/2") for i in range(4)))
    sets["sra"] = ("".join(rec(f"SRR001.{i} orig{i}/1") for i in range(3)),
                   "".join(rec(f"SRR001.{i} orig{i}/2") for i in range(3)))
    sets["tags"] = (
        "".join(rec(f"t{i}\tBC:Z:ACGT\tRG:Z:x\tXI:i:{i}") for i in range(3)),
        "".join(rec(f"t{i}\tBC:Z:ACGT\tBZ:Z:IIII") for i in range(3)))
    sets["umi"] = (
        "".join(rec(f"M1:8:FC:1:1101:{i}:{i}:ACGTAC") for i in range(3)),
        "".join(rec(f"M1:8:FC:1:1101:{i}:{i}:ACGTAC") for i in range(3)))
    # Index reads, as for --i1/--i2
    sets["index"] = ("".join(rec(f"M1:8:FC:1:1101:{1000+i}:{2000+i}", n=6) for i in range(4)),
                     "".join(rec(f"M1:8:FC:1:1101:{1000+i}:{2000+i}", n=8) for i in range(4)))
    # Mismatched pairs: different counts, and a truncated record
    sets["unpaired"] = ("".join(rec(f"u{i}") for i in range(4)),
                        "".join(rec(f"u{i}") for i in range(2)) + "@u9\nACGT\n+\n")
    return {k: (a.encode(), b.encode()) for k, (a, b) in sets.items()}


def alt_file():
    """A bwa-style .alt file (SAM-like lines) for dict --alt."""
    return ("@SQ\tSN:chr1\tLN:900\n"
            "seq1\t0\tchr1\t100\t60\t50M\t*\t0\t0\t*\t*\n"
            "seq2\t16\tchr1\t300\t60\t10S40M\t*\t0\t0\t*\t*\n"
            "absent\t0\tchrQ\t1\t60\t5M\t*\t0\t0\t*\t*\n").encode()


def checksum_report():
    """A checksum report in the format "samtools checksum" writes, for -m."""
    return ("# Checksum for file: x.bam\n"
            "# Aux tags:          BC,FI,QT,RT,TC\n"
            "# BAM flags:         PAIRED,READ1,READ2\n"
            "\n"
            "# Group\tQC\tcount\tflag+seq\t+name\t+qual\t+aux\tcombined\n"
            "all\tall\t10\t71169bbb\t633fd9f7\t2a2e693f\t71169bbb\t09d03ed4\n"
            "all\tpass\t9\t1a2b3c4d\t5e6f7a8b\t9c0d1e2f\t3a4b5c6d\t7e8f9a0b\n"
            "rg1\tall\tnotanumber\tzz\n").encode()

# ---- Deep coverage with variants, for consensus --------------------------

def deep_reads():
    ref = REF["chr1"]
    alt = lambda p: rng.choice([b for b in BASES if b != ref[p].upper()])
    hom = {170: ("snp", alt(170)), 205: ("ins", "AT")}
    het = dict(hom)
    het.update({150: ("snp", alt(150)), 190: ("del", 2)})

    recs = []
    for i in range(32):
        edits = het if i % 2 else hom
        p0 = rng.randint(99, 194)
        while p0 in edits and edits[p0][0] == "del":
            p0 += 1
        seq, ops, p = [], [], p0

        def op(c, n=1):
            if ops and ops[-1][1] == c:
                ops[-1][0] += n
            else:
                ops.append([n, c])

        while len(seq) < 25:
            e = edits.get(p)
            if e and e[0] == "del":
                op("D", e[1])
                p += e[1]
                continue
            b = e[1] if e and e[0] == "snp" else ref[p].upper()
            if rng.random() < 0.01:
                b = rng.choice(BASES)
            seq.append(b)
            op("M")
            if e and e[0] == "ins":
                for c in e[1][:25 - len(seq)]:
                    seq.append(c)
                    op("I")
            p += 1
        cigar = "".join(f"{n}{c}" for n, c in ops)
        recs.append(aligned(f"d{i:02d}", 16 if i % 3 == 0 else 0, "chr1",
                            p0 + 1, cigar, mapq=rng.choice([0, 20, 40, 60]),
                            seq="".join(seq), rg="rg1"))
    return by_pos(recs)


# ---------------------------------------------------------------------------
# Datasets and seeds

def build_sam():
    mixed = mixed_reads()
    unsorted = [mixed[i] for i in (7, 2, 20, 0, 15, 4, 11, 22)]
    d = {
        "mixed":       sam(header(), mixed),
        "pairs":       sam(header(), by_pos(pair_reads())),
        "deep":        sam(header(), deep_reads()),
        "collated":    sam(header("unsorted", "query"), collated_reads()),
        "basemods":    sam(header("unsorted", "query"), basemod_reads()),
        "unsorted":    sam(header("unsorted"), unsorted),
        "dups":        sam(header(), dup_reads()),
        "basemods_sorted": sam(header(), by_pos(basemod_reads())),
        "barcodes":    sam(header("unsorted", "query"), barcode_reads()),
        "amplicons":   sam(header(), amplicon_reads()),
        "minimal":     sam(header(), mixed[:1]),
        "header_only": header(),
    }
    return {k: v.encode() for k, v in d.items()}


def binary_versions(name, text):
    """BAM and CRAM copies via samtools, if available."""
    out = {"sam": text}
    st = shutil.which("samtools")
    if not st:
        return out
    path = f"{OUT}/sam/{name}.sam"
    for fmt, args in (("bam", ["-b"]), ("cram", ["-O", "cram,no_ref"])):
        r = subprocess.run([st, "view", "--no-PG", *args, path],
                           capture_output=True)
        if r.returncode == 0:
            out[fmt] = r.stdout
        else:
            print(f"warning: samtools {fmt} conversion of {name} failed:",
                  r.stderr.decode().strip(), file=sys.stderr)
    return out


def seed(fuzzer, name, prefix, payload):
    write(f"{OUT}/seeds/{fuzzer}/{name}", bytes(prefix) + payload)


# Control-byte prefixes: all-zero (defaults) plus a few hand-picked sets.
# None of the view sets enable its file-writing options (-o, -U,
# --save-counts).
PREFIX = {
    "view":      {"0": [0] * 5, "a": [0x00, 0x14, 0x01, 0x41, 0x00]},
    # a: -f fastq -a --show-del yes;  b: -m simple -A -q
    "consensus": {"0": [0] * 6, "a": [0x12, 0x01, 0, 0, 0, 0],
                  "b": [0x44, 0, 0x10, 0, 0, 0]},
    # a: -e -r;  b: -A -r with ref.fa.gz;  c: ref_short.fa
    "calmd":     {"0": [0, 0], "a": [0x42, 0], "b": [0x41, 0x01], "c": [0, 0x02]},
    # a: -c -m;  b: -M;  c: -r -z off
    "fixmate":   {"0": [0, 0], "a": [0x0c, 0], "b": [0x10, 0], "c": [0x01, 0x01]},
    # a: -t --duplicate-count -d 100 -s
    # b: -S --use-read-groups --barcode-name --move-umi-to-tag
    # c: -r -d 2500 --json -f /dev/null, Illumina --read-coords, txy order
    "markdup":   {"0": [0] * 4, "a": [0x18, 0x12, 0, 0],
                  "b": [0x02, 0, 0x20, 0x06], "c": [0x01, 0x63, 0x05, 0]},
    # a: -n -O sam;  b: --template-coordinate;  c: -I ref.fa -M -w 1
    "sort":      {"0": [0] * 3, "a": [0x01, 0x04, 0], "b": [0x06, 0, 0],
                  "c": [0x87, 0x01, 0]},
    # default, json, tsv
    "flagstat":  {"0": [0], "a": [2], "b": [3]},
    # a: -r ref.fa -t targets.txt -S RG;  b: -d -x -p -c 1,10,1 --ref-stats
    "stats":     {"0": [0] * 4, "a": [0x80, 0, 0x40, 0x01],
                  "b": [0x61, 0x01, 0, 0x04]},
    # a: -aa -b regions.bed, two inputs;  b: -H -J -s -q 0 -Q 0
    "depth":     {"0": [0] * 3, "a": [0x07, 0, 0x20], "b": [0xc8, 0x05, 0]},
    # a: -f ref.fa -l regions.bed --output-QNAME, two inputs, and
    #    --output-extra QNAME,FLAG,RNAME + RG,NM,MD,XB,XH
    # b: -A -x -a -M
    # c: --output-extra with every field and tag, plus invalid names
    "mpileup":   {"0": [0] * 7, "a": [0, 0x05, 0, 0x0c, 0, 0x0f, 0x37],
                  "b": [0x02, 0, 0x30, 0, 0x40, 0, 0],
                  "c": [0, 0, 0, 0x18, 0, 0xfe, 0xff]},
    # a: -t, -1/-2;  b: fasta;  c: -i --index-format i* --i1/--i2
    "fastq":     {"0": [0] * 4, "a": [0x10, 0, 0x01, 0], "b": [0x01, 0, 0, 0],
                  "c": [0x80, 0, 0x40, 0x01]},
    # a: -f -r 1;  b: -n 64 --no-PG
    "collate":   {"0": [0] * 2, "a": [0x80, 0x01], "b": [0x60, 0x04]},
    # one per read-group spec, plus -m overwrite_all -w
    "addreplacerg": dict([(f"r{i}", [i, 0]) for i in range(8)] +
                         [("m", [0x2b, 0x01])]),
    # a: -x ^ --reject-PG bwa --no-RG;  b: --keep-tag RG,BC --dupflag -O sam
    "reset":     {"0": [0] * 2, "a": [0x52, 0], "b": [0x04, 0x03]},
    # a: -s -1;  b: -T ref.fa -u
    "depad":     {"0": [0], "a": [0x11], "b": [0x28]},
    # a: -A -F -k 13;  b: -b prefix -q 0 -Q 0
    "phase":     {"0": [0] * 2, "a": [0x0b, 0], "b": [0x50, 0x01]},
    # a: -Q 13 -i 5;  b: -f ref.fa, negative emissions
    "targetcut": {"0": [0] * 2, "a": [0x0a, 0], "b": [0xa0, 0x06]},
    # header lines / record counts
    "head":      {"0": [0], "a": [0x05], "b": [0x0f]},
    # -s;  -S
    "rmdup":     {"0": [0], "a": [0x01], "b": [0x02]},
}

USES = {
    "view":      ["mixed", "pairs", "deep", "collated", "basemods",
                  "unsorted", "dups", "barcodes", "minimal", "header_only"],
    # consensus rejects unsorted input up front
    "consensus": ["mixed", "pairs", "deep", "minimal", "header_only"],
    "calmd":     ["mixed", "pairs", "deep", "unsorted", "minimal"],
    "fixmate":   ["collated", "basemods", "unsorted", "minimal", "header_only"],
    # pairs has MC but no ms tags; markdup rejects unsorted input up front
    "markdup":   ["dups", "pairs", "mixed", "minimal", "header_only"],
    "sort":      ["mixed", "pairs", "deep", "collated", "basemods", "unsorted",
                  "dups", "barcodes", "minimal", "header_only"],
    "flagstat":  ["mixed", "pairs", "collated", "dups", "minimal", "header_only"],
    "stats":     ["mixed", "pairs", "deep", "dups", "unsorted", "barcodes",
                  "minimal", "header_only"],
    # depth and mpileup reject unsorted input up front
    "depth":     ["mixed", "pairs", "deep", "dups", "minimal", "header_only"],
    "mpileup":   ["mixed", "pairs", "deep", "dups", "basemods_sorted",
                  "minimal", "header_only"],
    "fastq":     ["collated", "basemods", "barcodes", "pairs", "mixed",
                  "minimal", "header_only"],
    "collate":   ["mixed", "pairs", "dups", "basemods", "unsorted", "minimal",
                  "header_only"],
    "addreplacerg": ["mixed", "pairs", "barcodes", "minimal", "header_only"],
    "reset":     ["mixed", "pairs", "dups", "basemods", "unsorted", "minimal",
                  "header_only"],
    "depad":     ["mixed", "pairs", "minimal", "header_only"],
    "phase":     ["deep", "mixed", "pairs"],
    "targetcut": ["deep", "mixed", "pairs"],
    "head":      ["mixed", "pairs", "deep", "minimal", "header_only"],
    "rmdup":     ["pairs", "dups", "mixed"],
}


# File-backed fuzzers (fuzz_one_file_input in fuzz_common.h).  Inputs are
# [file byte][tool bytes][16 region bytes][data].
#
# File byte variants: BAM + default index + 1 region; CRAM + 4 regions;
# BGZF SAM + alternate index name (via ##idx## / -X) + 3 regions.
FILE_BYTES = {"bam": 0x00, "cram": 0x31, "samgz": 0x2a}

# name:beg-end on ref 0; {name}:beg-end on ref 1; "*"; name:beg- on ref 2
REGION_BYTES = bytes([0x02 | 0 << 3, 0x01, 0xC0 | 40, 0x00,
                      0x04 | 1 << 3, 0xC0 | 10, 0x40, 0x01,
                      0x07,          0x00, 0x00, 0x00,
                      0x03 | 2 << 3, 0xC0 | 60, 0x00, 0x00])

FILE_TOOL = {
    # a: -c;  b: -c -m 14 -o alt;  c: -M
    "index":    {"0": [0], "a": [0x02], "b": [0x13], "c": [0x20]},
    # a: -m -A -w 1 -r;  b: -H -b list
    "coverage": {"0": [0] * 3, "a": [0x40, 0x09, 0x10], "b": [0, 0x02, 0x40]},
    # a: -j -c -d 0;  b: -X with two inputs
    "bedcov":   {"0": [0] * 2, "a": [0xc0, 0x01], "b": [0, 0x30]},
    # One per command (view, view -M, view -L, view -P, depth, mpileup,
    # consensus, stats), plus view -M -X and depth -b with -X
    "regions":  dict([(f"t{i}", [i]) for i in range(8)] +
                     [("x", [0x09]), ("y", [0x2c])]),
    # a: -q -r;  b: -e (embedded reference)
    "reference": {"0": [0], "a": [0x06], "b": [0x01]},
    # a: -v;  b: -v -e
    "cram_size": {"0": [0], "a": [0x01], "b": [0x03]},
    # a: -h -i -T SM;  b: -X with two inputs, -f ref.fa
    "samples":   {"0": [0], "a": [0x07], "b": [0x70]},
    # cut 0, 1, 27, 28 (the EOF block) bytes; -v -v
    "quickcheck": {"0": [0], "a": [0x20], "b": [0x60], "c": [0x83]},
    # text at a region;  HTML, width 40, with ref.fa
    "tview":     {"0": [0], "a": [0x02], "b": [0x0f]},
}

FILE_USES = {
    "index":    ["mixed", "pairs", "deep", "dups", "minimal", "header_only"],
    "coverage": ["mixed", "pairs", "deep", "dups", "minimal"],
    "bedcov":   ["mixed", "pairs", "deep", "dups", "minimal"],
    "regions":  ["mixed", "pairs", "deep", "dups", "basemods_sorted"],
    "reference": ["mixed", "pairs", "deep"],
    "cram_size": ["mixed", "pairs", "deep", "header_only"],
    "samples":   ["mixed", "pairs", "deep", "dups"],
    "quickcheck": ["mixed", "minimal", "header_only", "deep"],
    "tview":     ["mixed", "deep", "pairs"],
}


# File-backed fuzzers that take an extra payload (FUZZ_EXTRA):
# [file byte][tool bytes][16 region bytes][2-byte length][extra][data]
EXTRA_TOOL = {
    # a: -r -c;  b: -p, extra as data: URL, -R, -O sam;  c: three inputs, -h
    "merge":         {"0": [0] * 3, "a": [0x18, 0, 0], "b": [0x20, 0x09, 0x01],
                      "c": [0, 0x07, 0]},
    # a: main twice;  b: main + extra + main, -h
    "cat":           {"0": [0], "a": [0x01], "b": [0x06]},
    # a: -d RG -u;  b: -d BX -u -h --output-fmt sam
    "split":         {"0": [0] * 2, "a": [0x05, 0], "b": [0x0e, 0x01]},
    # a: -P;  b: -i (in place; CRAM only)
    "reheader":      {"0": [0], "a": [0x01], "b": [0x02]},
    # a: hard clip, both ends, strand, lengths/tolerance, side outputs
    # b: --clipped --fail with --rejects-file
    "ampliconclip":  {"0": [0] * 3, "a": [0x07, 0x55, 0x0e], "b": [0x18, 0, 0x02]},
    # a: -d 1 -s;  b: two inputs, -S
    "ampliconstats": {"0": [0] * 3, "a": [0, 0x40, 0x04], "b": [0, 0, 0x88]},
    # a: -a GRCh38 -A -s "Homo sapiens" --alt;  b: -H, odd -u
    "dict":          {"0": [0] * 2, "a": [0x15, 0x01], "b": [0xc8, 0]},
    # a: -1/-2 -i;  b: --i1/--i2 -i -U, -r with tabs, --order rn:12
    # c: -T * -N, -R "";  d: positional pair -O sam
    "import":        {"0": [0] * 3, "a": [0x0c, 0, 0], "b": [0x1d, 0x42, 0],
                      "c": [0x23, 0x0f, 0x05], "d": [0x01, 0x01, 0x01]},
    # a: -a -O;  b: -m on the report;  c: second input, -T -v
    "checksum":      {"0": [0] * 3, "a": [0x11, 0, 0], "b": [0, 0, 0x10],
                      "c": [0x40, 0x5b, 0x28]},
}
EXTRA_FILE_BYTES = {"bam": 0x00, "cram": 0x31}


def faidx_sources():
    seqs = [("seq1", rand_seq(150)), ("seq2", rand_seq(61)),
            ("chr3:1-10", rand_seq(95)), ("seq4", rand_seq(7))]
    desc = {"seq1": "first sequence", "chr3:1-10": "colon in name"}
    irregular = (">a\nACGTACGT\nACG\nACGTACGTAC\n>empty\n>a dup\nACGT\n\n"
                 ">b\nAC GT\r\nACGT")  # no final newline
    fq = [(f"q{i} comment", rand_seq(n), rand_qual(n))
          for i, n in enumerate((20, 35, 1, 40, 28, 33), 1)]
    fqw = [(f"w{i}", rand_seq(n), rand_qual(n)) for i, n in enumerate((31, 45, 15), 1)]
    return {
        "multi.fa":     fasta(seqs, 60, desc=desc).encode(),
        "wrap10.fa":    fasta(seqs, 10).encode(),
        "oneline.fa":   fasta(seqs[:2], 0).encode(),
        "crlf.fa":      fasta(seqs, 20, eol="\r\n").encode(),
        "reads.fq":     fastq(fq).encode(),
        "wrapped.fq":   fastq(fqw, 15).encode(),
        "irregular.fa": irregular.encode(),   # build-only; not indexable
    }


def faidx_prefixes(is_fq):
    fq = 0x08 if is_fq else 0
    # 4 regions: name:beg-end, {name}:beg-end, name:beg, name:beg-
    regions = bytes([0x02 | 0 << 3, 0x01, 0xC0 | 40, 0,
                     0x04 | 1 << 3, 0xC0 | 10, 0x40, 1,
                     0x01 | 2 << 3, 0x80 | 2, 0, 0,
                     0x03 | 3 << 3, 0xC0 | 60, 0, 0])
    return {
        "0": bytes([fq, 0, 0]) + bytes(16),
        "a": bytes([fq | 0x03, 0, 0]) + regions,
        # regions via -r file, -c, -i, --mark-strand sign, -n 7
        "b": bytes([fq | 0x87, 0x01 | 0x04 | 0x18, 0]) + regions,
    }


def main():
    # Runtime files
    ref_fa = fasta(list(REF.items()), 60).encode()
    short = [(n, s[:len(s) // 2]) for n, s in REF.items() if n != "chrM"]
    short_fa = fasta(short, 60).encode()
    ref_gz, ref_gzi = bgzf(ref_fa, 1000)
    rt = f"{OUT}/runtime"
    write(f"{rt}/ref.fa", ref_fa)
    write(f"{rt}/ref.fa.fai", make_fai(ref_fa))
    write(f"{rt}/ref.fa.gz", ref_gz)
    write(f"{rt}/ref.fa.gz.fai", make_fai(ref_fa))  # uncompressed offsets
    write(f"{rt}/ref.fa.gz.gzi", ref_gzi)
    write(f"{rt}/ref_short.fa", short_fa)
    write(f"{rt}/ref_short.fa.fai", make_fai(short_fa))
    write(f"{rt}/qcal.txt", "# QUAL value substitution undercall overcall\n" +
          "".join(f"QUAL {q} {q} {max(q - 5, 0)} {max(q - 8, 0)}\n"
                  for q in range(0, 45, 5)))
    write(f"{rt}/qname.txt", "r01\nr05\np01\nd03\nnot_present\n")
    write(f"{rt}/RG.txt", "rg1\n")
    write(f"{rt}/TAG.txt", "0\n1\n2\n")
    write(f"{rt}/skip.bed", "chr1\t100\t300\nchr2\t0\t50\nchrM\t0\t120\n")
    # depth -b, mpileup -l (BED: 0-based, half open)
    write(f"{rt}/regions.bed",
          "chr1\t90\t200\nchr2\t0\t100\nchr3:1-10\t10\t20\nchrX\t0\t10\n")
    # stats -t (chr, from, to; 1-based inclusive)
    write(f"{rt}/targets.txt", "chr1\t100\t300\nchr2\t1\t50\nchrM\t1\t120\n")
    # mpileup -G (read groups to exclude)
    write(f"{rt}/rg_exclude.txt", "rg2\n")
    # (ampliconclip / ampliconstats take their primer BED from the fuzz
    # input, so there is no runtime primers file)

    # SAM-based seeds
    data = build_sam()
    for name, text in data.items():
        write(f"{OUT}/sam/{name}.sam", text)
    versions = {name: binary_versions(name, text) for name, text in data.items()}
    for fuzzer, names in USES.items():
        for name in names:
            for fmt, payload in versions[name].items():
                for v, prefix in PREFIX[fuzzer].items():
                    seed(fuzzer, f"{name}.{fmt}-{v}", prefix, payload)

    # File-backed seeds.  Only SAM payloads: the driver re-encodes them as
    # BAM/CRAM/BGZF SAM itself (the "raw" form needs BAM/CRAM payloads).
    for fuzzer, names in FILE_USES.items():
        for name in names:
            for fv, fb in FILE_BYTES.items():
                for tv, tool in FILE_TOOL[fuzzer].items():
                    seed(fuzzer, f"{name}-{fv}-{tv}",
                         [fb] + tool + list(REGION_BYTES), data[name])
            for fmt in ("bam", "cram"):
                if fmt in versions[name]:   # raw form, file byte 0x03
                    seed(fuzzer, f"{name}.{fmt}-raw", [0x03] +
                         list(FILE_TOOL[fuzzer].values())[0] +
                         list(REGION_BYTES), versions[name][fmt])

    # Extra-payload seeds: (fuzzer, main dataset, extra label, extra bytes)
    bed = primers_bed().encode()
    bad_bed = b"chr1\t250\t90\tbackwards\t60\t+\nchrQ\t1\t2\nnot a bed line\n"
    hdr_only = (header() + "@RG\tID:rg3\tSM:s3\n").encode()
    extras = [("merge", "mixed", "pairs", data["pairs"]),
              ("merge", "pairs", "conflict", conflict_sam().encode()),
              ("merge", "dups", "mixed", data["mixed"]),
              ("cat", "mixed", "pairs", data["pairs"]),
              ("cat", "dups", "dups", data["dups"])]
    extras += [("split", m, "hdr", hdr_only)
               for m in ("mixed", "pairs", "barcodes", "dups")]
    extras += [("reheader", m, lbl, h.encode())
               for m in ("mixed", "minimal")
               for lbl, h in reheader_variants().items()]
    extras += [(f, m, lbl, x)
               for f in ("ampliconclip", "ampliconstats")
               for m in ("amplicons", "pairs")
               for lbl, x in (("primers", bed), ("badbed", bad_bed))]
    # Group 4: (fuzzer, label, main payload, extra payload)
    fq = fastq_sets()
    fa = faidx_sources()
    report = checksum_report()
    g4x = [("dict", f"{n}+alt", fa[n], alt_file())
           for n in ("multi.fa", "wrap10.fa", "crlf.fa", "irregular.fa")]
    g4x += [("dict", "multi.fa.gz+alt", bgzf(fa["multi.fa"], 100)[0], alt_file())]
    g4x += [("import", f"{k}", a, b) for k, (a, b) in fq.items()]
    g4x += [("checksum", f"{m}+report", data[m], report)
            for m in ("mixed", "pairs", "barcodes", "collated")]
    g4x += [("checksum", "pairs+mixed", data["pairs"], data["mixed"])]
    g4x += [("checksum", f"fq_{k}+report", a, report) for k, (a, b) in fq.items()
            if k in ("casava", "tags")]

    for fuzzer, main, lbl, x in extras:
        g4x.append((fuzzer, f"{main}+{lbl}", data[main], x))
    for fuzzer, label, main_bytes, x in g4x:
        for fv, fb in EXTRA_FILE_BYTES.items():
            for tv, tool in EXTRA_TOOL[fuzzer].items():
                seed(fuzzer, f"{label}-{fv}-{tv}",
                     [fb] + tool + list(REGION_BYTES) +
                     list(struct.pack("<H", len(x))) + list(x), main_bytes)

    # quickcheck on raw text and on damaged non-alignment data
    for tv, tool in FILE_TOOL["quickcheck"].items():
        for name in ("mixed", "header_only"):
            seed("quickcheck", f"{name}-raw-{tv}", [0x03] + tool +
                 list(REGION_BYTES), data[name])
        seed("quickcheck", f"notsam-raw-{tv}", [0x03] + tool +
             list(REGION_BYTES), b"this is not alignment data\n")

    # FASTA/FASTQ seeds
    src = faidx_sources()
    for name in ("multi.fa", "reads.fq"):
        src[name + ".gz"] = bgzf(src[name], 100)[0]
    for name, text in src.items():
        write(f"{OUT}/fasta/{name}", text)

    for name, text in src.items():
        for v, prefix in faidx_prefixes(".fq" in name).items():
            seed("faidx", f"{name}-{v}", prefix, text)

    def idx_seed(label, prefix, fai, gzi, payload):
        seed("faidx_idx", label, prefix,
             struct.pack("<HH", len(fai), len(gzi)) + fai + gzi + payload)

    for name, text in src.items():
        if name == "irregular.fa" or name.endswith(".gz"):
            continue
        plain_fai = make_fai(text)
        for v, prefix in faidx_prefixes(".fq" in name).items():
            idx_seed(f"{name}-{v}", prefix, plain_fai, b"", text)
            if name in ("multi.fa", "reads.fq"):
                gz, gzi = bgzf(text, 100)
                idx_seed(f"{name}.gz-{v}", prefix, plain_fai, gzi, gz)

    # Deliberately mismatched data and indices
    p = faidx_prefixes(False)["a"]
    idx_seed("mismatch-fai", p, make_fai(src["wrap10.fa"]), b"", src["multi.fa"])
    idx_seed("mismatch-fmt", faidx_prefixes(True)["a"],
             make_fai(src["multi.fa"]), b"", src["reads.fq"])
    gz, _ = bgzf(src["multi.fa"], 100)
    _, other_gzi = bgzf(src["multi.fa"], 37)
    idx_seed("mismatch-gzi", p, make_fai(src["multi.fa"]), other_gzi, gz)

    write(f"{OUT}/README.txt", __doc__.strip() + "\n")
    for d in sorted(os.listdir(f"{OUT}/seeds")):
        print(f"seeds/{d}: {len(os.listdir(f'{OUT}/seeds/{d}'))} files")
    if not shutil.which("samtools"):
        print("samtools not found: SAM seeds only (rerun with it on PATH "
              "for BAM/CRAM)")


if __name__ == "__main__":
    main()
