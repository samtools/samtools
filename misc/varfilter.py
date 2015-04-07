#!/usr/bin/env python
#
#    Copyright (C) 2009, 2010 Genome Research Ltd.
#
#    Author: Aylwyn Scally <as6@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

# Author: lh3, converted to python and modified to add -C option by Aylwyn Scally
#
# About:
#   varfilter.py is a port of Heng's samtools.pl varFilter script into
#   python, with an additional -C INT option. This option sets a minimum
#   consensus score, above which the script will output a pileup line
#   wherever it _could have_ called a variant, even if none is actually
#   called (i.e. hom-ref positions). This is important if you want to
#   subsequently merge the calls with those for another individual to get a
#   synoptic view of calls at each site. Without this option, and in all
#   other respects, it behaves like samtools.pl varFilter.
#
#   Aylwyn Scally as6@sanger.ac.uk


# Filtration code:
#
# C low CNS quality (hom-ref only)
# d low depth
# D high depth
# W too many SNPs in a window (SNP only)
# G close to a high-quality indel (SNP only)
# Q low RMS mapping quality (SNP only)
# g close to another indel with higher quality (indel only)
# s low SNP quality (SNP only)
# i low indel quality (indel only)


import sys
import getopt

def usage():
    print '''usage: varfilter.py [options] [cns-pileup]

Options: -Q INT minimum RMS mapping quality for SNPs
         -q INT minimum RMS mapping quality for gaps
         -d INT minimum read depth
         -D INT maximum read depth
         -S INT minimum SNP quality
         -i INT minimum indel quality
         -C INT minimum consensus quality for hom-ref sites

         -G INT min indel score for nearby SNP filtering
         -w INT SNP within INT bp around a gap to be filtered

         -W INT window size for filtering dense SNPs
         -N INT max number of SNPs in a window

         -l INT window size for filtering adjacent gaps

         -p print filtered variants'''

def varFilter_aux(first, is_print):
    try:
        if first[1] == 0:
            sys.stdout.write("\t".join(first[4:]) + "\n")
        elif is_print:
            sys.stderr.write("\t".join(["UQdDWGgsiCX"[first[1]]] + first[4:]) + "\n")
    except IOError:
        sys.exit()

mindepth = 3
maxdepth = 100
gapgapwin = 30
minsnpmapq = 25
mingapmapq = 10
minindelscore = 25
scorefactor = 100
snpgapwin = 10
densesnpwin = 10
densesnps = 2
printfilt = False
minsnpq = 0
minindelq = 0
mincnsq = 0

try:
    options, args = getopt.gnu_getopt(sys.argv[1:], 'pq:d:D:l:Q:w:W:N:G:S:i:C:', [])
except getopt.GetoptError:
    usage()
    sys.exit(2)
for (oflag, oarg) in options:
    if oflag == '-d': mindepth = int(oarg)
    if oflag == '-D': maxdepth = int(oarg)
    if oflag == '-l': gapgapwin = int(oarg)
    if oflag == '-Q': minsnpmapq = int(oarg)
    if oflag == '-q': mingapmapq = int(oarg)
    if oflag == '-G': minindelscore = int(oarg)
    if oflag == '-s': scorefactor = int(oarg)
    if oflag == '-w': snpgapwin = int(oarg)
    if oflag == '-W': densesnpwin = int(oarg)
    if oflag == '-C': mincnsq = int(oarg)
    if oflag == '-N': densesnps = int(oarg)
    if oflag == '-p': printfilt = True
    if oflag == '-S': minsnpq = int(oarg)
    if oflag == '-i': minindelq = int(oarg)

if len(args) < 1:
    inp = sys.stdin
else:
    inp = open(args[0])

# calculate the window size
max_dist = max(gapgapwin, snpgapwin, densesnpwin)

staging = []
for t in (line.strip().split() for line in inp):
    (flt, score) = (0, -1)
    # non-var sites
    if t[3] == '*/*':
        continue
    is_snp = t[2].upper() != t[3].upper()
    if not (is_snp or mincnsq):
        continue
    # clear the out-of-range elements
    while staging:
        # Still on the same chromosome and the first element's window still affects this position?
        if staging[0][4] == t[0] and int(staging[0][5]) + staging[0][2] + max_dist >= int(t[1]):
            break
        varFilter_aux(staging.pop(0), printfilt)

    # first a simple filter
    if int(t[7]) < mindepth:
        flt = 2
    elif int(t[7]) > maxdepth:
        flt = 3
    if t[2] == '*': # an indel
        if minindelq and minindelq > int(t[5]):
            flt = 8
    elif is_snp:
        if minsnpq and minsnpq> int(t[5]):
            flt = 7
    else:
        if mincnsq and mincnsq > int(t[4]):
            flt = 9

    # site dependent filters
    dlen = 0
    if flt == 0:
        if t[2] == '*': # an indel
            # If deletion, remember the length of the deletion
            (a,b) = t[3].split('/')
            alen = len(a) - 1
            blen = len(b) - 1
            if alen>blen:
                if a[0] == '-': dlen=alen
            elif b[0] == '-': dlen=blen

            if int(t[6]) < mingapmapq:
                flt = 1
            # filtering SNPs
            if int(t[5]) >= minindelscore:
                for x in (y for y in staging if y[3]):
                    # Is it a SNP and is it outside the SNP filter window?
                    if x[0] >= 0 or int(x[5]) + x[2] + snpgapwin < int(t[1]):
                        continue
                    if x[1] == 0:
                        x[1] = 5

            # calculate the filtering score (different from indel quality)
            score = int(t[5])
            if t[8] != '*':
                score += scorefactor * int(t[10])
            if t[9] != '*':
                score += scorefactor * int(t[11])
            # check the staging list for indel filtering
            for x in (y for y in staging if y[3]):
              # Is it a SNP and is it outside the gap filter window
                if x[0] < 0 or int(x[5]) + x[2] + gapgapwin < int(t[1]):
                    continue
                if x[0] < score:
                    x[1] = 6
                else:
                    flt = 6
                    break
        else: # a SNP or hom-ref
            if int(t[6]) < minsnpmapq:
                flt = 1
            # check adjacent SNPs
            k = 1
            for x in (y for y in staging if y[3]):
                if x[0] < 0 and int(x[5]) + x[2] + densesnpwin >= int(t[1]) and (x[1] == 0 or x[1] == 4 or x[1] == 5):
                    k += 1

            # filtering is necessary
            if k > densesnps:
                flt = 4
                for x in (y for y in staging if y[3]):
                    if x[0] < 0 and int(x[5]) + x[2] + densesnpwin >= int(t[1]) and x[1] == 0:
                        x[1] = 4
            else: # then check gap filter
                for x in (y for y in staging if y[3]):
                    if x[0] < 0 or int(x[5]) + x[2] + snpgapwin < int(t[1]):
                        continue
                    if x[0] >= minindelscore:
                        flt = 5
                        break

    staging.append([score, flt, dlen, is_snp] + t)

# output the last few elements in the staging list
while staging:
    varFilter_aux(staging.pop(0), printfilt)
