#! /usr/bin/env python3
import sys
import os
import pandas as pd
import numpy as np
import math
import argparse
import pyranges as pr
from tqdm import tqdm
import gzip
import re
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.collections as mc
import random
"""
Find SVs segregating between two genomes, using gene order
"""

"""
Interpret command-line arguments.
"""
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--match_score", "-m", type=int, help='Match score', default=20)
    parser.add_argument("--skip_penalty", "-s", type=int, help="Skip penalty", default=20)
    parser.add_argument('--break_penalty', '-b', type=int, help="Break penalty", default=100)
    parser.add_argument('--translocation_penalty', '-t', type=int, help="Translocation penalty", default=100)
    parser.add_argument("--max_skip", "-S", type=int, help="Maximum number of genes that can be skipped at once", default=100)
    parser.add_argument("--anno1", "-1", help="Input annotation file 1 (gtf/gff3(.gz))", required=True)
    parser.add_argument("--anno2", "-2", help="Input annotation file 2 (gtf/gff3(.gz))", required=True)
    parser.add_argument("--gene_name1", "-n1", help="Gene name field in annotation file 1 (default = gene_name)", 
            default="gene_name", required=False)
    parser.add_argument("--gene_name2", "-n2", help="Gene name field in annotation file 2 (default = gene_name)", 
            default="gene_name", required=False)
    parser.add_argument("--gene_id1", "-i1", help="Gene ID field in annotation file 1 (default = gene_id)",
            default="gene_id", required=False)
    parser.add_argument("--gene_id2", "-i2", help="Gene ID field in annotation file 2 (default = gene_id)",
            default="gene_id", required=False)
    parser.add_argument("--output", "-o", help="Output file name prefix", required=True)
    parser.add_argument("--exclude_chr", "-e", help="Regex pattern (pipe separated) of sequence names to exclude (e.g. \
chrUn|random)", required=False, default="chrUn|random|alt|qp|Scaffold|chrM|MT")
    parser.add_argument("--include_chr", "-i", help="Regex pattern (pipe separated) of sequence names to include (e.g. \
chr)", required=False)
    parser.add_argument("--fai1", "-f1", help="FAI file (from samtools faidx) for ref genome for species 1. \
This is required to make plot-friendly data, but not otherwise.", required=False)
    parser.add_argument("--fai2", "-f2", help="FAI file (from samtools faidx) for ref genome for species 2. \
This is required to make plot-friendly data, but not otherwise.", required=False)
    parser.add_argument("--species1", "-s1", help="Name of species 1 (to be added to plot axis label)",
        required=False, default="Species 1")
    parser.add_argument("--species2", "-s2", help="Name of species 2 (to be added to plot axis label)",
        required=False, default="Species 2")
    parser.add_argument("--mito1", "-m1", help="Name of mitochondrial sequence in species 1 (will allow \
these genes to pass through without filtering)", required=False, default='chrM')
    parser.add_argument("--mito2", "-m2", help="Name of mitochondrial sequence in species 2 (will allow \
these genes to pass through without filtering)", required=False, default='chrM')
    parser.add_argument("--enrichment", "-E", help="Provide a number of random samples for enrichment testing, \
if you wish to test for overrepresentation of genes in inverted regions. e.g. 1000. Defalt = no enrichment testing.",
        required=False, default=None, type=int)
    parser.add_argument("--nofilt", "-nf", help="Make BED files and plot, then quit without filtering \
annotations. Default behavior is to filter both annotations to remove genes not along the expected \
synteny path", action="store_true")
    parser.add_argument("--filt_files", "-F", help="Write out separate files for each set of removed \
genes, to attempt to rescue later", action='store_true')
    return parser.parse_args()

"""
Score emission function for HMM state where both sequences
are moving in the same direction.
"""
def forward_emit(diff, match_score):
    if diff == 0:
        return 0
    elif diff > 0:
        #return match_score - diff
        return match_score/diff
    else:
        return match_score*diff

"""
Score emission function for HMM state where sequences are
moving in opposite directions.
"""
def rev_emit(diff, match_score):
    if diff == 0:
        return 0
    elif diff > 0:
        return -match_score*diff
    else:
        #return match_score + diff
        return -match_score/diff

"""
Runs the HMM to find the likeliest path along a genome, in aligned
sequence blocks.
"""
def find_path(df, match_score, skip_penalty, 
    break_penalty, translocation_penalty, max_skip):
    # Get sizes of chroms
    df = df.merge(-df['chrom1'].value_counts(), 
            left_on='chrom1', right_on='chrom1').rename({'count': 'chrom1size'}, axis=1).merge(
                    -df['chrom2'].value_counts(),
                    left_on='chrom2', right_on='chrom2').rename({'count': 'chrom2size'}, axis=1)

    df = df.sort_values(['chrom2size', 'start2'])
    df['idx2'] = range(0, df.shape[0])
    
    df = df.sort_values(['chrom1size', 'start1'])
    df['idx1'] = range(0, df.shape[0])
    
    # Row 0: forward
    # Row 1: reverse
    probs = np.zeros((2,df.shape[0]))
    # Idx of previous state
    prevs = np.zeros((2, df.shape[0]), dtype=int)

    # Do not penalize starting in forward or reverse
    chrom1 = df['chrom1'].to_list()
    start1 = df['start1'].to_list()
    end1 = df['end1'].to_list()

    chrom2 = df['chrom2'].to_list()
    start2 = df['start2'].to_list()
    end2 = df['end2'].to_list()

    idx2 = df['idx2'].to_list()
    
    print("Forward pass...", file=sys.stderr)

    for pos in tqdm(range(1, df.shape[0])):
    
        newchrom = False
        if chrom2[pos-1] != chrom2[pos]:
            newchrom = True
        #diff = idx2[pos] - idx2[pos-1]
        diff = start2[pos] - start2[pos-1]

        bp = -break_penalty
        if newchrom:
            bp += -translocation_penalty

        maxprevf = None
        maxprevind_f = None
        maxprevr = None
        maxprevind_r = None

        for prevind in range(pos-1, max(pos-1-max_skip-1, -1), -1):
            #skips = pos - prevind - 1
            skips = start1[pos] - start1[prevind] - 1

            thisdiff = start1[pos] - start1[prevind]
            fe = forward_emit(thisdiff, match_score)
            re = rev_emit(thisdiff, match_score)

            score_f = probs[0,prevind] + -skips*skip_penalty
            score_r = probs[1,prevind] + -skips*skip_penalty
            
            score_ff = score_f
            score_fr = score_r
            
            score_rf = score_f
            score_rr = score_r

            if chrom1[prevind] != chrom1[pos]:
                # Started a new chromosome
                # No penalty; no emission based on prev site
                pass
            else:
                if chrom2[prevind] != chrom2[pos]:
                    # Translocation
                    score_fr += -translocation_penalty
                    score_ff += -translocation_penalty
                    score_rf += -translocation_penalty
                    score_rr += -translocation_penalty
                else:
                    # Break
                    score_fr += -break_penalty
                    score_rf += -break_penalty

                    # No break
                    #score_ff += forward_emit(idx2[pos]-idx2[prevind], match_score)
                    #score_rr += rev_emit(idx2[pos]-idx2[prevind], match_score)
                    
                    score_ff += forward_emit(start2[pos]-start2[prevind], match_score)
                    score_rr += rev_emit(start2[pos]-start2[prevind], match_score)

            if score_ff > score_fr:
                if maxprevf is None or score_ff > maxprevf:
                    maxprevf = score_ff
                    maxprevind_f = prevind
            else:
                if maxprevf is None or score_fr > maxprevf:
                    maxprevf = score_fr
                    maxprevind_f = -prevind

            if score_rr > score_rf:
                if maxprevr is None or score_rr > maxprevr:
                    maxprevr = score_rr
                    maxprevind_r = -prevind
            else:
                if maxprevr is None or score_rf > maxprevr:
                    maxprevr = score_rf
                    maxprevind_r = prevind
        

        probs[0,pos] = maxprevf
        prevs[0,pos] = maxprevind_f
        probs[1,pos] = maxprevr
        prevs[1,pos] = maxprevind_r
        
        #if chrom1[pos] == 'chr5':
        #    print("{}\t{}\t{}\t{}\t{}".format(pos, maxprevf, maxprevind_f, maxprevr, maxprevind_r))

        #if pos > 29250:
        #    exit(1)


    print("Backtrace...", file=sys.stderr)
    
    # Get best ending position
    endp = probs.shape[1]-1
    endmod = 0
    endscore = probs[0,endp]
    if probs[1,endp] > endscore:
        endmod = 1
        endscore = probs[1,endp]
    endp_new = None
    
    for prevind in range(endp-1, max(endp-1-max_skip-1, -1), -1):
        skips = endp - prevind - 1
        score_f = probs[0,prevind] + -skips*skip_penalty
        score_r = probs[1,prevind] + -skips*skip_penalty
        if score_f > score_r:
            if score_f > endscore:
                endscore = score_f
                endmod = 0
                endp_new = prevind
        else:
            if score_r > endscore:
                endscore = score_r
                endmod = 1
                endp_new = prevind
    
    if endp_new is not None:
        endp = endp_new

    # Backtrace
    segs1 = []
    segs2 = []
    dirs = []

    seg1 = [chrom1[endp], start1[endp], end1[endp]]
    seg2 = [chrom2[endp], start2[endp], end2[endp]]
    seg_dir = '+'
    if endmod == 1:
        seg_dir = '-'

    while endp > 0:
        prevp = prevs[endmod,endp]

        if chrom1[abs(prevp)] != seg1[0] or chrom2[abs(prevp)] != seg2[0]:
            # Ended a chromosome.
            segs1.append(seg1)
            segs2.append(seg2)
            dirs.append(seg_dir)
            
            if prevp < 0:
                endmod = 1
                seg_dir = '-'
            else:
                endmod = 0
                seg_dir = '+'
            
            prevp = abs(prevp)
            seg1 = [chrom1[prevp], start1[prevp], end1[prevp]]
            seg2 = [chrom2[prevp], start2[prevp], end2[prevp]]
        
        elif (prevp < 0 and endmod == 0) or (prevp > 0 and endmod == 1):
            
            segs1.append(seg1)
            segs2.append(seg2)
            dirs.append(seg_dir)
            
            if prevp < 0:
                endmod = 1
                seg_dir = '-'
                prevp = -prevp
            else:
                endmod = 0
                seg_dir = "+"
            
            seg1 = [chrom1[prevp], start1[prevp], end1[prevp]]
            seg2 = [chrom2[prevp], start2[prevp], end2[prevp]]

        else:
            prevp = abs(prevp)

            # Extend
            if endmod == 0:
                seg1[1] = start1[prevp]
                seg2[1] = start2[prevp]
            else:
                seg1[1] = start1[prevp]
                seg2[2] = end2[prevp]
    
        endp = prevp

    segs1.append(seg1)
    segs2.append(seg2)
    dirs.append(seg_dir)

    segs1 = segs1[::-1]
    segs2 = segs2[::-1]
    dirs = dirs[::-1]
    
    return (segs1, segs2, dirs)

"""
After running the HMM on species1 coordinates, and then running again
on species2 coordinates, takes the results of both runs and finds a
consensus set of mappings of sequence blocks.
"""
def reconcile_dfs(segs1, segs1b, segs2, segs2b, dirs, dirsb):

    colnames = ['Chromosome', 'Start', 'End']
    s1a = pd.DataFrame(segs1, columns=colnames)
    s1a['Strand'] = dirs

    s1b = pd.DataFrame(segs1b, columns=colnames)
    s1b['Strand'] = dirsb
    
    s2a = pd.DataFrame(segs2, columns=colnames)
    s2a['Strand'] = dirs

    s2b = pd.DataFrame(segs2b, columns=colnames)
    s2b['Strand'] = dirsb
    
    s1a['id'] = range(0, s1a.shape[0])
    s1b['id'] = range(0, s1b.shape[0])
    s2a['id'] = range(0, s2a.shape[0])
    s2b['id'] = range(0, s2b.shape[0])
    
    gr1a = pr.PyRanges(s1a)
    gr1b = pr.PyRanges(s1b)
    
    gr2a = pr.PyRanges(s2a)
    gr2b = pr.PyRanges(s2b)
    
    gr_s1 = gr1a.join(gr1b, strandedness='same').df.rename({"Chromosome": "chrom1a", 
                                                            "Start": "start1a", 
                                                            "End": "end1a", 
                                                            "Strand": "orientation1a", 
                                                            "id": "ida", 
                                                            "Start_b": "start1b", 
                                                            "End_b": "end1b",
                                                           "Strand_b": "orientation1b", 
                                                            "id_b": "idb"}, axis=1)

    gr_s2 = gr2a.join(gr2b, strandedness='same').df.rename({"Chromosome": "chrom2a", 
                                                            "Start": "start2a", 
                                                            "End": "end2a", 
                                                            "Strand": "orientation2a", 
                                                            "id": "ida", 
                                                            "Start_b": "start2b", 
                                                            "End_b": "end2b",
                                                            "Strand_b": "orientation2b", 
                                                            "id_b": "idb"}, axis=1)
    
    gr_s1['jacc'] = gr_s1.apply(lambda x: 2*(min(x['end1a'], x['end1b']) \
            - max(x['start1a'], x['start1b']))/(x['end1a'] - x['start1a'] + x['end1b'] - x['start1b']), axis=1)
    gr_s2['jacc'] = gr_s2.apply(lambda x: 2*(min(x['end2a'], x['end2b']) \
            - max(x['start2a'], x['start2b']))/(x['end2a'] - x['start2a'] + x['end2b'] - x['start2b']), axis=1)
    
    grboth = gr_s1.merge(gr_s2, left_on=['ida', 'idb'], right_on=['ida', 'idb'])
    grboth['jsum'] = grboth['jacc_x'] + grboth['jacc_y']
    mj1 = grboth.groupby('ida', observed=True)['jsum'].agg('max').reset_index().rename({'jsum': 'maxjacca'}, axis=1)
    mj2 = grboth.groupby('idb', observed=True)['jsum'].agg('max').reset_index().rename({'jsum': 'maxjaccb'}, axis=1)
    grboth = grboth.merge(mj1).merge(mj2)
    grboth = grboth.loc[(grboth['jsum'] == grboth['maxjacca']) & (grboth['jsum'] == grboth['maxjaccb']),:]
    
    grboth['len1A'] = grboth['end1a'] - grboth['start1a']
    grboth['len1B'] = grboth['end1b'] - grboth['start1b']
    #grboth['len1C'] = grboth['end1a'] - grboth['start1b']
    #grboth['len1D'] = grboth['end1b'] - grboth['start1a']

    grboth['len2A'] = grboth['end2a'] - grboth['start2a']
    grboth['len2B'] = grboth['end2b'] - grboth['start2b']
    #grboth['len2C'] = grboth['end2a'] - grboth['start2b']
    #grboth['len2D'] = grboth['end2b'] - grboth['start2a']

    coords = grboth.apply(findmin, axis=1)
    starts1 = []
    ends1 = []
    starts2 = []
    ends2 = []
    for s1, e1, s2, e2 in coords:
        starts1.append(s1)
        ends1.append(e1)
        starts2.append(s2)
        ends2.append(e2)

    grboth['start1'] = starts1
    grboth['end1'] = ends1
    grboth['start2'] = starts2
    grboth['end2'] = ends2

    species1 = pd.DataFrame({'chrom1': grboth['chrom1a'],
        'start1': grboth['start1'],
        'end1': grboth['end1'],
        'orientation': grboth['orientation1a'],
        'chrom2': grboth['chrom2a'],
        'start2': grboth['start2'],
        'end2': grboth['end2']})
    
    return species1

"""
Check whether a file is gzipped.
"""
def is_gz(file):
    with open(file, 'rb') as test:
        return test.read(2) == b'\x1f\x8b'

def get_gz_lines(file, name_field, id_field, exclude_chr=None, include_chr=None, 
                names_all=None):
    f = None
    gz = False
    if is_gz(file):
        f = gzip.open(file)
        gz = True
    else:
        f = open(file)
    
    first = True
    gff3 = False
    
    excl_reg = None
    if exclude_chr is not None:
        excl_reg = re.compile(exclude_chr)
    incl_reg = None
    if include_chr is not None:
        incl_reg = re.compile(include_chr)

    for line in tqdm(f):
        if gz:
            line = line.decode().rstrip()
        else:
            line = line.rstrip()
        if line[0] != "#":
            dat = line.split('\t')
            if dat[2] == 'gene':
                if excl_reg is not None and excl_reg.search(dat[0]):
                    continue
                if incl_reg is not None and not incl_reg.search(dat[0]):
                    continue
                vals = dat[8]
                if first:
                    if '=' in vals:
                        gff3 = True
                    else:
                        gff3 = False
                    first = False
                gid = None
                gname = None
                for elt in vals.split(';'):
                    elt = elt.strip()
                    k = None
                    v = None
                    if gff3:
                        kv = elt.split('=')
                        if len(kv) == 2:
                            k, v = kv
                    else:
                        kv = elt.split()
                        if len(kv) == 2:
                            k, v = kv
                            v = v.strip('"')
                    if names_all is not None:
                        names_all.add(k)
                    if k == name_field:
                        gname = v
                    elif k == id_field:
                        gid = v
                
                if gid is not None and gname is not None:
                    yield {'chrom': dat[0], 'start': int(dat[3])-1, 
                           'end': int(dat[4]), 'id': gid, 'name': gname}
    f.close()

def get_gz_lines2(file, name_field, id_field, exclude_chr=None, include_chr=None, 
                names_all=None):
    f = None
    gz = False
    if is_gz(file):
        f = gzip.open(file)
        gz = True
    else:
        f = open(file)
    
    first = True
    gff3 = False
    
    excl_reg = None
    if exclude_chr is not None:
        excl_reg = re.compile(exclude_chr)
    incl_reg = None
    if include_chr is not None:
        incl_reg = re.compile(include_chr)

    for line in tqdm(f):
        if gz:
            line = line.decode().rstrip()
        else:
            line = line.rstrip()
        if line[0] != "#":
            dat = line.split('\t')
            if excl_reg is not None and excl_reg.search(dat[0]):
                continue
            if incl_reg is not None and not incl_reg.search(dat[0]):
                continue
            vals = dat[8]
            if first:
                if '=' in vals:
                    gff3 = True
                else:
                    gff3 = False
                first = False
            gid = None
            gname = None
            for elt in vals.split(';'):
                elt = elt.strip()
                k = None
                v = None
                if gff3:
                    kv = elt.split('=')
                    if len(kv) == 2:
                        k, v = kv
                else:
                    kv = elt.split()
                    if len(kv) == 2:
                        k, v = kv
                        v = v.strip('"')
                if names_all is not None:
                    names_all.add(k)
                if k == name_field:
                    gname = v
                elif k == id_field:
                    gid = v
            
            yield {'line': line, 'id': gid, 'name': gname, 'seq': dat[0]}
    f.close()

"""
Extract genes and coordinates from a GTF/GFF3, optionally gzipped.
"""
def get_genes(file, first_file=True, name_field="gene_name", id_field="gene_id", 
    exclude_chr=None, include_chr=None):
    if exclude_chr == "":
        exclude_chr = None
    if include_chr == "":
        include_chr = None
    
    gid2name = {}

    chroms = []
    starts = []
    ends = []
    genes = []

    if first_file:
        print("Load GTF/GFF3 file 1...", file=sys.stderr)
    else:
        print("Load GTF/GFF3 file 2...", file=sys.stderr)
    
    names_all = set([])
    for entry in get_gz_lines(file, name_field, id_field, exclude_chr, include_chr,
        names_all):
        chroms.append(entry['chrom'])
        starts.append(entry['start'])
        ends.append(entry['end'])
        genes.append(entry['name'])
        gid2name[entry['id']] = entry['name']

    if len(chroms) == 0:
        if first_file:
            print("ERROR: no genes found in --anno1/-1. Check --gene_name1/-n1", file=sys.stderr)
            exit(1)
        else:
            print("ERROR: no genes found in --anno2/-2. Check --gene_name2/-n2", file=sys.stderr)
        
        print("Candidate fields:", file=sys.stderr)
        for n in sorted(list(names_all)):
            print(n, file=sys.stderr)
        exit(1)

    dat = None
    if first_file:
        dat = {'chrom1': chroms, 'start1': starts, 'end1': ends, 'gene': genes}
    else:
        dat = {"chrom2": chroms, "start2": starts, "end2": ends, "gene": genes}

    return (pd.DataFrame(dat), gid2name)

"""
Extract chromosome/sequence lengths from a fasta index (fai).
"""
def get_chromlens(fai, excl_chr=None, incl_chr=None):
    if excl_chr == "":
        excl_chr = None
    if incl_chr == "":
        incl_chr = None
    excl_reg = None
    incl_reg = None
    if excl_chr is not None:
        excl_reg = re.compile(excl_chr)
    if incl_chr is not None:
        incl_reg = re.compile(incl_chr)

    f = open(fai, 'r')
    chroms = []
    lens = []
    for line in f:
        line = line.rstrip()
        dat = line.split('\t')
        if excl_reg is not None and excl_reg.search(dat[0]):
            continue
        if incl_reg is not None and not incl_reg.search(dat[0]):
            continue
        chroms.append(dat[0])
        lens.append(int(dat[1]))
    f.close()
    return (chroms, lens)

"""
Plot rectangles marking inversions/translocations before plotting gene coordinates.
"""
def add_rectangles(ax, df, col, maxX, maxY):
    for _, row in df.iterrows():
        rect = patches.Rectangle((row['chrstart1'] + row['start1'], 0),
                                 row['end1'] - row['start1'],
                                 maxY,
                                 linewidth=0,
                                 facecolor=col,
                                 alpha=1)
        ax.add_patch(rect)
        """
        rect = patches.Rectangle((row['chrstart1'] + row['start1'], row['chrstart2'] + row['start2']),
                                 row['end1'] - row['start1'],
                                 row['end2'] - row['start2'],
                                 linewidth=0,
                                 facecolor=col,
                                 alpha=1)

        ax.add_patch(rect)
        """
        r2 = patches.Rectangle((0, row['chrstart2'] + row['start2']),
                               maxX,
                               row['end2'] - row['start2'],
                               linewidth=0,
                               facecolor=col,
                               alpha=1)
        ax.add_patch(r2)

"""
From HMM output, match up chromosomes in species1 and species2, determining which are
homologous. For each chrom pair, the Jaccard index is computed as 
2*(intersecting bases)/(species1 bases + species2 bases).
Then, compute the maximum Jaccard index per species1 chromosome and the maximum per
species2 chromosome.
Chrom pairs that produce the maximum Jaccard index for the species1 chrom as well
as the species2 chrom are chosen as homologues.
"""
def get_homologous_chroms(species):
    species1 = species.copy()

    species1['size1'] = species1['end1'] - species1['start1']
    species1['size2'] = species1['end2'] - species1['start2']
    agg1 = species1.groupby(['chrom1', 'chrom2'], observed=True)['size1'].agg('sum').reset_index()
    agg2 = species1.groupby(['chrom1', 'chrom2'], observed=True)['size2'].agg('sum').reset_index()
    tot1 = species1.groupby('chrom1', observed=True)['size1'].agg('sum').reset_index().rename({'size1': 'tot1'}, axis=1)
    tot2 = species1.groupby('chrom2', observed=True)['size2'].agg('sum').reset_index().rename({'size2': 'tot2'}, axis=1)
    jacc = agg1.merge(agg2, left_on=['chrom1', 'chrom2'], right_on=['chrom1', 'chrom2'])
    jacc = jacc.merge(tot1, left_on='chrom1', right_on='chrom1')
    jacc = jacc.merge(tot2, left_on='chrom2', right_on='chrom2')
    jacc['jacc'] = (jacc['size1'] + jacc['size2'])/(jacc['tot1'] + jacc['tot2'])
    maxjacc1 = jacc.groupby('chrom1', observed=True)['jacc'].agg('max').reset_index().rename({'jacc': 'maxjacc1'}, axis=1)
    maxjacc2 = jacc.groupby('chrom2', observed=True)['jacc'].agg('max').reset_index().rename({'jacc': 'maxjacc2'}, axis=1)
    jacc = jacc.merge(maxjacc1, left_on='chrom1', right_on='chrom1')
    jacc = jacc.merge(maxjacc2, left_on='chrom2', right_on='chrom2')
    jacc = jacc.loc[(jacc['jacc'] == jacc['maxjacc1']) & (jacc['jacc'] == jacc['maxjacc2']),:]
    jacc = jacc.drop(['size1', 'size2', 'tot1', 'tot2', 'jacc', 'maxjacc1', 'maxjacc2'], axis=1)
    jacc['homologous'] = 1
    return jacc

"""
Used as a pd.DataFrame.apply() function by the reconcile_dfs() function.
After running the HMM of species2 coords on species1 coords, then also
running it on species1 coords on species2 coords, need to reconcile the
output of both.
We find all sets of (species1a, species2a), (species1b, species2b) 
intervals, where a denotes the first run and b denotes the second run,
and species1a/species1b intersect, and species2a/species2b also intersect.
We then need to choose a single (species1, species2) interval pair as
consensus.

For this, we choose the (species1, species2) interval where the length
on species1 is the closest to the length on species2.
"""
def findmin(x):
    keys = {'len1A': ('start1a', 'end1a'),
            'len1B': ('start1b', 'end1b'),
            'len1C': ('start1a', 'end1b'),
            'len1D': ('start1b', 'end1a'),
            'len2A': ('start2a', 'end2a'),
            'len2B': ('start2b', 'end2b'),
            'len2C': ('start2a', 'end2b'),
            'len2D': ('start2b', 'end2a')}

    key1 = ['len1A', 'len1B', 'len1C', 'len1D']
    key2 = ['len2A', 'len2B', 'len2C', 'len2C']
    
    min1 = None
    min2 = None
    minval = None

    for i in range(0, len(key1)):
        k1 = key1[i]
        for j in range(0, len(key2)):
            k2 = key2[j]
            if k1 in x and k2 in x:
                diff = abs(x[k1] - x[k2])
                if minval is None or diff < minval:
                    minval = diff
                    min1 = k1
                    min2 = k2
    
    return [x[keys[min1][0]], x[keys[min1][1]], x[keys[min2][0]], x[keys[min2][1]]]

"""
Create plot at the end.
"""
def make_plot(species1, df, chroms1, lens1, chroms2, lens2, 
    species1name, species2name, output):
    
    # Make plottable data.
    l1df = pd.DataFrame({'chrom1': chroms1, 'len1': lens1})
    missing1 = set(df['chrom1'].to_list()).difference(set(l1df['chrom1'].to_list()))
    if len(missing1) > 0:
        print("ERROR: chromosomes in --anno1 missing in --fai1:", file=sys.stderr())
        for chrom in sorted(list(missing1)):
            print(chrom, file=sys.stderr)
        print("Check/adjust files, or adjust --exclude_chr/--include_chr", file=sys.stderr)
        exit(1)
    l1df = l1df.loc[l1df['chrom1'].isin(set(df['chrom1'].to_list())),:]

    l2df = pd.DataFrame({'chrom2': chroms2, 'len2': lens2})
    missing2 = set(df['chrom2'].to_list()).difference(set(l2df['chrom2'].to_list()))
    if len(missing2) > 0:
        print("ERROR: chromosomes in --anno2 missing in --fai2:", file=sys.stderr)
        for chrom in sorted(list(missing2)):
            print(chrom, file=sys.stderr)
        print("Check/adjust files, or adjust --exclude_chr/--include_chr", file=sys.stderr)
        exit(1)
    
    # Learn which chromosomes are homologous, via reciprocal maximum Jaccard index
    jacc = get_homologous_chroms(species1) 

    l2df = l2df.loc[l2df['chrom2'].isin(set(df['chrom2'].to_list())),:]
    l1df = l1df.sort_values('len1', ascending=False)
    l2df = l2df.sort_values('len2', ascending=False)
    
    # Introduce homologous chroms to l2df
    l2df = l2df.merge(jacc, left_on='chrom2', right_on='chrom2', how='left')
    nonhom = l2df.loc[l2df['chrom1'].isna(),:]
    l2df = l2df.loc[~l2df['chrom1'].isna(),:]

    # Sort l2df according to chrom1 order in l1df
    l2df = pd.concat([l1df.drop(['len1'], axis=1).merge(l2df, left_on='chrom1', right_on='chrom1', 
                                             sort=False).drop(['chrom1', 'homologous'], axis=1),
                      nonhom.drop(['chrom1', 'homologous'], axis=1)], axis=0)

    l1df['chrstart1'] = l1df['len1'].cumsum().shift(1, fill_value=0)
    l2df['chrstart2'] = l2df['len2'].cumsum().shift(1, fill_value=0)
    
    df = df.merge(l1df, left_on='chrom1', right_on='chrom1').merge(l2df, 
                                                                    left_on='chrom2', 
                                                                   right_on='chrom2')
    
    l1df['even'] = range(0, l1df.shape[0])
    l1df['even'] = l1df['even'] % 2
    l2df['even'] = range(0, l2df.shape[0])
    l2df['even'] = l2df['even'] % 2

    col_line = "#0F0E0E"
    col_line_trans = "#522B47"

    col_bg_inv = "#8DAA9D"
    col_bg_trans = "#CB807D"
    col_bg_invtrans = "#9C7EA0"

    maxX = (l1df['chrstart1'] + l1df['len1']).max()
    maxY = (l2df['chrstart2'] + l2df['len2']).max()

    fig, ax = plt.subplots(figsize=(10,10))
    
    # Disable borders
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Add colored rectangles
    species1 = species1.merge(l1df, left_on='chrom1', right_on='chrom1')
    species1 = species1.merge(l2df, left_on='chrom2', right_on='chrom2')
    
    species1 = species1.merge(jacc, left_on=['chrom1', 'chrom2'], 
                              right_on=['chrom1', 'chrom2'], how='left')
    species1.loc[species1['homologous'].isna(),'homologous'] = 0
    
    invs = species1.loc[(species1['homologous'] == 1) & (species1['orientation'] == '-'),:]
    trans = species1.loc[(species1['homologous'] == 0) & (species1['orientation'] == '+'),:]
    invtrans = species1.loc[(species1['homologous'] == 0) & (species1['orientation'] == '-'),:]
    
    add_rectangles(ax, invs, col_bg_inv, maxX, maxY)
    add_rectangles(ax, trans, col_bg_trans, maxX, maxY)
    add_rectangles(ax, invtrans, col_bg_invtrans, maxX, maxY)

    # Add vertical lines
    ax.vlines(l1df['chrstart1'], 0, maxY, lw=0.5, colors='black', linestyle='dotted')
    max1 = l1df.loc[l1df['chrstart1'] == l1df['chrstart1'].max(),:]
    ax.vlines(max1['chrstart1'] + max1['len1'], 0, maxY, lw=0.5, colors='black', linestyle='dotted')
    ax.hlines(l2df['chrstart2'], 0, maxX, lw=0.5, colors='black', linestyle='dotted')
    max2 = l2df.loc[l2df['chrstart2'] == l2df['chrstart2'].max(),:]
    ax.hlines(max2['chrstart2'] + max2['len2'], 0, maxX, lw=0.5, colors='black', linestyle='dotted')
    
    # Add chrom names
    text_bump_y = 125e6
    text_bump_x = 100e6

    for _, row in l1df.iterrows():
        ax.text(row['chrstart1'] + row['len1']/2, 
                -text_bump_x + 0.5*text_bump_x*row['even'],
                row['chrom1'], fontsize=5, ha='center')
    for _, row in l2df.iterrows():
        ax.text(-text_bump_y + 0.5*text_bump_y*row['even'],
            row['chrstart2'] + row['len2']/2,
            row['chrom2'], fontsize=5, rotation=90, va='center')

    # Plot genes (as points)
    ax.scatter(df['chrstart1'] + (df['start1'] + df['end1'])/2,
        df['chrstart2'] + (df['start2'] + df['end2'])/2,
        color='#000000',
        marker=",",
        linewidths=0,
        s=1)

    ax.set_xlabel(species1name)
    ax.set_ylabel(species2name)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    
    plt.savefig("{}.pdf".format(output), format='pdf')
    
    return (l1df, l2df)

"""
Save lists of genes as text files.
"""
def write_genelists(genes, svs, output_prefix):
    jacc = get_homologous_chroms(svs)
    
    gcpy = pr.PyRanges(genes.copy().drop(['chrom2', 'start2', 'end2'], axis=1)\
            .rename({"chrom1": "Chromosome", "start1": "Start", "end1": "End"}, axis=1))

    svcpy = svs.copy().merge(jacc, how='left').drop(['chrom2', 'start2', 'end2'], axis=1).\
            rename({"chrom1": "Chromosome", "start1": "Start", "end1": "End"}, axis=1)
    
    invs = pr.PyRanges(svcpy.loc[(svcpy['homologous'] == 1) & (svcpy['orientation'] == '-'),:]\
            .drop(['orientation'], axis=1).rename({"chrom1": "Chromosome", 
                                                   "start1": "Start", "end1": "End"}))
    tx = pr.PyRanges(svcpy.loc[svcpy['homologous'] == 0,:].drop(['orientation'], axis=1))
    
    g_inv = gcpy.join(invs).df
    g_tx = gcpy.join(tx).df
    
    if g_inv.shape[0] > 0:
        #df = pd.DataFrame({'gene': g_inv['gene']})
        df = g_inv.drop(['homologous'], axis=1)
        df = df.sort_values(['gene'])
        df.to_csv("{}.inv_genes.txt".format(output_prefix), sep='\t', index=False, header=False)
    if g_tx.shape[0] > 0:
        #df = pd.DataFrame({'gene': g_tx['gene']})
        df = g_tx.drop(['homologous'], axis=1)
        df = df.sort_values(['gene'])
        df.to_csv("{}.trans_genes.txt".format(output_prefix), sep='\t', index=False, header=False)
    
    # Write out set of all genes (to use as background in enrichment tests)
    #df = pd.DataFrame({'gene': genes['gene']}).sort_values(['gene'])
    df = genes.sort_values(['chrom1', 'start1', 'end1'])
    df.to_csv('{}.all_genes.txt'.format(output_prefix), sep='\t', index=False, header=False)

def samp_region(genes_gr, l1df, size, totsize):
    success = False
    while not success:
        base = int(round(random.random() * totsize))
        endbase = base + size
        chrom = l1df.loc[(l1df['chrstart1'] <= base) & (l1df['chrstart1'] + l1df['len1'] > endbase),:]
        if chrom.shape[0] == 0:
            continue
        else:
            base_chrom = base - chrom['chrstart1'].to_list()[0]
            endbase_chrom = base_chrom + size
            gr = pr.PyRanges(pd.DataFrame({'Chromosome': chrom['chrom1'], 
                              'Start': base_chrom, 'End': endbase_chrom}, index=[0]))
            inter = gr.join(genes_gr).df
            success = True
            if inter.shape[0] > 0:
                return len(inter['gene'].unique())
            else:
                return 0

def enrich_test(species1, l1df, l2df, df, nsamp):
    invs = species1.loc[(species1['orientation'] == '-'),:]
    genes_gr = pr.PyRanges(pd.DataFrame({'Chromosome': df['chrom1'],
                                         'Start': df['start1'],
                                         'End': df['end1'],
                                         'gene': df['gene']}))
    
    totsize = l1df['chrstart1'][l1df.shape[0]-1] + l1df['len1'][l1df.shape[0]-1]
    
    ngenes = []
    mean_genes_trials = []
    enrich_genes_p = []
    for _, row in invs.iterrows():
        size = row['end1'] - row['start1']
        gr1 = pr.PyRanges(pd.DataFrame({'Chromosome': row['chrom1'],
                                        'Start': row['start1'],
                                        'End': row['end1']}, index=[0]))
        ngenes_true = len(gr1.join(genes_gr).df['gene'].unique())
        ngenes.append(ngenes_true)
        nmore = 0
        trials = nsamp
        nums = []
        for i in tqdm(range(0, trials)):
            ngenes_trial = samp_region(genes_gr, l1df, size, totsize)
            nums.append(ngenes_trial)
            if ngenes_trial >= ngenes_true:
                nmore += 1
        mean_genes_trials.append(np.mean(np.array(nums)))
        enrich_genes_p.append(nmore/trials)
    
    invs.loc[:,'ngenes'] = ngenes_true
    invs.loc[:,'mean_ngenes_trials'] = mean_genes_trials
    invs.loc[:,'enrich_genes_p'] = enrich_genes_p
    
    return invs

def filter_annotations(outbase,
                       species1, 
                       species2, 
                       gene_df, 
                       annofile1, 
                       gid1, 
                       annofile2, 
                       gid2,
                       name_field1,
                       id_field1,
                       name_field2,
                       id_field2,
                       exclude_chr,
                       include_chr,
                       mito1,
                       mito2,
                       filt_files=False):
    
    g1 = pr.PyRanges(gene_df.copy().drop(['chrom2', 'start2', 'end2'], axis=1)\
        .rename({"chrom1": "Chromosome", "start1": "Start", "end1": "End"}, axis=1))
    g2 = pr.PyRanges(gene_df.copy().drop(['chrom1', 'start1', 'end1'], axis=1)\
        .rename({"chrom2": "Chromosome", "start2": "Start", "end2": "End"}, axis=1))
    
    species1['id'] = range(0, species1.shape[0])
        
    s1 = pr.PyRanges(species1.copy().drop(['chrom2', 'start2', 'end2', 'orientation'], axis=1)\
            .rename({"chrom1": "Chromosome", "start1": "Start", "end1": "End"}, axis=1))
    s2 = pr.PyRanges(species1.copy().drop(['chrom1', 'start1', 'end1', 'orientation'], axis=1)\
            .rename({"chrom2": "Chromosome", "start2": "Start", "end2": "End"}, axis=1))
    
    gs1 = g1.join(s1).df.drop(['Chromosome', 'Start', 'End', 'Start_b', 'End_b'], axis=1)
    gs2 = g2.join(s2).df.drop(['Chromosome', 'Start', 'End', 'Start_b', 'End_b'], axis=1)
    
    gsmerge = gs1.merge(gs2, left_on='gene', right_on='gene')
    gsmerge = gsmerge.merge(species1.drop(['start1', 'end1', 
                                           'start2', 'end2', 'orientation'],
                                          axis=1), left_on='id_x', right_on='id').drop(
                                                  ['id'], axis=1).rename({
                                                          'chrom1': "chr1x",
                                                          'chrom2': "chr2x"
                                                          }, axis=1)
    gsmerge = gsmerge.merge(species1.drop(['start2', 'end2',
                                           'start1', 'end1', 'orientation'],
                                          axis=1), left_on='id_y', right_on='id').drop(
                                                  ['id'], axis=1).rename({
                                                      'chrom1': "chr1y",
                                                      'chrom2': "chr2y"
                                                      }, axis=1)

    # Strict
    gkeep = set(gsmerge.loc[(gsmerge['id_x'] == gsmerge['id_y']),'gene'].unique())
    # Loose
    #gkeep = set(gsmerge.loc[(gsmerge['id_x'] == gsmerge['id_y']) | 
    #                        ((gsmerge['chr1x'] == gsmerge['chr1y']) &
    #                          (gsmerge['chr2x'] == gsmerge['chr2y'])),'gene'].unique())
    
    # Map discarded genes to the ID of region on species 1 to which we expect them to belong
    # (if treating species 1 as truth)
    g2seg1 = {}
    # Track all synteny tracks that are meant to contain genes from species2 that are getting
    # discarded
    segs1withgenes = set([])
    # Output file for GTF entries from species2 expected to fall in synteny tracks above
    segs1_gtf_out = {}
    if filt_files:
        for _, row in gsmerge.loc[~gsmerge['gene'].isin(gkeep),:].iterrows():
            g2seg1[row['gene']] = row['id_x']
            segs1withgenes.add(row['id_x'])
        for seg_id in segs1withgenes:
            outn1 = '{}.remap.segment.{}.bed'.format(outbase, seg_id)
            outn2 = '{}.remap.segment.{}.gtf'.format(outbase, seg_id)
            segs1_gtf_out[seg_id] = open(outn2, 'w')
            # Write out BED segment
            bed_seg = species1.loc[species1['id'] == seg_id,:].drop(['chrom1', 'start1', 'end1', 'orientation', 'id'], axis=1)
            # Write out species2 expected segment in BED format
            bed_seg.to_csv(outn1, sep='\t', header=False, index=False)

    g_all = set(gsmerge['gene'].unique())
    g_rm = g_all.difference(gkeep)
    print("Remove genes ({}):".format(len(g_rm)), file=sys.stderr)
    for g in sorted(list(g_rm)):
        print(g, file=sys.stderr)
    
    annobase1 = annofile1.split('/')[-1]
    annoext1 = ""
    m = re.match(r'(.+)\.(gff|gff3|gtf)(.gz)?', annobase1)
    if m:
        annobase1 = m.group(1)
        annoext1 = m.group(2)
    annobase2 = annofile2.split('/')[-1]
    annoext2 = ""
    m = re.match(r'(.+)\.(gff|gff3|gtf)(.gz)?', annobase2)
    if m:
        annobase2 = m.group(1)
        annoext2 = m.group(2)
    f1 = gzip.open('{}.{}.{}.gz'.format(outbase, annobase1, annoext1), 'wt')
    f1r = gzip.open('{}.removed.{}.{}.gz'.format(outbase, annobase1, annoext1), 'wt')
    f1rl = open('{}.removed.{}.genes'.format(outbase, annobase1), 'w')
    print("Filter annotation 1...", file=sys.stderr)
    f1rm = set([])
    for entry in get_gz_lines2(annofile1, name_field1, id_field1, exclude_chr,
                              include_chr):
        if entry['seq'] == mito1 or (
            entry['id'] is not None and entry['id'] in gid1 and \
            gid1[entry['id']] in gkeep ):
                print(entry['line'], file=f1)
        else:
            print(entry['line'], file=f1r)
            if entry['name'] is not None:
                f1rm.add(entry['name'])
    f1.close()
    f1r.close()
    for rm in sorted(list(f1rm)):
        print(rm, file=f1rl)
    f1rl.close()

    print("Filter annotation 2...", file=sys.stderr)
    f2 = gzip.open('{}.{}.{}.gz'.format(outbase, annobase2, annoext2), 'wt')
    f2r = gzip.open('{}.removed.{}.{}.gz'.format(outbase, annobase2, annoext2), 'wt')
    f2rl = open('{}.removed.{}.genes'.format(outbase, annobase2), 'w')
    f2rm = set([])
    for entry in get_gz_lines2(annofile2, name_field2, id_field2, exclude_chr,
                               include_chr):
        if entry['seq'] == mito2 or (
            entry['id'] is not None and entry['id'] in gid2 and \
            gid2[entry['id']] in gkeep ):
                print(entry['line'], file=f2)
        else:
            print(entry['line'], file=f2r)
            if entry['name'] is not None:
                f2rm.add(entry['name'])
            if filt_files and entry['name'] in g2seg1:
                # Find which output file this belongs to
                seg = g2seg1[entry['name']]
                print(entry['line'], file=segs1_gtf_out[seg])
            
    if filt_files:
        for seg in segs1_gtf_out:
            segs1_gtf_out[seg].close()
    for rm in sorted(list(f2rm)):
        print(rm, file=f2rl)

    f2.close()
    f2r.close()
    f2rl.close()

def main(args):
    options = parse_args()
    
    if options.exclude_chr is not None:
        if options.exclude_chr[0] != '(':
            options.exclude_chr = '(' + options.exclude_chr
        if options.exclude_chr[-1] != ')':
            options.exclude_chr += ')'
    
    # Parse annotation files
    df1, gid1 = get_genes(options.anno1, first_file=True, name_field=options.gene_name1,
            id_field=options.gene_id1, exclude_chr=options.exclude_chr, 
            include_chr=options.include_chr)
    df2, gid2 = get_genes(options.anno2, first_file=False, name_field=options.gene_name2,
            id_field=options.gene_id2, exclude_chr=options.exclude_chr, 
            include_chr=options.include_chr)
    
    # Merge together to compare position by gene
    df = df1.merge(df2, left_on='gene', right_on='gene')
    
    # Find path through species1-centric HMM
    segs1, segs2, dirs = find_path(df, options.match_score, options.skip_penalty,
            options.break_penalty, options.translocation_penalty, options.max_skip)

    df2 = df.copy().rename({'chrom1': 'chrom2', 'start1': 'start2', 'end1': 'end2', 'chrom2': 'chrom1',
        'start2': 'start1', 'end2': 'end1'}, axis=1)
    
    # Find path through species2-centric HMM
    segs2b, segs1b, dirsb = find_path(df2, options.match_score, options.skip_penalty,
            options.break_penalty, options.translocation_penalty, options.max_skip)
    
    # Take both paths and find mappings of (genome1, genome2) common to both
    # Return a set of (genome1, genome2) segments with relative orientation
    species1 = reconcile_dfs(segs1, segs1b, segs2, segs2b, dirs, dirsb) 

    # Flip column order to create second output file (write in species2 coordinates 
    # instead of species1 coordinates)
    species2 = pd.DataFrame({'chrom2': species1['chrom2'], 'start2': species1['start2'], 
                           'end2': species1['end2'], 'orientation': species1['orientation'],
                           'chrom1': species1['chrom1'], 'start1': species1['start1'],
                           'end1': species1['end1']})

    species1.to_csv('{}.{}.bed'.format(options.output, "_".join(options.species1.split(' '))), 
        header=False, index=False, sep='\t')
    species2.to_csv('{}.{}.bed'.format(options.output, "_".join(options.species2.split(' '))), 
        header=False, index=False, sep='\t')
    
    # Create lists of genes
    write_genelists(df, species1, options.output)

    # Check whether chrom lengths were provided, which allows plotting
    if options.fai1 is not None and options.fai2 is not None:
        print("Make plot...", file=sys.stderr)
        chroms1, lens1 = get_chromlens(options.fai1, options.exclude_chr, options.include_chr)
        chroms2, lens2 = get_chromlens(options.fai2, options.exclude_chr, options.include_chr)
        l1df, l2df = make_plot(species1, df, chroms1, lens1, chroms2, lens2, 
            options.species1, options.species2, options.output)
        
        if not options.nofilt and options.enrichment is not None and options.enrichment > 0:
            print("Enrichment testing...", file=sys.stderr)
            inv_enrich, trans_enrich = enrich_test(species1, l1df, l2df, 
                                                   df, options.enrichment)
    if not options.nofilt:
        # Filter annotations to exclude genes mapped to wrong regions
        print("Filtering annotations...", file=sys.stderr)
        filter_annotations(options.output, species1, species2, df, 
                           options.anno1, gid1, options.anno2, gid2,
                           options.gene_name1, options.gene_id1,
                           options.gene_name2, options.gene_id2,
                           options.exclude_chr, options.include_chr,
                           options.mito1, options.mito2,
                           options.filt_files)

if __name__ == '__main__':
    sys.exit(main(sys.argv))
