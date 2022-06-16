#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import collections
import operator
import pandas as pd
import string

##################################seq_id_fetcher##########################################
def seq_id_fetcher (IN1,IN2,IN3,OUT,LOG):

    corRdogtaglist = []
    corRdogtag2line = {}
    corFdict = collections.defaultdict(list)
    corRdogtagdict = collections.defaultdict(list)
    corRdogtag2real = {}
    newlines = []

    with open(IN1, 'r') as file2in1:
	    header = file2in1.readline()
	    header = header.rstrip('\n')
	    for linein in file2in1:
		    linein = linein.rstrip('\n')
		    eles = re.split(r'\t', linein)
		    readcount = eles[0]
		    corF = eles[1] + eles[2] + eles[3]
		    corR = eles[1] + eles[2] + eles[3] + '&' + eles[7]
		    dogtag = eles[8]
		    corFdogtag = corF + '&' + dogtag
		    corRdogtag = corR + '&' + dogtag
		    corRdogtaglist.append(corRdogtag)
		    corRdogtag2line[corRdogtag] = linein
		    corFdict[corF].append(linein)
		
    with open(IN2, 'r') as file2in2:
	    for linein in file2in2:
		    linein = linein.rstrip('\n')
		    eles = re.split(r'\t', linein)
		    seqid = eles[0]
		    corF = eles[1] + eles[2] + eles[3]
		    corF2 = corF + '*'
		    if eles[6]:
		        length = abs(int(eles[6])-int(eles[3]))
		    else:
		        length = ''
		    dogtag = eles[7]
		    corRdogtag = eles[1] + eles[2] + eles[3] + '&' + str(length) + '&' + eles[7]
		    corRdogtagdict[corRdogtag].append(seqid)

    with open(IN3, 'r') as file2in3:
	    for linein in file2in3:
		    linein = linein.rstrip('\n')
		    eles = re.split(r'\t', linein)
		    corR = eles[0] + eles[1] + eles[2] + '&' + eles[6]
		    real_dogtag = eles[7]
		    dogtag = eles[11]
		    corRdogtag = corR + '&' + dogtag
		    corRdogtag2real[corRdogtag] = real_dogtag
		
    for corRdogtag in corRdogtaglist:
        linein = corRdogtag2line[corRdogtag]
        corRdogtag_eles = re.split(r'&', corRdogtag)
        corF = corRdogtag_eles[0]
        corR = corRdogtag_eles[0] + '&' + corRdogtag_eles[1]
        BPcount = len(corFdict[corF])
        id2count = {}
        real_dogtag = corRdogtag2real[corRdogtag]
        corRdogtag_real = corR + '&' + real_dogtag
        for seqid_all in corRdogtagdict[corRdogtag_real]:
            seqid, readcount = seqid_all.split('#',2)
            id2count[seqid] = readcount
        major_id = max(id2count.items(), key=operator.itemgetter(1))[0]
        newline = linein + '\t' + str(BPcount) + '\t' + major_id
        newlines.append(newline)
     
    df = pd.DataFrame([l.split('\t') for l in newlines])
    df[0] = df[0].astype(int)
    df[7] = df[7].astype(int)
    df[9] = df[9].astype(int)
    df_sorted = df.sort_values(by=[9,1,2,3,7,0], ascending = [False, True, True, True, False, False])
    df_sorted = df_sorted.astype(str)
    newlineslist = df_sorted.values.tolist()

    with open(OUT, 'w') as file2out:
        file2out.write(header + '\tbreakpoint_count\texample_seq_ID\n')
        for newlinelist in newlineslist:
            newline = '\t'.join([str(ele) for ele in newlinelist])
            file2out.write(newline + '\n')
        
    with open(LOG, 'a') as file2log:
        file2log.write('sequence id fetching: done\n')
##########################################################################################



##################################ref_seq_fetcher#########################################
def ref_seq_fetcher (IN,OUT,LOG):

    def getseq_fast(chrom, seqstart, seqend):
        refpath = '/Users/guos2/Genomes/hg19/hg19_by_chroms/'
        with open(refpath + chrom + '.fa', 'r') as ref:
            if seqend < seqstart:
                raise Exception('Incorrect coordinates!\n')
            header = ref.readline()
            header = header.rstrip('\n')
            if not re.match(r'^>', header):
                raise Exception('The reference file is not in fasta format!\n')
        
            headeroffset = len(header)
            lineref = ref.readline()
            basesperline = len(lineref) - 1

            # compute various offsets required...
            lfstostart = int((seqstart/basesperline)) #how many lfeed bytes
            lfswithin = int(seqend/basesperline) - lfstostart #line feed bytes within sequence
            startbyte = headeroffset + lfstostart + seqstart #fasta header + linefeeder + start
            bytestoread = seqend + lfswithin - seqstart
        
            # go get the sequence...
            ref.seek(startbyte, 0)
            seq = ref.read(bytestoread)
            seq = re.sub('\s', '', seq)
            return seq

    file2out = open (OUT, 'w')
    with open(IN, 'r') as file2in:
	    header = file2in.readline()
	    header = header.rstrip('\n')
	    file2out.write(header + '\tupstream_30bp\tdownstream_30bp\n')
	    for linein in file2in:
		    linein = linein.rstrip('\n')
		    eles = re.split(r'\t', linein)
		    chrom = eles[1]
		    strand = eles[2]
		    cor = int(eles[3])
		    if 'chr' in chrom and 'GL' not in chrom:
		        if strand == '+':
		            seqstart = cor - 30
		            seqend = cor + 30 + 1
		            seq = getseq_fast(chrom, seqstart, seqend)
		            seq = seq.upper()
		        elif strand == '-':
		            seqstart = cor - 30 + 1
		            seqend = cor + 30 + 2
		            seq = getseq_fast(chrom, seqstart, seqend)
		            seq = seq.upper()
		            seq = seq[::-1]
		            seq = seq.translate(str.maketrans({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}))
		        sequp = seq[1:31]
		        seqdown = seq[31:]
		        file2out.write(linein + '\t' + sequp + '\t' + seqdown + '\n')

    file2out.close()

    with open(LOG, 'a') as file2log:
        file2log.write('reference sequence fetching: done\n')
##########################################################################################



######################################seq_fetcher#########################################
def seq_fetcher (IN1,IN2,OUT,LOG):

    id2line = {}

    with open(IN1, 'r') as file2in1:
	    for linein in file2in1:
		    linein = linein.rstrip('\n')
		    eles = re.split(r'\t', linein)
		    seqid, rest = eles[0].split('&',2)
		    newline = eles[1] + '\t' + eles[2] + '\t' + eles[4] + '\t' + eles[5]
		    id2line[seqid] = newline

    file2out = open(OUT, 'w')
    with open(IN2, 'r') as file2in2:
	    header = file2in2.readline()
	    header = header.rstrip('\n')
	    file2out.write(header + '\tread1_index\tread1_trimmed\tread2_index\tread2_trimmed\tspike_in_control\n')
	    for linein in file2in2:
		    linein = linein.rstrip('\n')
		    eles = re.split(r'\t', linein)
		    cor = eles[1] + eles[2] + eles[3]
		    seqid = eles[10]
		    seq = id2line[seqid]
		    mark = 'pass'
		    if cor == 'chr17+3618824':
		        mark = 'spike_in'
		    elif eles[2] == eles[5]:
		        mark = 'strand_error'
		    file2out.write(linein + '\t' + seq + '\t' + mark + '\n')
    file2out.close()

    with open(LOG, 'a') as file2log:
        file2log.write('read sequence fetching: done\n')
##########################################################################################



######################################misprime_inspector##################################
def misprime_inspector (primer_seq,IN,OUT,LOG):
    loci2score1 = {}
    loci2score2 = {}

    file2out = open(OUT, 'w')
    with open(IN, 'r') as file2in:
        header = file2in.readline()
        header = header.rstrip('\n')
        file2out.write(header + '\tmisprime\n')
        for linein in file2in:
            linein = linein.rstrip('\n')
            eles = re.split(r'\t', linein)
            eles.append('pass')
            primer_len = len(primer_seq)
            upstream_seq = eles[11]
            upstream_len = len(upstream_seq)
            len_diff = abs(upstream_len - primer_len)
            genomic_seq = eles[13]
            loci1 = loci2 = score1 = score2 = seq_match = max1 = max2 = 0
            constant_primer_seq = primer_seq
            constant_upstream_seq = upstream_seq
            constant_genomic_seq = genomic_seq
            constant_primer_len = primer_len
            constant_upstream_len = upstream_len
        
            # move upstream seq around primer seq, loci=0(left side aligned), truncate 1 nucleotide on upstream seq after each alignment, until the right side is aligned    
            while not loci1 > len_diff:
                len_range = min(upstream_len, primer_len)
                for i in range(0, len_range-1): # at each loci, match the nucleotides from first to last nucleotide in range
                    if upstream_seq[i] == primer_seq[i]:
                        score1 += 1 # record the matched nucleotide number
                loci2score1[loci1] = str(loci1) + '#' + str(score1) # assign score to each loci
                loci1 += 1
                upstream_seq = constant_upstream_seq[loci1:] # truncate upstream seq after each loci scan
                upstream_len = len(upstream_seq)
                score1 = 0
        
            # reset upstream seq and length for next cycle
            upstream_seq = constant_upstream_seq
            upstream_len = constant_upstream_len
        
            # move primer seq around upstream seq, loci=0(left side aligned), truncate 1 nucleotide on upstream seq after each alignment, until there is no way to get a higher score
            while loci2 < 8: # primer length(12) * lowest percent mismatch(70%) = 8.4
                for i in range(0, primer_len-1):
                    if upstream_seq[i] == primer_seq[i]:
                        score2 += 1 # record the matched nucleotide number
                loci2score2[loci2] = str(loci2) + '#' + str(score2) # assign score to each loci
                loci2 += 1
                primer_seq = constant_primer_seq[loci2:] # truncate primer seq after each loci scan
                primer_len = len(primer_seq)
                score2 = 0

            # reset primer seq and length for next cycle
            primer_seq = constant_primer_seq
            primer_len = constant_primer_len
        
            # get both of the best loci for previous two ways of scanning   
            score1s = []
            score2s = []
            for score1 in loci2score1.values():
                rest, real_score1 = score1.split('#',2)
                score1s.append(int(real_score1))
            max1 = max(score1s)
            for score2 in loci2score2.values():
                rest, real_score2 = score2.split('#',2)
                score2s.append(int(real_score2))
            max1 = max(score2s)

            # when max1 >= max2, at each best loci, scan the tail seqs
            best_loci = []
            if not max1 < max2:
                score1_str = '#' + str(max1)
                for loci, score1 in loci2score1.items():
                    if score1_str in score1:
                        best_loci.append(loci)
                for i in best_loci:
                    upstream_tail = constant_upstream_seq[i+constant_primer_len:]
                    genomic_head, genomic_tail = genomic_seq.split(primer_seq,2)
                    for i in range(0, min(len(upstream_tail), len(genomic_tail))): # scan upstream tail and genomic tail at the best loci
                        if upstream_tail[i] == genomic_tail[i]:
                            seq_match += 1
                    if seq_match > 10 and max1 > 8:
                        eles[-1] = 'misprime!(upstream_loci:' + str(best_loci) + ',tail_match:' + str(seq_match) + ',primer_mismatch:' + str(max1) + ')' # match threshold is adjustable
                    seq_match = 0
                
            # when max1 <= max2, at each loci, scan the tail seqs
            best_loci = []
            if not max1 > max2:
                score2_str = '#' + str(max2)
                for loci, score2 in loci2score2.items():
                    if score2_str in score2:
                        best_loci.append(loci)
                for i in best_loci:
                    upstream_tail = constant_upstream_seq[constant_primer_len-i:]
                    genomic_head, genomic_tail = genomic_seq.split(primer_seq,2)
                    for i in range(0, min(len(upstream_tail), len(genomic_tail))): # scan upstream tail and genomic tail at the best loci
                        if upstream_tail[i] == genomic_tail[i]:
                            seq_match += 1
                    if seq_match > 10 and max2 > 8:
                        eles[-1] = 'misprime!(primer_loci:' + str(best_loci) + ',tail_match:' + str(seq_match) + ',primer_mismatch:' + str(max2) + ')' # match threshold is adjustable
                    seq_match = 0

            newline = '\t'.join(eles)
            file2out.write(newline + '\n')        
    file2out.close()

    with open(LOG, 'a') as file2log:
        file2log.write('mispriming checking: done\n')
##########################################################################################
