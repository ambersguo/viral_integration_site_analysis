#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re

###################internal_filter##############################################
def internal_filter (IN,OUT,LOG):
    file2in = open(IN, 'r')
    file2out = open(OUT, 'w')

    internal = noninternal = 0
    for line in file2in:
        eles = re.split(r'\t', line.rstrip('\n'))
        if re.search('GTGGCGCC|GTCCCCCC', eles[3]):
            internal += 1
        else:
            file2out.write(line)
        noninternal += 1

    file2in.close()
    file2out.close()

    with open(LOG, 'a') as file2log:
        file2log.write('reads internal: {}\n'.format(internal))
        file2log.write('reads noninternal: {}\n'.format(noninternal))
################################################################################



###################len20_filter#################################################
def len20_filter (IN,OUT,LOG):
    file2in = open(IN, 'r')
    file2out = open(OUT, 'w')

    short = long = 0
    for line in file2in:
        eles = re.split(r'\t', line.rstrip('\n'))
        seq = eles[3]
        pre_linker_len = 100
        if re.search('AGTCCTCTAAG', seq):
            pre_linker, post_linker = seq.split('AGTCCTCTAAG',1)
            pre_linker_len = len(pre_linker)
            if pre_linker_len < 20:
                short += 1
            else:
                file2out.write(line)
                long += 1
        else:
            file2out.write(line)
            long += 1

    file2in.close()
    file2out.close()

    with open(LOG, 'a') as file2log:
        file2log.write('reads noninternal length less than 20bp: {}\n'.format(short))
        file2log.write('reads noninternal length equal or more than 20bp: {}\n'.format(long))
################################################################################



###################q20_filter###################################################
def q20_filter (IN,OUT,LOG):
    file2in = open(IN, 'r')
    file2out = open(OUT, 'w')

    qscore_cutoff = 20
    total = total_pass = 0
    average_ltr_qscore20 = average_linker_qscore20 = 0

    for line in file2in:
        sum = sum_linker = 0
        total += 1
        eles = re.split(r'\t', line.rstrip('\n'))
        pre_ltr = eles[2]
        ltr_qscore = eles[4]
        pre_ltr_len = len(pre_ltr)
        ltr_qscore20 = ltr_qscore[pre_ltr_len: pre_ltr_len+20] #get quality score of the first 20 bp of the LTR side sequence
        ltr_qscore20_list = list(ltr_qscore20)
        for n_asc in ltr_qscore20_list: #convert each ascii to numerical value and calculate average
            n = ord(n_asc)
            sum = sum + n
        average_ltr_qscore20 = sum/20-33 #start from ASCII 33 in Illumina new version

        #get quality score of linker side sequence
        linker_qscore= eles[17]
        linker_qscore20 = linker_qscore[25:45]
        linker_qscore20_list = list(linker_qscore20)
        for l_asc in linker_qscore20_list:
            l = ord(l_asc)
            sum_linker = sum_linker + l
        average_linker_qscore20 = sum/20-33
        if average_ltr_qscore20 > qscore_cutoff and average_linker_qscore20 > qscore_cutoff:
            file2out.write(line)
            total_pass += 1

    file2in.close()
    file2out.close()

    with open(LOG, 'a') as file2log:
        file2log.write('pre q20 filter: {}\n'.format(total))
        file2log.write('post q20 filter: {}\n'.format(total_pass))
################################################################################



#######################merger###################################################
def merger (IN,OUT,LOG):
    file2in = open(IN, 'r')
    file2out = open(OUT, 'w')

    uniq = []
    seq2line = {}
    seq2count = {}
    newline = ''

    for line in file2in:
        eles = re.split(r'\t', line.rstrip('\n'))
        seqR = eles[15]
        if seqR not in uniq:
            uniq.append(seqR)
            seq2line[seqR] = line
            seq2count[seqR] = 0
        seq2count[seqR] += 1

    for seq in uniq:
        line = seq2line[seq]
        count = seq2count[seq]
        eles = re.split(r'\t', line)
        seqF_prejunc = eles[2]
        seqF_postjunc = eles[3]
        seqF_qual = eles[4]
        idR = eles[13]
        linker = eles[16]
        seqR_qual = eles[17]
        seqR_pre_linker = eles[15].replace(linker,'')
        newline = idR + '#' + str(count) + '\t' + seqF_prejunc + '\t' + seqF_postjunc + '\t' + seqF_qual + '\t' + \
                  seqR_pre_linker + '\t' + linker + '\t' + seqR_qual + '\n'
        file2out.write(newline)

    file2in.close()
    file2out.close()

    with open(LOG, 'a') as file2log:
        file2log.write('merge forward and reverse reads: done\n')
################################################################################



#######################fastq_generator##########################################
def fastq_generator (IN,OUTFfq,OUTRfq,OUTFfa,OUTRfa,LOG):
    file2in = open(IN, 'r')
    file2outFfq = open(OUTFfq, 'w')
    file2outRfq = open(OUTRfq, 'w')
    file2outFfa = open(OUTFfa, 'w')
    file2outRfa = open(OUTRfa, 'w')

    uniq = []
    seq2line = {}
    seq2count = {}
    newlineF = ''
    newlineR = ''

    for line in file2in:
        eles = re.split(r'\t', line.rstrip('\n'))
        seqR = eles[15]
        if seqR not in uniq:
            uniq.append(seqR)
            seq2line[seqR] = line
            seq2count[seqR] = 0
        seq2count[seqR] += 1

    for seq in uniq:
        line = seq2line[seq]
        count = seq2count[seq]
        eles = re.split(r'\t', line)
        idF_fa = eles[0]
        idF_fq = eles[5]
        seqF = eles[3]
        seqF_qual = eles[4]
        idR_fa = eles[13]
        idR_fq = eles[18]
        seqR = eles[16]
        seqR_qual = eles[17]
        newlineF_fq = idF_fq + '\n' + seqF + '\n' + '+\n' + seqF_qual + '\n'
        newlineR_fq = idR_fq + '\n' + seqR + '\n' + '+\n' + seqR_qual + '\n'
        newlineF_fa = '>' + idF_fa + '#' + str(count) + '\n' + seqF + '\n'
        newlineR_fa = '>' + idR_fa + '#' + str(count) + '\n' + seqR + '\n'
        file2outFfq.write(newlineF_fq)
        file2outRfq.write(newlineR_fq)
        file2outFfa.write(newlineF_fa)
        file2outRfa.write(newlineR_fa) 

    file2in.close()
    file2outFfq.close()
    file2outRfq.close()
    file2outFfa.close()
    file2outRfa.close()

    with open(LOG, 'a') as file2log:
        file2log.write('generate fasta files for alignment: done\n')

################################################################################
