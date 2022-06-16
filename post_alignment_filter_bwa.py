#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import pandas as pd

##################################bwa_filter##############################################
def bwa_filter (INext,OUText,LOG):
    
    path = './'
    filelist = os.listdir(path)

    for file in filelist:
        if file.endswith(INext):
            with open(path+file, 'r') as file2in:
                eles = []
                lines = []
                for line in file2in:
                    line = line.rstrip('\n')
                    eles = re.split(r'\t', line)
                    eles = [int(s) if s.isdigit() else s for s in eles]

                    if not eles[0].startswith('@'):
                        mapQ = eles[4]
                        cigar = eles[5]
                        mdtag = eles[12]
                        if mapQ > 0 and mapQ != 255 and cigar != '*' and mdtag.startswith('MD'):
                            val = ''
                            val_match = val_non_match = 0
                    
                            for c in cigar:
                                if c.isdigit():
                                    val += c
                                else:
                                    if c == 'M':
                                        val_match = int(val)
                                        break
                                    else:
                                        val_non_match += int(val)
                                        val = ''

                            mdtag_new = ''
                            identity = 0
                            for c in mdtag:
                                if not c.isdigit():
                                    c = ','
                                mdtag_new += c
                            eles_mdtag = re.split(r',', mdtag_new)
                            for ele in eles_mdtag:
                                if ele.isdigit():
                                    identity += int(ele)

                            if mapQ > 10 and val_non_match < 4 and identity/val_match > 0.95:
                                newline = line + '\t' + str(val_match)
                                lines.append(newline)

                df = pd.DataFrame([l.split('\t') for l in lines])
                newfile = file
                newfile = newfile.replace(INext,OUText)
                df.to_csv(newfile, header=None, index=None, sep='\t')

    with open(LOG, 'a') as file2log:
        file2log.write('post alignment filter(mapping start in first 3bp and mapping quality larger than 10 and percent identity larger than 95%: done\n')
##########################################################################################     



##################################dogtag_fetcher##########################################
def dogtag_fetcher (REF,INext,OUText,LOG):

    path2in = './'
    filelist = os.listdir(path2in)

    id2dogtag = {}
    id2readcount = {}

    with open(REF, 'r') as file2ref:
        for linein in file2ref:
            linein = linein.rstrip('\n')
            eles = re.split(r'\t', linein)
            id_dogtag_count = eles[0]
            (readid,dogtag_count) = id_dogtag_count.split('&')
            (dogtag,readcount) = dogtag_count.split('#')
            id2dogtag[readid] = dogtag
            id2readcount[readid] = readcount

    for file2in in filelist:
        if file2in.endswith(INext):
            newfile = file2in
            file2outname = newfile.replace(INext,OUText)
            newlines = []    
            with open(path2in+file2in, 'r') as file2in:
                for linein in file2in:
                    linein = linein.rstrip('\n')
                    linein = linein.rstrip('\t')
                    eles = re.split(r'\t', linein)
                    readid = eles[0]
                    dogtag = id2dogtag[readid]
                    readcount = id2readcount[readid]
                    newline = linein + '\t' + dogtag + '\t' + readcount + '\n'
                    newlines.append(newline)
        
            with open(file2outname,'w') as file2out:
                for newline in newlines:
                    file2out.write(newline)
    with open(LOG, 'a') as file2log:
        file2log.write('read count and dogtag fetching: done\n')    
##########################################################################################



##################################IS_BP_fetcher###########################################
def IS_BP_fetcher (INF,INR,OUT,LOG):

    id2corF = {}
    id2corR = {}
    id2dogtagF = {}
    id2dogtagR = {}
    bestmatches = []

    with open('reads_debarcoded_noninternal_linker_long_q20_F_bwa_filter_dogtag.txt', 'r') as file2inF:
        for linein in file2inF:
            linein = linein.rstrip('\n')
            eles = re.split(r'\t', linein)
            readid = eles[0] + '#' + eles[-1]
            val_match = int(eles[-3])
            dogtag = eles[-2]
            id2dogtagF[readid] = dogtag
            if eles[1] == '16':
                chrcor = 'chr' + eles[2] + '\t-\t' + str(int(eles[3])-1+val_match)
            else:
                chrcor = 'chr' + eles[2] + '\t+\t' + str(int(eles[3])-1)
            id2corF[readid] = chrcor

    with open('reads_debarcoded_noninternal_linker_long_q20_R_bwa_filter_dogtag.txt', 'r') as file2inR:
        for linein in file2inR:
            linein = linein.rstrip('\n')
            eles = re.split(r'\t', linein)
            readid = eles[0] + '#' + eles[-1]
            val_match = int(eles[-3])
            dogtag = eles[-2]
            id2dogtagR[readid] = dogtag
            if eles[1] == '16':
                chrcor = 'chr' + eles[2] + '\t-\t' + str(int(eles[3])-1+val_match)
            else:
                chrcor = 'chr' + eles[2] + '\t+\t' + str(int(eles[3])-1)
            id2corR[readid] = chrcor

    for key, value in id2corF.items():
        bestmatch = key + '\t' + value
        bestmatches.append(bestmatch)

    df = pd.DataFrame([l.split('\t') for l in bestmatches])
    df[3] = df[3].astype(int)
    df_sorted = df.sort_values(by=[1,3], ascending = [True, True])
    bestlists_sorted = df_sorted.values.tolist()

    with open('reads_debarcoded_noninternal_linker_long_q20_F_R_bwa.txt', 'w') as file2out:
        for bestlist in bestlists_sorted:
            bestmatch = '\t'.join([str(ele) for ele in bestlist])
            readid = bestlist[0]
            if readid in id2dogtagR:
                file2out.write(bestmatch + '\t' + id2corR[readid] + '\t' + id2dogtagR[readid] + '\n')
            else:
                file2out.write(bestmatch + '\t\t\t\t' + id2dogtagF[readid] + '\n')
    
    with open(LOG, 'a') as file2log:
        file2log.write('IS and breakpoints fetching: done\n')
##########################################################################################



##################################unique_BP_dogtag_fetcher################################
def unique_BP_dogtag_fetcher (IN,OUT,LOG):

    cordogtag2count = {}
    uniq = []

    with open(IN, 'r') as file2in:
        for line in file2in:
            line = line.rstrip('\n')
            line = line.replace('#','\t')
            eles = re.split(r'\t', line)
            if eles[7]:
                len = abs(int(eles[7]) - int(eles[4]))
                if len > 10000:
                    len_str = ''
                else:
                    len_str = str(len)
            else:
                len_str = ''
            cordogtag = eles[2] + '\t' + eles[3] + '\t' + eles[4] + '\t' + eles[5] + '\t' + eles[6] + '\t' + eles[7] + '\t' + len_str + '\t' + eles[8]
            if cordogtag not in uniq:
                uniq.append(cordogtag)
                cordogtag2count[cordogtag] = 0
            cordogtag2count[cordogtag] = cordogtag2count[cordogtag] + int(eles[1])

    with open(OUT, 'w') as file2out:        
        for cordogtag in uniq:
            file2out.write(cordogtag + '\t' + str(cordogtag2count[cordogtag]) + '\n')
            
    with open(LOG, 'a') as file2log:
        file2log.write('unique IS-BP-dogtag pair fetching: done\n')
##########################################################################################



##################################nonBP_merger############################################
def nonBP_merger(IN, OUT, LOG):

    corFdogtag2count = {}
    lines = []
    uniq = []

    with open(IN, 'r') as file2in:
        for line in file2in:
            line = line.rstrip('\n')
            eles = re.split(r'\t', line)
            corF = eles[0] + eles[1] + eles[2]
            corF2 = corF + '*'
            length = eles[6]
            dogtag = eles[7]
            count = int(eles[8])
            corFdogtag = corF2 + dogtag
            if length:
                lines.append(line)
            else:
                if corFdogtag not in uniq:
                    uniq.append(corFdogtag)
                    corFdogtag2count[corFdogtag] = 0
                corFdogtag2count[corFdogtag] = corFdogtag2count[corFdogtag] + count

    with open(OUT, 'w') as file2out: 
	    for line in lines:
		    eles = re.split(r'\t', line)
		    corF = eles[0] + eles[1] + eles[2]
		    corF2 = corF + '*'
		    length = eles[6]
		    dogtag = eles[7]
		    count = int(eles[8])
		    corFdogtag = corF2 + dogtag
		    if corFdogtag in corFdogtag2count:
			    count = count + corFdogtag2count[corFdogtag]
		    newline = eles[0] + '\t' + eles[1] + '\t' + eles[2] + '\t' + eles[3] + '\t' + eles[4] + '\t' + eles[5] + '\t' + eles[6] + '\t' + eles[7] + '\t' + str(count) + '\n'
		    file2out.write(newline)
		    
    with open(LOG, 'a') as file2log:
        file2log.write('non BP IS-dogtag pair merging: done\n')
##########################################################################################
