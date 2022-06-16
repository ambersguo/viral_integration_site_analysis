#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import pandas as pd
import Levenshtein
import collections
import operator

##################################dogtag_interpreter######################################
def dogtag_interpreter (IN,OUT,LOG):

    corFlist = []
    corFdogtaglist = []
    corFdict = collections.defaultdict(list)
    corF2dict = collections.defaultdict(list)
    corRdict = collections.defaultdict(list)
    corFdogtag2count = {}
    corRdogtag2count = {}
    lines = []
    corF2major = {}
    corR2major = {}
    newlines = []

    with open(IN, 'r') as file2in:
        for line in file2in:
            line = line.rstrip('\n')
            eles = re.split(r'\t', line)
            corF = eles[0] + eles[1] + eles[2]
            corF2 = corF + '*'
            corR = eles[0] + eles[1] + eles[2] + '&' + eles[6]
            dogtag = eles[7]
            corFdogtag = corF + '&' + dogtag
            corRdogtag = corR + '&' + dogtag
            count = int(eles[8])
        
            corFdict[corF].append(corR)
            if dogtag not in corF2dict[corF2]:
                corF2dict[corF2].append(dogtag)
            corRdict[corR].append(dogtag)
            corRdogtag2count[corRdogtag] = count

            if corF not in corFlist:
                corFlist.append(corF)
            if corFdogtag not in corFdogtaglist:
                corFdogtaglist.append(corFdogtag)
                corFdogtag2count[corFdogtag] = 0
            corFdogtag2count[corFdogtag] = corFdogtag2count[corFdogtag] + count
            lines.append(line)

    for corF in corFlist:
        single_corFdogtag2count = {}
        for corR in corFdict[corF]:
            single_corRdogtag2count = {}
            for dogtag in corRdict[corR]:
                corRdogtag = corR + '&' + dogtag
                count = corRdogtag2count[corRdogtag];
                single_corRdogtag2count[corRdogtag] = count;
            major_corRdogtag = max(single_corRdogtag2count.items(), key=operator.itemgetter(1))[0]
            corR2major[corR] = major_corRdogtag[-10:]
  
        corF2 = corF + '*'
        for dogtag in corF2dict[corF2]:
            corFdogtag = corF + '&' + dogtag
            single_corFdogtag2count[corFdogtag] = corFdogtag2count[corFdogtag] 
        major_corFdogtag = max(single_corFdogtag2count.items(), key=operator.itemgetter(1))[0]
        corF2major[corF] = major_corFdogtag[-10:]  

    for line in lines:
        eles = re.split(r'\t', line)
        corF = eles[0] + eles[1] + eles[2]
        corR = eles[0] + eles[1] + eles[2] + '&' + eles[6]
        major_dogtag = corF2major[corF]
        major_dogtag_corR = corR2major[corR]
        newline = line + '\t' + major_dogtag_corR + '\t' + major_dogtag
        newlines.append(newline)
   
    df = pd.DataFrame([l.split('\t') for l in newlines])
    df[8] = df[8].astype(int)
    df_sorted = df.sort_values(by=[0,1,2,6,8], ascending = [True, True, True, False, False])
    df_sorted = df_sorted.astype(str)
    newlineslist = df_sorted.values.tolist()
        
    with open(OUT, 'w') as file2out:
        file2out.write('chr\tstrand\tIS\tchr\tstrand\tBP\tlength\tdogtag\tcount\tBP_major_dogtag\tIS_major_dogtag\n')
        for newlinelist in newlineslist:
            newline = '\t'.join([str(ele) for ele in newlinelist])
            file2out.write(newline + '\n')

    with open(LOG, 'a') as file2log:
        file2log.write('dogtag interpretation: done\n')
##########################################################################################



##################################dogtag_correction#######################################
def dogtag_correction (IN,OUT,LOG):

    file2out = open(OUT, 'w')
    with open(IN, 'r') as file2in:
	    header = file2in.readline()
	    header = header.rstrip('\n')
	    file2out.write(header + '\tcorrected_dogtag\tmark\n')
	    for line in file2in:
		    line = line.rstrip('\n')
		    eles = re.split(r'\t', line)
		    BPdogtag = eles[9]
		    ISdogtag = eles[10]
		    if BPdogtag == ISdogtag:
			    dogtag_corrected = BPdogtag
			    mark = 'BP_major_dogtag_same_as_IS_major_dogtag'
		    else:
			    dogtag_diff = Levenshtein.distance(BPdogtag, ISdogtag)
			    if dogtag_diff < 3 or dogtag_diff == 3:
				    dogtag_corrected = ISdogtag
				    mark = 'BP_dogtag_corrected_due_to_similarity_with_IS_major_dogtag'
			    else:
				    dogtag_corrected = BPdogtag
				    mark = 'BP_dogtag_different_from_IS_major_dogtag'
		    file2out.write(line + '\t' + dogtag_corrected + '\t' + mark + '\n')
    file2out.close()

    with open(LOG, 'a') as file2log:
        file2log.write('dogtag correction: done\n')
##########################################################################################



##################################dogtag_merger1##########################################
def dogtag_merger1 (IN,OUT,LOG):

    corRdogtag2count = {}
    with open(IN, 'r') as file2in:
	    firstline = file2in.readline()
	    for line in file2in:
		    line = line.rstrip('\n')
		    eles = re.split(r'\t', line)
		    corR = eles[0] + eles[1] + eles[2] + '&' + eles[6]
		    count = int(eles[8])
		    dogtag_corrected = eles[11]
		    corRdogtag = corR + '&' + dogtag_corrected
		    if corRdogtag not in corRdogtag2count:
			    corRdogtag2count[corRdogtag] = 0
		    corRdogtag2count[corRdogtag] = corRdogtag2count[corRdogtag] + count

    file2out = open(OUT, 'w')
    with open(IN, 'r') as file2in:
	    header = file2in.readline()
	    header = header.rstrip('\n')
	    file2out.write(header + '\tcount_unique_dogtag\n')
	    for line in file2in:
		    line = line.rstrip('\n')
		    eles = re.split(r'\t', line)    
		    corR = eles[0] + eles[1] + eles[2] + '&' + eles[6]
		    count = int(eles[8])
		    BPdogtag = eles[9]
		    ISdogtag = eles[10]
		    dogtag_corrected = eles[11]
		    mark = eles[12]
		    corRdogtag = corR + '&' + dogtag_corrected
		    if 'different' in mark:
			    if corRdogtag2count[corRdogtag] > 10:
				    file2out.write(line + '\t' + str(corRdogtag2count[corRdogtag]) + '\n')
		    else:
			    file2out.write(line + '\t' + str(corRdogtag2count[corRdogtag]) + '\n')
    file2out.close()

    with open(LOG, 'a') as file2log:
        file2log.write('dogtag merger1: done\n')
##########################################################################################



##################################dogtag_merger2##########################################
def dogtag_merger2 (IN,OUT,LOG):

    corRdogtaglist = []
    corRdogtag2count = {}
    with open(IN, 'r') as file2in:
	    firstline = file2in.readline()
	    for line in file2in:
		    line = line.rstrip('\n')
		    eles = re.split(r'\t', line)
		    corR = eles[0] + '\t' + eles[1] + '\t' + eles[2] + '\t' + eles[3] + '\t' + eles[4] + '\t' + eles[5] + '\t' + eles[6]
		    dogtag_corrected = eles[11]
		    count = int(eles[13])
		    corRdogtag = corR + '\t' + dogtag_corrected
		    if corRdogtag not in corRdogtaglist:
		        corRdogtaglist.append(corRdogtag)
		        corRdogtag2count[corRdogtag] = count

    with open(OUT, 'w') as file2out:
        file2out.write('count\tchr\tstrand\tIS\tchr\tstrand\tBP\tlength\tdogtag\n')
        for corRdogtag in corRdogtaglist:
            file2out.write(str(corRdogtag2count[corRdogtag]) + '\t' + corRdogtag + '\n')

    with open(LOG, 'a') as file2log:
        file2log.write('dogtag merger2: done\n')
##########################################################################################



##################################dogtag_merger3##########################################
def dogtag_merger3 (IN,OUT,LOG):

    corFdogtaglist = []
    corRdogtag2count = {}
    corRdogtag2line = {}
    corFdict = collections.defaultdict(list)

    with open(IN, 'r') as file2in:
	    header = file2in.readline()
	    for line in file2in:
		    line = line.rstrip('\n')
		    eles = re.split(r'\t', line)
		    count = int(eles[0])
		    corF = eles[1] + eles[2] + eles[3]
		    corR = eles[1] + eles[2] + eles[3] + '&' + eles[7]
		    dogtag = eles[8]
		    corFdogtag = corF + '&' + dogtag
		    corRdogtag = corR + '&' + dogtag
		    if corFdogtag not in corFdogtaglist:
		        corFdogtaglist.append(corFdogtag)
		    corRdogtag2count[corRdogtag] = count
		    newline = eles[1] + '\t' + eles[2] + '\t' + eles[3] + '\t' + eles[4] + '\t' + eles[5] + '\t' + eles[6] + '\t' + eles[7] + '\t' + eles[8]
		    corRdogtag2line[corRdogtag] = newline
		    corFdict[corF].append(corR)

    with open(OUT, 'w') as file2out:
        file2out.write(header)
        for corFdogtag in corFdogtaglist:
            corF, dogtag = corFdogtag.split('&',2)
            single_corRdogtag2count = {}
            total = 0
            for corR in corFdict[corF]:
                corRdogtag = corR + '&' + dogtag
                if corRdogtag in corRdogtag2count.keys():
                    single_corRdogtag2count[corRdogtag] = corRdogtag2count[corRdogtag]  
                    total = total + corRdogtag2count[corRdogtag]
            major_corRdogtag = max(single_corRdogtag2count.items(), key=operator.itemgetter(1))[0]
            newline = corRdogtag2line[major_corRdogtag]
            file2out.write(str(total) + '\t' + newline + '\n')

    with open(LOG, 'a') as file2log:
        file2log.write('dogtag merger3: done\n')
##########################################################################################
