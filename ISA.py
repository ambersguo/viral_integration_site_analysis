#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import operator
import subprocess
import re
from pre_alignment_filter_bwa import internal_filter, len20_filter, q20_filter, merger,fastq_generator
from post_alignment_filter_bwa import bwa_filter, dogtag_fetcher, IS_BP_fetcher, unique_BP_dogtag_fetcher, nonBP_merger
from UMI_interpreter import dogtag_interpreter, dogtag_correction, dogtag_merger1, dogtag_merger2, dogtag_merger3
from post_UMI_filter import seq_id_fetcher, ref_seq_fetcher, seq_fetcher, misprime_inspector

current_path = os.path.abspath(os.getcwd())
path_eles = re.split(r'/', current_path)
LTR = path_eles[-1]
sample = path_eles[-2]
setting_path_eles = path_eles[:-2]
setting_path = '/'.join(setting_path_eles) + '/settings/setting.txt'

with open(setting_path, 'r') as setting:
    for linein in setting:
        if sample in linein:
            eles = re.split(r'\t', linein)
            i7_barcode_rc = eles[5]
            if LTR == '3LTR':
                i5_barcode = eles[1]
                inline_barcode = eles[3]
                primer_seq = eles[-2][0:12]
            elif LTR == '5LTR':
                i5_barcode = eles[2]
                inline_barcode = eles[4]
                primer_seq = eles[-1][0:12]

debarcoded2open = '1' + inline_barcode + '2' + i7_barcode_rc + '.txt'

def main():
    #pre alignment processing
    #internal_filter(debarcoded2open,'reads_debarcoded_noninternal_linker.txt','log.txt')
    #len20_filter('reads_debarcoded_noninternal_linker.txt','reads_debarcoded_noninternal_linker_long.txt','log.txt')
    #q20_filter('reads_debarcoded_noninternal_linker_long.txt','reads_debarcoded_noninternal_linker_long_q20.txt','log.txt')
    #merger('reads_debarcoded_noninternal_linker_long_q20.txt','reads_debarcoded_noninternal_linker_long_q20_merged.txt','log.txt')
    #fastq_generator('reads_debarcoded_noninternal_linker_long_q20.txt','reads_debarcoded_noninternal_linker_long_q20_F.fastq','reads_debarcoded_noninternal_linker_long_q20_R.fastq','reads_debarcoded_noninternal_linker_long_q20_F.fa','reads_debarcoded_noninternal_linker_long_q20_R.fa','log.txt')
    
    #align to hg19
    #subprocess.call("bwa mem -M -t 16 /Users/guos2/Genomes/bwa_index/hg19/hg19.fa reads_debarcoded_noninternal_linker_long_q20_F.fastq > reads_debarcoded_noninternal_linker_long_q20_F_bwa.sam",shell=True)
    #subprocess.call("bwa mem -M -t 16 /Users/guos2/Genomes/bwa_index/hg19/hg19.fa reads_debarcoded_noninternal_linker_long_q20_R.fastq > reads_debarcoded_noninternal_linker_long_q20_R_bwa.sam",shell=True)
    
    #post alignment processing
    #bwa_filter('_bwa.sam','_bwa_filter.txt','log.txt')
    #dogtag_fetcher('reads_debarcoded_noninternal_linker_long_q20_merged.txt','_bwa_filter.txt','_bwa_filter_dogtag.txt','log.txt')
    #IS_BP_fetcher('reads_debarcoded_noninternal_linker_long_q20_F_bwa_filter_dogtag.txt','reads_debarcoded_noninternal_linker_long_q20_R_bwa_filter_dogtag.txt','reads_debarcoded_noninternal_linker_long_q20_F_R_bwa.txt','log.txt')
    #unique_BP_dogtag_fetcher('reads_debarcoded_noninternal_linker_long_q20_F_R_bwa.txt','reads_debarcoded_noninternal_linker_long_q20_F_R_bwa_summary.txt','log.txt')
    #nonBP_merger('reads_debarcoded_noninternal_linker_long_q20_F_R_bwa_summary.txt','reads_debarcoded_noninternal_linker_long_q20_F_R_bwa_non_BP_merged.txt','log.txt')
    
    #UMI/dogtag interpretation
    #dogtag_interpreter('reads_debarcoded_noninternal_linker_long_q20_F_R_bwa_non_BP_merged.txt','reads_debarcoded_noninternal_linker_long_q20_F_R_bwa_major_dogtags.txt','log.txt')
    #dogtag_correction('reads_debarcoded_noninternal_linker_long_q20_F_R_bwa_major_dogtags.txt','reads_debarcoded_noninternal_linker_long_q20_F_R_bwa_major_dogtags_corrected.txt','log.txt')
    #dogtag_merger1('reads_debarcoded_noninternal_linker_long_q20_F_R_bwa_major_dogtags_corrected.txt','reads_debarcoded_noninternal_linker_long_q20_F_R_bwa_major_dogtags_corrected_merge1.txt','log.txt')
    #dogtag_merger2('reads_debarcoded_noninternal_linker_long_q20_F_R_bwa_major_dogtags_corrected_merge1.txt','reads_debarcoded_noninternal_linker_long_q20_F_R_bwa_major_dogtags_corrected_merge2.txt','log.txt')
    #dogtag_merger3('reads_debarcoded_noninternal_linker_long_q20_F_R_bwa_major_dogtags_corrected_merge2.txt','reads_debarcoded_noninternal_linker_long_q20_F_R_bwa_major_dogtags_corrected_merge3.txt','log.txt')
    
    #post UMI/dogtag processing
    seq_id_fetcher('reads_debarcoded_noninternal_linker_long_q20_F_R_bwa_major_dogtags_corrected_merge3.txt','reads_debarcoded_noninternal_linker_long_q20_F_R_bwa.txt','reads_debarcoded_noninternal_linker_long_q20_F_R_bwa_major_dogtags_corrected.txt','reads_debarcoded_noninternal_linker_long_q20_F_R_bwa_major_dogtags_corrected_merged_formatted1.txt','log.txt')
    ref_seq_fetcher('reads_debarcoded_noninternal_linker_long_q20_F_R_bwa_major_dogtags_corrected_merged_formatted1.txt','reads_debarcoded_noninternal_linker_long_q20_F_R_bwa_major_dogtags_corrected_merged_formatted2.txt','log.txt')
    seq_fetcher('reads_debarcoded_noninternal_linker_long_q20_merged.txt','reads_debarcoded_noninternal_linker_long_q20_F_R_bwa_major_dogtags_corrected_merged_formatted2.txt','reads_debarcoded_noninternal_linker_long_q20_F_R_bwa_major_dogtags_corrected_merged_formatted3.txt','log.txt')
    misprime_inspector(primer_seq,'reads_debarcoded_noninternal_linker_long_q20_F_R_bwa_major_dogtags_corrected_merged_formatted3.txt','reads_debarcoded_noninternal_linker_long_q20_F_R_bwa_major_dogtags_corrected_merged_formatted4.txt','log.txt')



if __name__ == '__main__':
    sys.exit(main())
