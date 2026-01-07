#!/usr/bin/env python3
##################################################
#
#    --- siRNA off-target Discovery Pipeline --- 
#  off-targeting potential and off-target prediction tool
#
#  Copyright 2016 Ferhat Alkan <ferro@rth.dk>
#
#  This file is part of siRNA off-target Discovery Pipeline.
#
#  siRNA off-target Discovery Pipeline is free software: you can 
#  redistribute it and/or modify it under the terms of the 
#  GNU General Public License as published by the 
#  Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  siRNA off-target Discovery Pipeline is distributed in the hope
#  that it will be useful, but WITHOUT ANY WARRANTY; without even 
#  the implied warranty of MERCHANTABILITY or 
#  FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with siRNA off-target Discovery Pipeline, see file COPYING.
#  If not, see <http://www.gnu.org/licenses/>.
#################################################

import sys
import os
import subprocess
import gzip
import re
from math import exp
import time
import argparse
import random
import string
import array

__author__ = "Ferhat Alkan: ferro@rth.dk"
# This is the general pipeline to analyze off-target candidates that are predicted by RIsearch2


TEMP_FILES = set()

def pipeline_exit(failed=False):
    for temp_file in TEMP_FILES:
        try:
            os.remove(temp_file)
        except:
            sys.stderr.write("#WARNING: "+temp_file+" could not be removed. It might not exist.\n")
    if failed:
        sys.stdout.write("#TERMINATED: CHECK STANDARD ERROR OUTPUT!\n")
        sys.stderr.write("#TERMINATED: CHECK STANDARD ERROR OUTPUT!\n")
        exit(1)
    else:
        sys.stderr.write("#DONE: Pipeline has completed running successfully\n")
        sys.stdout.write("#DONE: Pipeline has completed running successfully\n")
    exit()

# Create the argument parser for the pipeline
def get_parser():
    parser = argparse.ArgumentParser(description="siRNA off-target Discovery Pipeline ")
    parser.add_argument("-r",
                        metavar='<file>',
                        help="RIsearch2 interaction prediction results for given siRNA (obtained in -p2 or default format) [REQUIRED]",
                        type=str,
                        required=True)
    parser.add_argument("-type",
                        metavar='<gw/tw>',
                        help="Type of RIsearch2 predictions. gw for genome-wide and tw for transcriptome-wide (default: gw)",
                        type=str,
                        default="gw")
    parser.add_argument("-os",
                        metavar='<file>',
                        help="FASTA file that contains the sequence for siRNA/miRNA antisense(guide) strand [REQUIRED]",
                        type=str,
                        required=True)
    parser.add_argument("-q",
                        metavar='<str>',
                        help="ID of the query siRNA/miRNA [REQUIRED]",
                        type=str,
                        required=True)
    parser.add_argument("-o",
                        metavar='<file>',
                        help="write output to given file [REQUIRED]",
                        type=str, required=True)
    parser.add_argument("-rx",
                        metavar='<command>',
                        help="RIsearch2 executable (default: risearch2.x)",
                        type=str,
                        default='risearch2.x')
    parser.add_argument("-t",
                        metavar='<file>',
                        help="transcriptome file in .gtf, .gff or .bed format to assign expression values. See MANUAL for more info!",
                        type=str,
                        default='')
    parser.add_argument("-feature",
                        metavar='<str>',
                        help="specify the feature to select lines from .gtf, .gff files for transcriptome intersection (default: exon)",
                        type=str,
                        default='exon')
    parser.add_argument("-expmetric",
                        metavar='<str>',
                        help="select which field from .gtf, .gff files determines the expression abundance level (default: FPKM)",
                        type=str,
                        default='FPKM')
    parser.add_argument("-alpha",
                        metavar='<alpha>',
                        help="alpha parameter for partition function. Seperate multiple values by ';' (default:0.8)",
                        type=str,
                        default='0.8')
    parser.add_argument("-gamma",
                        metavar='<gamma>',
                        help="gamma parameter for partition function. Seperate multiple values by ';' (default:0.8)",
                        type=str,
                        default='0.8')
    parser.add_argument("-p",
                        metavar='<path>',
                        help="path to the folder where accessibility (opening energy) files reside (default: no opening energies)",
                        type=str,
                        default='')
    parser.add_argument("--offPs",
                        help="report off-target probabilities for individual transcripts.",
                        action="store_true")
    parser.add_argument("--sort",
                        help="Sort RIsearch2 siRNA-target interaction prediction file by strand, chromosome and location for increased performance while reading opening energies",
                        action="store_true")

    group_pos = parser.add_argument_group('# on-target Option Group 1 #\nDetect the on-target with given transcript ID or genomic location')
    group_pos.add_argument("-oi",
                        metavar='<str>',
                        help="on-target transcript ID",
                        type=str,
                        default='None')#,required=True)
    group_pos.add_argument("-on",
                        metavar='chr;sp;ep;str',
                        help="genomic location info for on-target region",
                        type=str,
                        default='')

    group_on_pre = parser.add_argument_group('# on-target Option Group 2 #\n# Predict on-targets and compute accessibility on given fasta file')
    group_on_pre.add_argument("-of",
                        metavar='<file>',
                        help="FASTA sequence file for on-target transcript",
                        type=str,
                        default='')
    group_on_pre.add_argument("-rp",
                        metavar='s,X;l,X;e,X',
                        help="RIsearch2 parameters for ontarget interaction predictions (default: s,1:12/6;e,-10;l,20;p2)",
                        type=str,
                        default='')
    group_on_pre.add_argument("-ap",
                        metavar='L,X;W,X',
                        help="parameters for computing RNAplfold opening energies (default: L,40;W,80)",
                        type=str,
                        default='')
    group_on_pre.add_argument("-oexp",
                        metavar='<float>',
                        help="abundance estimate (expression level) of ontarget transcript",
                        type=float,
                        default=0.0)

    group_on = parser.add_argument_group('# on-target Option Group 3 #\nUse precomputed ontarget predictions and accessibilities')
    group_on.add_argument("-op",
                        metavar='<file>',
                        help="RIsearch2 on-target prediction file for given siRNA",
                        type=str,
                        default='')
    group_on.add_argument("-oa",
                        metavar='<file>',
                        help="precomputed binary accessibility file for the opening energies of on-target region",
                        type=str,
                        default='')
    group_on.add_argument("-oexp2",
                        metavar='<float>',
                        help="abundance estimate (expression level) of ontarget transcript",
                        type=float,
                        default=0.0)

    group_experimental = parser.add_argument_group('Other useful/experimental options')
    group_experimental.add_argument("--less",
                        help="limit RIsearch2 predictions with sense strand, ignore interactions predicted on antisense strand",
                        action="store_true")
    group_experimental.add_argument("-chr",
                        metavar='<chromosome>',
                        help="limit RIsearch2 predictions by chromosome",
                        type=str,
                        default='')
    group_experimental.add_argument("-loc",
                        metavar='chrID;spos;epos',
                        help="limit RIsearch2 predictions by genomic location",
                        type=str,
                        default='')
    group_experimental.add_argument("-thr",
                        metavar='<float>',
                        help="introduce fixed hybridization energy threshold value (RIsearch2 interaction energy)",
                        type=float, default=None)
    group_experimental.add_argument("-tp",
                        metavar='<F>',
                        help="introduce dynamic hybridization energy threshold for interaction predictions (F * E_min)",
                        type=float, default=None)
    group_experimental.add_argument("--thrAfterOpEn", 
                        help="use free energies instead of hybridization energies when filtering. Free energy is the sum of opening energy and hybridization energy.",
                        action="store_true")
    group_experimental.add_argument("--intersection", 
                        help="do not delete the intersection file which contains the expressed RIsearch2 binding sites",
                        action="store_true")
    group_experimental.add_argument("-sa",
                        metavar='<file>',
                        help="saves opening energies of the expressed RIsearch2 binding sites in a given file",
                        type=str,
                        default='')
    group_experimental.add_argument("-theta",
                        metavar='<float>',
                        help="experimental theta parameter, separate multiple values by ';'",
                        type=str,
                        default='')
    return parser


# Random id generator
def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

R_SEMICOLON = re.compile(r'\s*;\s*')
R_COMMA     = re.compile(r'\s*,\s*')
R_KEYVALUE  = re.compile(r'(\s+|\s*=\s*)')

# Parse a single GTF/GFF line and return a dictionary #
def parse_gtf_gff_line(line):
    result = {}
    fields = line.rstrip().split('\t')
    for i, col in enumerate(['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame']):
        # Strip double and single quotes.
        value = fields[i].strip('"\'')
        # Return a list if the value has a comma.
        if ',' in value:
            value = re.split(R_COMMA, value)
        # These values are equivalent to None.
        elif value in ['', '.', 'NA']:
            value = None
        result[col] = value

    # INFO field consists of "key1=value;key2=value;...".
    infos = [x for x in re.split(R_SEMICOLON, fields[8]) if x.strip()]

    for i, info in enumerate(infos, 1):
        # It should be key="value".
        try:
            key, _, value = re.split(R_KEYVALUE, info, 1)
        # But sometimes it is just "value".
        except ValueError:
            key = 'INFO{}'.format(i)
            value = info

        # Strip double and single quotes.
        value = value.strip('"\'')
        # Return a list if the value has a comma.
        value = re.split(R_COMMA, value) if ',' in value else value
        if value not in ['', '.', 'NA',None]:
            result[key] = value

    return result

# Function that splits the line whether from RIsearch output or transcriptome intersection output
def split_prediction_line(given_line):
    given_cols = given_line.rstrip().split()
    if len(given_cols) == 9 or len(given_cols) == 8:  # RIsearch output file
        sirna, si_spos, si_epos, chromosome, spos, epos, strand, energy = given_cols[:8]
        return chromosome, int(spos), int(epos), strand, float(energy), 1.0, chromosome, chromosome
    elif len(given_cols) == 10:  # RIsearch output file with expression
        sirna, si_spos, si_epos, chromosome, spos, epos, strand, energy, structure, expression = given_cols
        return chromosome, int(spos), int(epos), strand, float(energy), float(expression), chromosome, chromosome
    elif len(given_cols) == 7:  # tw RIsearch output file intersected with expression
        chromosome, spos, epos, extra, energy, strand, expression = given_cols
        return chromosome, int(spos), int(epos), strand, float(energy), float(expression), chromosome, chromosome
    elif len(given_cols) > 12:  # Transcriptome output file
        ex_chr, ex_spos, ex_epos, ex_tid, expression, ex_str, ex_gid = given_cols[:7]
        ta_chr, ta_spos, ta_epos, si_l, score, strand, energy, length = given_cols[-8:]

        if ex_str == strand:
            return ta_chr, int(ta_spos), int(ta_epos), strand, float(energy), float(expression), ex_gid, ex_tid
        else:
            return None
    else:
        sys.stderr.write('\n#ERROR: Problem while reading intersection file. column size:' +
                         str(len(given_cols)) + '\n#TERMINATING\n')
        pipeline_exit(failed=True)

def fn_open(filepath):
    if filepath.endswith('.gz'):
        return gzip.open(filepath, 'rt')
    else:
        return open(filepath, 'rt')

# Transcriptome overlap function
def get_transcriptome_intersect(ris_file, trans_file, feature, score_type, Rtype, sort=False, keep=False):
    # Intersects transcriptome file with RIsearch2 predictions
    ris = ris_file.split('/')[-1]
    trans = trans_file.split('/')[-1]

    # Bed files to be used for intersectBed
    rand_id = id_generator()
    ris_bed = ris + '.' + rand_id + '.temp.bed'
    trans_bed = trans + '_for_' + ris + '.' + rand_id + '.temp.bed'
    intersection_file = ris + '_intersected_' + trans + '.targets'

    # Convert risearch output to bed format
    # RIsearch2: sirnaid sirnaspos sirnaepos chromosome(TranscriptID) spos epos strand enrgy str(only with p2)
    # BED: chromosome spos epos length energy strand
    TEMP_FILES.add(ris_bed)
    ris_in_f = fn_open(ris_file)
    with open(ris_bed, 'wt') as ris_out_f:
        for ris_line in ris_in_f:
            ris_cols = ris_line.split()
            ris_out_f.write('\t'.join([ris_cols[3], ris_cols[4], ris_cols[5],
                                       ris_cols[3]+'_'+str(int(ris_cols[2]) - int(ris_cols[1]) + 1),
                                       "1", ris_cols[6],ris_cols[7] + '\n']))
    ris_in_f.close()

    # Convert any transcriptome file to bed formatted file
    # chromosome spos epos tid score strand gid
    delete_trans_temp_bed_file = False
    if '.gtf' in trans_file or '.gff' in trans_file or '.GTF' in trans_file or '.GFF' in trans_file:
        in_f = fn_open(trans_file)
        with open(trans_bed, 'wt') as out_f:
            for line in in_f:
                t_dic = parse_gtf_gff_line(line)
                # We can add a line that sanity checks the dictionary keys for feature, seqname,start,end and score_type(FPKM)
                if "feature" in t_dic.keys() and t_dic["feature"] == feature:
                    if "seqname" in t_dic.keys() and "start" in t_dic.keys() and "end" in t_dic.keys() and "strand" in t_dic.keys() and score_type in t_dic.keys():
                        g_id = t_dic["gene_id"] if "gene_id" in t_dic.keys() else t_dic["seqname"]
                        t_id = t_dic["transcript_id"] if "transcript_id" in t_dic.keys() else g_id
                        if float(t_dic[score_type])>0.0:
                            out_f.write('\t'.join([t_dic["seqname"],t_dic["start"],t_dic["end"],t_id,"1",t_dic["strand"],t_dic[score_type],(g_id+'\n')]))
                    else:
                        sys.stderr.write('WARNING: Skipping line in .gtf/.gff transcriptome file, one or more fields missing (required fields: seqname, start, end, strand, '+score_type+')\n')
        delete_trans_temp_bed_file = True
        in_f.close()

    elif trans_file.endswith('.bed') or trans_file.endswith('.bed.gz') or trans_file.endswith('.BED') or trans_file.endswith('.BED.gz'):
        in_f = fn_open(trans_file)
        t_cols = in_f.readline().rstrip().split()
        if len(t_cols) < 6:
            sys.stderr.write('#ERROR: Transcriptome .bed file has less than 6 columns\n#TERMINATING\n')
            pipeline_exit(failed=True)
        in_f.close()
        trans_bed = trans_file
        delete_trans_temp_bed_file = False

    else:
        sys.stderr.write('ERROR: Unknown type for Transcriptome file (only .gtf,.gff or .bed)\n#TERMINATING\n')
        pipeline_exit(failed=True)

    if delete_trans_temp_bed_file:
        TEMP_FILES.add(trans_bed)

    if Rtype=="gw":
        # Intersect exon coordinates with risearch intervals
        TEMP_FILES.add(intersection_file + '.temp')
        with open(intersection_file + '.temp', "wt") as output:
            subprocess.call("bedtools intersect -wo -a %s -b %s" % (trans_bed, ris_bed), stdout=output, shell=True)

        if not keep:
            TEMP_FILES.add(intersection_file)
        with open(intersection_file + '.temp', 'rt') as int_in_f:
            with open(intersection_file,'wt') as int_out_f:
                prev = ('id', 'chr', 'spos', 'epos', 'str')
                for line in int_in_f:
                    cols = line.rstrip().split('\t')
                    cur = (cols[3], cols[-7], cols[-6], cols[-5], cols[-2])
                    if not (prev==cur):# or (prev[0] == cur[0]  and 'hap' in prev[1] and 'hap' in cur[1])):
                        int_out_f.write(line)
                    prev = cur

        if sort:
            with open(intersection_file + '.sorted', "wt") as output:
                subprocess.call("sort -k9,9 -k14,14 -k10n %s" % (intersection_file), stdout=output, shell=True)
            intersection_file = intersection_file + '.sorted'
            if not keep:
                TEMP_FILES.add(intersection_file)

    elif Rtype=="tw":
        expression = {}
        trans_in_f = fn_open(trans_bed) 
        for tr_line in trans_in_f:
            tr_cols=tr_line.split()
            expression[tr_cols[3]]=tr_cols[6]
        trans_in_f.close()

        if not keep:
            TEMP_FILES.add(intersection_file)

        ris_in_f = fn_open(ris_bed)
        with open(intersection_file, 'wt') as out_f:
            for ris_line in ris_in_f:
                ris_cols = ris_line.split()
                if ris_cols[0] in expression.keys():
                    out_f.write('\t'.join([ris_cols[0],ris_cols[1],ris_cols[2],ris_cols[3],ris_cols[6],ris_cols[5],expression[ris_cols[0]]+'\n']))
        ris_in_f.close()

        if sort:
            with open(intersection_file + '.sorted', "wt") as output:
                subprocess.call("sort -k1,1 -k6,6 -k2n %s" % (intersection_file), stdout=output, shell=True)
            intersection_file = intersection_file + '.sorted'
            if not keep:
                TEMP_FILES.add(intersection_file)

    else:
        sys.stderr.write("#ERROR: Wrong type selection for RIsearch2 results. Select gw or tw\n#TERMINATING\n")
        pipeline_exit(failed=True)

    return intersection_file

# Create the binary accessibility opening energy file for given sequence
def create_acc_bin(seq_file, seq_name, outfile, w=80, l=40):
    sf = fn_open(seq_file)
    seq = '-'
    for line in sf:
        if line[0]=='>':
            if line.rstrip()[1:]==seq_name:
                seq=''
            elif seq!='-':
                break
        elif seq!='-':
            seq = seq + line.rstrip()
    if len(seq)<w:
        w=len(seq)
    if len(seq)<l:
        l=len(seq)

    score_file = seq_name + '_openen'

    try:
        acc_cmd = "RNAplfold -W " + str(w) + " -L " + str(l) + " -u 30 -O < " + seq_file
        sub_output = subprocess.Popen(acc_cmd, shell=True, stdout=subprocess.PIPE).communicate()
    except:
        sys.stderr.write("#ERROR: RNAplfold could not be executed for on-target.\n#TERMINATING\n")
        sys.stderr.write("Unexpected error: "+str(sys.exc_info())+'\n')
        pipeline_exit(failed=True)
        
    with open(outfile, 'wt+') as out_f:
        with open(score_file, 'rt') as in_f:
            for sc_line in in_f:
                if sc_line[0] != '#' and sc_line[0] != ' ' and sc_line[0] != '':
                    sc_cols = sc_line.rstrip().split()[1:]
                    assert len(sc_cols) == 30
                    for num_str in sc_cols:
                        number = 0
                        if num_str != 'NA':
                            number = int(round(float(num_str) * 10.0))
                        if number > 255:
                            number = 255
                        elif number < 0:
                            number = 0
                        out_f.write(chr(number))
    TEMP_FILES.add(score_file)
    TEMP_FILES.add(seq_name + '_dp.ps')

# Function that returns the opening energy of the (spos-epos) region. Positions are inclusive (1-based)
def get_opening_energy(tid, spos, epos, strand, given_bin='', bin_folder='', acc_f=None):
    offset = epos - spos  # offset is length-1 where length is the binding site length
    value = -1.0
    #t1 = time.time()
    if strand == '+':
        filepath = bin_folder + '/' + tid + '.open.acc.bin'
        i_val_y = ((epos - 1) * 30) + offset # *30 because, RNAplfold parameter -u 30
        if given_bin != '':  # CHECK if bin is given
            filepath = given_bin

        if acc_f is None:
            try:
                acc_f = open(filepath, 'rb')
            except:
                sys.stderr.write("#ERROR: Opening energy file "+filepath+" can not be opened.\n#TERMINATING\n")
                sys.stderr.write("Unexpected error: "+str(sys.exc_info())+'\n')
                pipeline_exit(failed=True)
        elif filepath != acc_f.name:
            acc_f.close()
            try:
                acc_f = open(filepath, 'rb')
            except:
                sys.stderr.write("#ERROR: Opening energy file "+filepath+" can not be opened.\n#TERMINATING\n")
                sys.stderr.write("Unexpected error: "+str(sys.exc_info())+'\n')
                pipeline_exit(failed=True)

        acc_f.seek(i_val_y - acc_f.tell(),1)
        value = float(ord(acc_f.read(1))) / 10

    elif strand == '-':
        filepath = bin_folder + '/' + tid + '.rev.open.acc.bin'
        i_val_y = (30 * (spos-1)) + offset # *30 because, RNAplfold parameter -u 30
        if given_bin != '':  # CHECK if bin is given
            filepath = given_bin

        if acc_f is None:
            try:
                acc_f = open(filepath, 'rb')
            except:
                sys.stderr.write("#ERROR: Opening energy file "+filepath+" can not be opened.\n#TERMINATING\n")
                sys.stderr.write("Unexpected error: "+str(sys.exc_info())+'\n')
                pipeline_exit(failed=True)
        elif filepath != acc_f.name:
            acc_f.close()
            try:
                acc_f = open(filepath, 'rb')
            except:
                sys.stderr.write("#ERROR: Opening energy file "+filepath+" can not be opened.\n#TERMINATING\n")
                sys.stderr.write("Unexpected error: "+str(sys.exc_info())+'\n')
                pipeline_exit(failed=True)
        acc_f.seek(i_val_y - acc_f.tell(),1)
        value = float(ord(acc_f.read(1))) / 10

    else:
        sys.stderr.write('WARNING: False Strand Info:' + strand + '\n')

    return value, acc_f

# Returns the interactions from a given GW or TW input
def get_predictions(filepath, on_id='', on_chrom='', on_spos=-1, on_epos=-1, on_strand='',
                    threshold=None, thr_after_acc=False, bins='', less=False, chromosome_choice='', qid='', ris2='', seq_file=''):
    include_opening_energies = False
    if bins != '':
        include_opening_energies = True

    predictions = []
    on_target_predictions = []
    n = 0
    exc_open = 0
    exc_thr = 0
    exc_less = 0
    minimum = 0.0
    minimum_noaccs = 0.0
    
    # Minimum energy update
    if qid != '' and ris2 != '' and seq_file != '':
        seq = ''
        read = False
        with open(seq_file, 'rt') as f:
            for line in f:
                if line[0]=='>':
                    if read:
                        break
                    read = True if line.rstrip()[1:]==qid else False
                elif read:
                    seq += line.rstrip()
        if seq=="":
            sys.stderr.write('#ERROR: siRNA fasta file do not contain sequence with ID: '+args.q+'\n#TERMINATING\n')
            pipeline_exit(failed=True)
        seq = seq.replace('u','t').replace('U','T')
        temp = id_generator()

        TEMP_FILES.add(temp+'.fa')
        with open(temp+'.fa','w') as f:
            f.write('>'+temp+'\n'+seq+'\n')

        try:
            # GEt E_min
            temp_suf = temp+'.pksuf'
            ris_cmd = ris2 +" -c "+ temp +".fa -o "+ temp_suf
            TEMP_FILES.add(temp_suf)
            
            out_ris2_suf = subprocess.Popen(ris_cmd, shell=True, stdout=subprocess.PIPE).communicate()[0]

            ris_cmd = ris2 +" -q "+ temp +".fa -i "+ temp_suf +" -s "+str(len(seq)-1)+" -e 0"
            TEMP_FILES.add('risearch_'+temp+'.out.gz')
            out_ris2_pre = subprocess.Popen(ris_cmd, shell=True, stdout=subprocess.PIPE).communicate()[0]

            f = gzip.open('risearch_'+temp+'.out.gz', 'rt')
            minimum_noaccs = float(f.readline().rstrip().split('\t')[-1])
            for line in f:
                temp_eng = float(line.rstrip().split('\t')[-1])
                minimum_noaccs = temp_eng if temp_eng < minimum_noaccs else minimum_noaccs
            f.close()
            sys.stderr.write('#GETPREDICTIONS: siRNA sequence:'+seq+' minimum:'+str(minimum_noaccs)+'\n')
        except:
            sys.stderr.write('#ERROR: Could not run RIsearch2 to obtain E_min. Link to RIsearch2 executable might be broken.\n')
            sys.stderr.write("Unexpected error: "+str(sys.exc_info())+'\n')
            pipeline_exit(failed=True)
    # Minimum update finish
    
    temp_acc_f = None

    p_f = None
    try:
        p_f =  gzip.open(filepath, 'rt') if filepath[-3:] == '.gz' else open(filepath, 'rt')
    except:
        sys.stderr.write('#ERROR: Could not open the file '+filepath+"\n#TERMINATING\n")
        sys.stderr.write("Unexpected error: "+str(sys.exc_info())+'\n')
        pipeline_exit(failed=True)

    for p_line in p_f:
        # if n % 1000 == 0: # For tracking
        #    print n
        if split_prediction_line(p_line)==None:
            continue
        chromosome, spos, epos, strand, energy, expression, gid, tid = split_prediction_line(p_line)
        if (threshold is None) or (energy < threshold):
            if energy < minimum_noaccs:
                minimum_noaccs = energy
            if include_opening_energies:
                if ((less is False) or (strand == '+')) and ((chromosome_choice=='') or (chromosome==chromosome_choice)):
                    open_eng, temp_acc_f = get_opening_energy(chromosome, spos, epos, strand, bin_folder=bins, acc_f=temp_acc_f)
                    if open_eng >= 0.0:
                        if (thr_after_acc is False) or (threshold is None) or ((energy + open_eng) < threshold):
                            predictions.append((chromosome, spos, epos, strand, energy, open_eng, expression, gid, tid))
                            if (chromosome == on_chrom and on_strand == strand and on_spos < epos and on_epos > on_spos) or gid == on_id or tid == on_id:
                                on_target_predictions.append(n)
                            n += 1
                            if energy + open_eng < minimum:
                                minimum = energy + open_eng
                        else:
                            exc_thr += 1
                    else:
                        exc_open += 1
                        sys.stderr.write('Can not get opening energy for ' + ' '.join(
                            [chromosome, str(spos), str(epos), strand, gid, tid, '\n']))
                else:
                    exc_less += 1
            else:
                minimum = minimum_noaccs
                if ((less is False) or (strand == '+')) and ((chromosome_choice == '') or (chromosome==chromosome_choice)):
                    predictions.append((chromosome, spos, epos, strand, energy, 0.0, expression, gid, tid))
                    if (chromosome == on_chrom and on_strand == strand and on_spos < epos and on_epos > spos) or gid == on_id or tid == on_id:
                        on_target_predictions.append(n)
                    n += 1
                else:
                    exc_less += 1
        else:
            exc_thr += 1

    if temp_acc_f is not None:
        temp_acc_f.close()
    p_f.close()
    sys.stderr.write('#GETPREDICTIONS: '+str(exc_less) + ' predictions excluded (sense/antisense strand/chromosome choice); ' + str(exc_open) +
                     ' excluded (opening energy); ' + str(exc_thr) + ' excluded (energy threshold)\n')
    return predictions, on_target_predictions, minimum, minimum_noaccs


# Update the predictions list for on-targets
def add_update_on_targets(on_target_pred_file, on_target_acc_bin, predictions, on_target_predictions,
                          minimum, minimum_noaccs, oexp, threshold=None, thr_after_acc=False):
    include_opening_energies = False
    if on_target_acc_bin != '':
        include_opening_energies = True
    temp_acc_f = None
    for i in on_target_predictions:
        predictions[i] = ('', -1, -1, '', 0.0, 0.0, 0.0)
    # on_target_predictions = []  # optional
    p_f = gzip.open(on_target_pred_file, 'rt')
    for p_line in p_f:
        if split_prediction_line(p_line)==None:
            continue
        chromosome, spos, epos, strand, energy, expression, gid, tid = split_prediction_line(p_line)
        if (strand == '+') and ((threshold is None) or (energy < threshold)):
            if include_opening_energies:
                open_eng, temp_acc_f = get_opening_energy(chromosome, spos, epos, strand, given_bin=on_target_acc_bin, acc_f=temp_acc_f)
                if open_eng >= 0.0:
                    if (thr_after_acc is False) or (threshold is None) or ((energy + open_eng) < threshold):
                        on_target_predictions.append(len(predictions))
                        predictions.append((chromosome, spos, epos, strand, energy, open_eng, oexp, gid, tid))
                        if energy + open_eng < minimum:
                            minimum = energy + open_eng
                        if energy < minimum_noaccs:
                            minimum_noaccs = energy
                else:
                    sys.stderr.write('Can not get opening energy for on-target \n')
            else:
                if energy < minimum:
                    minimum = energy
                    minimum_noaccs = energy
                ontargets.append(len(predictions))
                predictions.append((chromosome, spos, epos, strand, energy, 0.0, oexp, gid, tid))
    p_f.close()

    if temp_acc_f is not None:
        temp_acc_f.close()
    # print on-targets
    return minimum, minimum_noaccs


def apply_tp_threshold(predictions, on_target_predictions, threshold_per, min, min_na, thr_after_openen=False):
    new_predictions = []
    new_on_targets = []
    for i in range(len(predictions)):
        if thr_after_openen is False:
            if predictions[i][4] <= min_na*threshold_per:
                if i in on_target_predictions:
                    new_on_targets.append(len(new_on_targets))
                new_predictions.append(predictions[i])
        else:
            if predictions[i][4]+predictions[i][5] <= min*threshold_per:
                if i in on_targets:
                    new_on_targets.append(len(new_on_targets))
                new_predictions.append(predictions[i])
    return new_predictions, new_on_targets


# Calculates the partition function (Value Z) for given predictions
def get_z(predictions, minimum, minimum_noaccs, parameter_list):
    kb = 0.001987
    t = 310.15
    par_beta = 1.0 / (kb * t)

    value_z = {}
    value_z_noaccs = {}
    for par in parameter_list:
        if type(par) is float or (type(par) is tuple and len(par) == 2):
            value_z[par] = 0.0
            value_z_noaccs[par] = 0.0
    for (chromosome, spos, epos, strand, energy, open_eng, expression, gid, tid) in predictions:
        for par in parameter_list:
            if type(par) is float:
                assigned_energy = ( ( par * (energy+10) ) - 10) + open_eng
                assigned_energy_noaccs = ( par * (energy+10) ) - 10
            elif type(par) is tuple and len(par) == 2:
                (alpha,gamma) = par
                assigned_energy = energy + open_eng
                assigned_energy_noaccs = energy
                if energy < alpha * minimum_noaccs:
                    assigned_energy = gamma * minimum_noaccs + open_eng
                    assigned_energy_noaccs = gamma * minimum_noaccs
            else:
                continue
            value_z[par] += expression * exp(-1.0 * par_beta * assigned_energy)
            value_z_noaccs[par] += expression * exp(-1.0 * par_beta * assigned_energy_noaccs)
    return value_z, value_z_noaccs


# Get probabilities for all predictions
def get_probabilities(predictions, on_targets, minimum, minimum_noaccs, parameter_list, off_ps=False):
    prob_scores = {}
    prob_scores_noaccs = {}
    kb = 0.001987
    t = 310.15
    par_beta = 1.0 / (kb * t)

    value_z, value_z_noaccs = get_z(predictions, minimum, minimum_noaccs, parameter_list)
    
    for par in parameter_list:
        if type(par) is float or (type(par) is tuple and len(par) == 2):
            prob_scores[par] = []
            prob_scores_noaccs[par] = []

    if off_ps:
        for (chromosome, spos, epos, strand, energy, open_eng, expression, gid, tid) in predictions:
            for par in parameter_list:
                if type(par) is float:
                    assigned_energy = ( ( par * (energy+10) ) - 10) + open_eng
                    assigned_energy_noaccs = ( par * (energy+10) ) - 10
                elif type(par) is tuple and len(par) == 2:
                    (alpha, gamma) = par
                    assigned_energy = energy + open_eng
                    assigned_energy_noaccs = energy
                    if energy < alpha * minimum_noaccs:
                        assigned_energy = gamma * minimum_noaccs + open_eng
                        assigned_energy_noaccs = gamma * minimum_noaccs
                else:
                    continue
                prob_scores[par].append((expression * exp(-1.0 * par_beta * assigned_energy)) / value_z[par])
                prob_scores_noaccs[par].append((expression * exp(-1.0 * par_beta * assigned_energy_noaccs)) / value_z_noaccs[par])

    on_p = {}
    on_p_noaccs = {}
    for par in parameter_list:
        on_target_prob = 0.0
        on_target_prob_noaccs = 0.0
        if type(par) is float:
            for i in on_targets:
                (chromosome, spos, epos, strand, energy, open_eng, expression, gid, tid) = predictions[i]
                assigned_energy = ( ( par * (energy+10) ) - 10) + open_eng
                assigned_energy_noaccs = ( ( par * (energy+10) ) - 10)
                on_target_prob += (expression * exp(-1.0 * par_beta * assigned_energy)) / value_z[par] 
                on_target_prob_noaccs += (expression * exp(-1.0 * par_beta * assigned_energy_noaccs)) / value_z_noaccs[par]
        elif type(par) is tuple and len(par) == 2:
            (alpha, gamma) = par
            for i in on_targets:
                (chromosome, spos, epos, strand, energy, open_eng, expression, gid, tid) = predictions[i]
                assigned_energy = energy + open_eng
                assigned_energy_noaccs = energy
                if energy < alpha * minimum_noaccs:
                    assigned_energy = gamma * minimum_noaccs + open_eng
                    assigned_energy_noaccs = gamma * minimum_noaccs
                on_target_prob += (expression * exp(-1.0 * par_beta * assigned_energy)) / value_z[par] 
                on_target_prob_noaccs += (expression * exp(-1.0 * par_beta * assigned_energy_noaccs)) / value_z_noaccs[par]
        else:
            continue
        on_p[par] =  on_target_prob
        on_p_noaccs[par] =  on_target_prob_noaccs
    return prob_scores, prob_scores_noaccs, value_z, value_z_noaccs, on_p, on_p_noaccs


# Prepare off-target dictionary from bed-file with min 6 columns
def get_off_target_dic(off_target_bed_path):
    off_dic = {}
    with open(off_target_bed_path, 'rt') as off_f:
        for off_line in off_f:
            off_cols = off_line.rstrip().split()
            if off_cols[3] in off_dic.keys():
                off_dic[off_cols[3]].append(off_cols[0], off_cols[1], off_cols[2], off_cols[5])
            else:
                off_dic[off_cols[3]] = [(off_cols[0], off_cols[1], off_cols[2], off_cols[5])]
    return off_dic


# Function to report results whether to file or screen
def report(predictions, prob_scores, prob_scores_noaccs, z_val, z_val_na, on_probs, on_probs_noaccs, parameter_list, query, to_file='', verbose=False):
    out = sys.stdout
    if to_file != '':
        out = open(to_file, 'w')

    # Report the on-target first
    out.write('# On-target info for '+query+' #\n')
    for par in parameter_list:
        if type(par) is float:
            out.write('# For theta='+str(par)+';')
        elif type(par) is tuple and len(par) == 2:
            (alpha, gamma) = par
            out.write('# For alpha='+str(alpha)+' and gamma='+str(gamma)+';')
        else:
            continue
        on_target_prob = on_probs[par]
        on_target_prob_noaccs = on_probs_noaccs[par]
        out.write(' Pon: ' + str(on_target_prob) + '; Pon_noacc: ' + str(on_target_prob_noaccs) +
                  '; Poff: ' + str(1.0-on_target_prob) + '; Poff_noacc: ' + str(1.0-on_target_prob_noaccs) +
                  '; Z: ' + str(z_val[par]) + '; Z_noacc: ' + str(z_val_na[par]) + 
                  '; Zoff: ' + str(z_val[par]*(1.0-on_target_prob)) + '; Zoff_noacc: ' + str(z_val_na[par]*(1.0-on_target_prob_noaccs))+ '\n')
    out.write('## End of on-target info ##\n')

    # Report off-target probabilities for default ids
    if verbose:
        out.write('# Off-target info#\n')
        off_probs = {}
        off_probs_noaccs = {}
        
        for i in range(len(predictions)):
            chrom, spos, epos, strand, energy, open_eng, expression, gid, tid = predictions[i]
            id_key = '\t'.join([gid, tid])
            if id_key not in off_probs.keys():
                off_probs[id_key] = {}
                off_probs_noaccs[id_key] = {}
                for par in parameter_list:
                    if type(par) is float or (type(par) is tuple and len(par) == 2):
                        off_probs[id_key][par] = 0.0
                        off_probs_noaccs[id_key][par] = 0.0

            for par in parameter_list:
                if type(par) is float or (type(par) is tuple and len(par) == 2):
                    off_probs[id_key][par] += prob_scores[par][i]
                    off_probs_noaccs[id_key][par] += prob_scores_noaccs[par][i]

        # Report off-target probabilities for default ids
        out.write('# ID1\tID2')
        for par in parameter_list:
            if type(par) is float:
                out.write('\ttheta:'+str(par)+'\ttheta:'+str(par)+'(NoAcc)')
            elif type(par) is tuple and len(par) == 2:
                out.write('\talpha:'+str(par[0])+',gamma:'+str(par[1])+'\t'+'alpha:'+str(par[0])+',gamma:'+str(par[1])+'(NoAcc)')
            else:
                continue
        out.write('\n')
        for id_key in off_probs.keys():
            out.write(id_key)
            for par in parameter_list:
                if type(par) is float or (type(par) is tuple and len(par) == 2):
                    out.write('\t'+str(off_probs[id_key][par])+'\t'+str(off_probs_noaccs[id_key][par]))
            out.write('\n')
        out.write('## End of off-target info\n\n')

    # Report off-target probabilities for given ids
    # LATER

    # Close the report file
    if to_file != '':
        out.close()
    return True


def run_pipeline(args):
    ont_id = args.oi
    on_info = args.on.split(';')
    ont_c, ont_s, ont_e, ont_str = ['', -1, -1, '']
    if len(on_info) == 4:
        ont_c = on_info[0]
        ont_s = int(on_info[1])
        ont_e = int(on_info[2])
        ont_str = on_info[3]
    elif len(on_info) > 1:
        sys.stderr.write('#ERROR: Wrong format for on-target location info\n')
        pipeline_exit(failed=True)

    off_predictions_file = args.r

    t1=time.time()

    # Step 2: Create the missing files and make necessary filepath changes
    # 2.1 Getting the on target files ready
    on_predictions_file = args.op
    on_accessibility_file = args.oa

    # Create if they are not precomputed
    if args.op == '':
        sys.stderr.write('#CHECKPOINT: On target prediction-file is not given\n')
        if args.of == '':
            sys.stderr.write('#CHECKPOINT: On-target fasta file is not given.' +
                             ' On-target info will be recognized with given Transcript ID or genomic location\n')
        elif args.os == '':
            sys.stderr.write('#CHECKPOINT: siRNA fasta file is not given.' +
                             ' On-target info will be recognized with given Transcript ID or genomic location\n')
        elif args.rx == '':
            sys.stderr.write('#ERROR: RIsearch2 executable command is missing\n#TERMINATING\n')
            pipeline_exit(failed=True)
        else:
            sys.stderr.write('#CHECKPOINT: On target prediction-file is computed with given sequence files\n')
            on_target_name = ''
            with open(args.of, 'rt') as of_f:
                for of_line in of_f:
                    if of_line[0] == '>':
                        on_target_name = of_line.split()[0].rstrip()[1:]
                        break
            if on_target_name == '':
                sys.stderr.write('#ERROR: wrong on-target sequence file (no sequence name)\n#TERMINATING\n')
                pipeline_exit(failed=True)

            ris2 = args.rx
            ris_cmd = ris2 + " -c " + args.of + " -o " + on_target_name + ".pksuf"
            try:
                # Create the .suf
                out_ris2_suf = subprocess.Popen(ris_cmd, shell=True, stdout=subprocess.PIPE).communicate()[0]
                # Create the remp fasta that includes only our siRNA sequence
                read=False
                seq=""
                with open(args.os, 'rt') as f:
                    for line in f:
                        if line[0]=='>':
                            if read:
                                break
                            read = True if line.rstrip()[1:]==args.q else False
                        elif read:
                            seq += line.rstrip()
                if seq=="":
                    sys.stderr.write('#ERROR: siRNA fasta file do not contain sequence with ID: '+args.q+'\n#TERMINATING\n')
                    pipeline_exit(failed=True)
                seq = seq.replace('u','t').replace('U','T')
                temp = args.q+'_only'

                TEMP_FILES.add(temp+'.fa')
                with open((temp+'.fa'),'wt') as f:
                    f.write('>'+temp+'\n'+seq+'\n')

                ris_cmd = ris2 + " -q " + temp + ".fa -i " + on_target_name + '.pksuf -s "1:12/6" -e -10 -l 20'
                if args.rp != '':
                    ris_cmd = ris2 + " -i " + on_target_name + ".pksuf"
                    for param in args.rp.split(';'):
                        par_cols = param.split(',')
                        if len(par_cols) > 2 or len(par_cols) == 0:
                            sys.stderr.write('#ERROR: wrong RIsearch2 parameters\n#TERMINATING\n')
                            TEMP_FILES.add(on_target_name + ".pksuf")
                            pipeline_exit(failed=True)
                        par_cols[0] = '-' + par_cols[0]
                        for s in par_cols:
                            ris_cmd = ris_cmd + " " + s

                out_ris2_pre = subprocess.Popen(ris_cmd, shell=True, stdout=subprocess.PIPE).communicate()[0]
                on_predictions_file = "risearch_" + temp + ".out.gz"
            except:
                sys.stderr.write('#ERROR: failed running RIsearch2 on on-target\n#TERMINATING\n')
                sys.stderr.write("Unexpected error: "+str(sys.exc_info())+'\n')
                pipeline_exit(failed=True)
                

    if args.oa == '':
        if on_predictions_file == '' or args.p == '':
            sys.stderr.write('#CHECKPOINT: On-target opening energies are not relevant/needed\n')
        elif args.of == '':
            sys.stderr.write('#ERROR: On-target fasta file is not given.\n#TERMINATING\n')
            pipeline_exit(failed=True)
        else:
            acc_w = 80
            acc_l = 40
            if args.ap != '':
                for inf in args.ap.split(';'):
                    if inf[0] == 'W':
                        acc_w = int(inf[2:])
                    elif inf[0] == 'L':
                        acc_l = int(inf[2:])
                sys.stderr.write('#CHECKPOINT: RNAplfold parameters W='+acc_w+' L='+acc_l+'\n')
            on_target_name = ''
            with open(args.of, 'rt') as f:
                for line in f:
                    if line[0] == '>':
                        on_target_name = line.split()[0].rstrip()[1:]
                        break
            if on_target_name == '':
                sys.stderr.write('#ERROR: wrong on-target sequence file (no sequence name)\n#TERMINATING\n')
                pipeline_exit(failed=True)

            on_accessibility_file = on_target_name + '.open.acc.bin'
            # Run RNAplfold and create binary opening energies file
            try:
                create_acc_bin(args.of, on_target_name, on_accessibility_file, acc_w, acc_l)
            except:
                sys.stderr.write('#ERROR: failed while computing accessibility on on-target\n#TERMINATING\n')
                sys.stderr.write("Unexpected error: "+str(sys.exc_info())+'\n')
                pipeline_exit(failed=True)
                
    t2=time.time()
    sys.stderr.write('#TIMETRACK: On-target files preparation '+str(t2-t1)+'\n')

    # 2.2 Transcriptome intersection to get expression levels
    if args.t != '':
        try:
            off_predictions_file = get_transcriptome_intersect(off_predictions_file, args.t, args.feature, args.expmetric, args.type, sort=args.sort, keep=args.intersection)
        except:
            sys.stderr.write('#ERROR: failed while transcriptome intersection\n#TERMINATING\n')
            sys.stderr.write("Unexpected error: "+str(sys.exc_info())+'\n')
            pipeline_exit(failed=True)
        sys.stderr.write('#CHECKPOINT: Predictions and Transcriptome intersection is done and saved to ' +
                         off_predictions_file + ' ' + on_predictions_file + '\n')

    t3=time.time()
    sys.stderr.write('#TIMETRACK: Transcriptome Intersection '+str(t3-t2)+'\n')

    # Step 3: Get the predictions
    try:
        pres, on_targets, min_p, min_na = get_predictions(off_predictions_file, on_id=ont_id, on_chrom=ont_c,
                                                          on_spos=ont_s, on_epos=ont_e, on_strand=ont_str,
                                                          threshold=args.thr, thr_after_acc=args.thrAfterOpEn,
                                                          bins=args.p, less=args.less, chromosome_choice=args.chr, 
                                                          qid = args.q, ris2 = args.rx, seq_file = args.os)
    except:
        sys.stderr.write("#ERROR: Cannot finish reading the predictions\n")
        sys.stderr.write("Unexpected error: "+str(sys.exc_info())+'\n')
        pipeline_exit(failed=True)

    t4=time.time()
    sys.stderr.write('#TIMETRACK: Get Predictions '+str(t4-t3)+'\n')

    # Step 4: Perform the on target update
    if on_predictions_file != '':
        oexp = args.oexp if args.oexp>0.0 else args.oexp2
        try:
            min_p, min_na = add_update_on_targets(on_predictions_file, on_accessibility_file, pres, on_targets, min_p, min_na, oexp,
                                                 threshold=args.thr, thr_after_acc=args.thrAfterOpEn)
        except:
            sys.stderr.write("#ERROR: Cannot finish updating ontarget predictions\n")
            sys.stderr.write("Unexpected error: "+str(sys.exc_info())+'\n')
            pipeline_exit(failed=True)

    t5=time.time()
    sys.stderr.write('#TIMETRACK: On target Update '+str(t5-t4)+'\n')

    # Step 5: Save interactions and opening energies in a file
    if args.sa != '':
        with open(args.sa, 'wt') as out_f:
            for (chromosome, spos, epos, strand, energy, open_eng, expression, gid, tid) in pres:
                out_f.write('\t'.join([chromosome, str(spos), str(epos), strand, str(energy), str(open_eng), str(expression), gid, tid])+'\n')

    t6=time.time()
    sys.stderr.write('#TIMETRACK: Save predictions and opening energies '+str(t6-t5)+'\n')

    # Step 6: Calculate partition function and probabilities for each prediction
    # Step 6.1: Get alpha beta parameter settings
    ab_list = [(1.0,1.0)]
    for alpha in args.alpha.split(';'):
        if float(alpha) == 1.0:
            continue
        else:
            for gamma in args.gamma.split(';'):
                if float(alpha) <= float(gamma) and (float(alpha), float(gamma)) not in ab_list:
                    ab_list.append((float(alpha), float(gamma)))
    if args.theta != '':
        for theta in args.theta.split(';'):
           if float(theta) not in ab_list:
               ab_list.append(float(theta))

    # Step 6.2: Apply distance from minimum energy threshold if given
    if args.tp is not None:
        pres, on_targets = apply_tp_threshold(pres, on_targets, args.tp, min_p, min_na, args.thrAfterOpEn)

    t7=time.time()
    sys.stderr.write('#TIMETRACK: Percent threshold and alpha beta parsing '+str(t7-t6)+'\n')

    # Step 6.3: Get probabilities
    try:
        probs, probs_na, z_val, z_val_na, on_probs, on_probs_na = get_probabilities(pres, on_targets, min_p, min_na, ab_list, args.offPs)
    except:
        sys.stderr.write("#ERROR: Cannot finish computing the partition function and probabilities\n")
        sys.stderr.write("Unexpected error: "+str(sys.exc_info())+'\n')
        pipeline_exit(failed=True)

    t8=time.time()
    sys.stderr.write('#TIMETRACK: Compute Z and probabilities '+str(t8-t7)+'\n')

    # Step 7: Report the results
    try:
        fin_prob = report(pres, probs, probs_na, z_val, z_val_na, on_probs, on_probs_na, ab_list, args.q, args.o, args.offPs)
    except:
        sys.stderr.write("#ERROR: Cannot finish reporting the results\n")
        sys.stderr.write("Unexpected error: "+str(sys.exc_info())+'\n')
        pipeline_exit(failed=True)

    t9=time.time()
    sys.stderr.write('#TIMETRACK: Report results '+str(t9-t8)+'\n')

# Main function for the pipeline
if __name__ == '__main__':

    parser = get_parser()

    # STEP 1: Get the necessary arguments
    args = parser.parse_args()
    sys.stderr.write('#START: arguments parsed\n')
    for k in sorted(args.__dict__):
        if (args.__dict__[k] is not None) and args.__dict__[k] != '':
            sys.stderr.write('#ARG: ' + k + ' = ' + str(args.__dict__[k]) + '\n')

    run_pipeline(args)
    pipeline_exit()


