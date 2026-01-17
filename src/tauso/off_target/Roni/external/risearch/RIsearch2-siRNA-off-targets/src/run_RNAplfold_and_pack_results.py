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

import os,sys,re,subprocess,re,string,random,gzip
import os.path

__author__ = "Ferhat Alkan: ferro@rth.dk"

if len(sys.argv)<4 or len(sys.argv)>5:
    print('python '+sys.argv[0]+' [input] [windowsize] [maxL] [-/+(strand info)]\n'+str(len(sys.argv))+' arg given, at least 4 needed')
    exit()

def id_generator(size=10, chars=string.ascii_uppercase + string.digits):
        return ''.join(random.choice(chars) for _ in range(size))

def saveplfoldscores(sequence,seqname,windowsize,maxl,ulength,strand):
    newws=windowsize
    newl=maxl
    infastafile=seqname+".fa.plFOLDtemp"
    newlen = len(sequence)
    if newlen<newws:
        newws=newlen
        newl=newlen
    input_f = open(infastafile,'wt')
    #tempid = id_generator()
    input_f.write('>'+seqname+'\n'+sequence+'\n')
    input_f.close()
    indscores=[]
    #run PLfold
    scorefile = seqname+'_openen'
    cmd = "RNAplfold -W "+str(newws)+" -L "+str(newl)+" -u "+str(ulength)+" -O < "+infastafile
    output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE ).communicate()
    if os.path.isfile(scorefile) and os.access(scorefile, os.R_OK):
        sf = open(scorefile,'r')
        for line in sf:
            if line[0]!='#' and line[0]!=' ' and line[0]!='':
                cols = line.rstrip().split()[1:]
                assert len(cols)==ulength
                for scstr in cols:
                    if scstr=='na' or scstr=='NA':
                        indscores.append(0.0)
                    else:
                        indscores.append(float(scstr))
        sf.close()
        try:
            os.remove(scorefile)
        except OSError as e:
            sys.stderr.write('ERRORCATCH delete temp RNAPLfold lunpfile\n') 
        try:
            os.remove(seqname+'_dp.ps')
        except OSError as e:
            sys.stderr.write('ERRORCATCH delete temp RNAPLfold psfile \n')
    else:
        sys.stderr.write('ERRORCATCH RNAPLfold do not give scores \n')
        indscores=[0.0]*newlen
    try:
        os.remove(infastafile)
    except OSError as e:
        sys.stderr.write('ERRORCATCH delete temp RNAPLfold inputfile')
    assert len(sequence)*ulength==len(indscores)
    if strand=='+':
        return indscores
    else:
        return indscores[::-1]

def complement(seq):
    c={'a':'t','c':'g','g':'c','t':'a','A':'T','C':'G','G':'C','T':'A','N':'N','n':'n','u':'a','U':'A'}
    l=len(seq)
    cseq=''
    for i in range(l):
        cur = seq[l-(i+1)]
        if cur in c:
            cur=c[cur]
        cseq+=cur
    return cseq

def checkgaps(seqs,len):
    newl=0
    newseqs={}
    for key in seqs.keys():
        newseqs[key]=''
    for i in range(0,len):
        check=False
        for key in seqs.keys():
            if seqs[key][i]!='-' and seqs[key][i]!='.':
                check=True
                break
        if check:
            newl+=1
            for key in seqs.keys():
                newseqs[key]+=seqs[key][i]
    return newseqs,newl

def runWithParts(infile,ws=int(sys.argv[2]),maxl=int(sys.argv[3]),ulength=30,strand='+'):
    in_f = gzip.open(infile,'r') if infile.endswith(".gz") else open(infile)
    seqchr=''
    sequence=''
    outf=None
    start=True
    end=False
    spos_i=0

    #Read aligned sequences
    for line in in_f:
        if line[0]=='>':
            if outf!=None and len(sequence)>0:
                namename='_'.join([seqchr,str(spos_i),str(spos_i+len(sequence))])
                sys.stderr.write(namename+'\n')
                scores = saveplfoldscores(sequence,namename,ws,maxl,ulength,strand)
                if not start:
                    scores=scores[ws*ulength:]
                for i in range(len(scores)):
                    number=int(round(scores[i]*10.0))
                    if number>255:
                        number=255
                    elif number<0:
                        number=0
                    outf.write(chr(number))
                outf.close()

            seqchr=line.split()[0].rstrip()[1:]
            outf=open(seqchr+'.open.acc.bin','wt+') if strand=='+' else open(seqchr+'.rev.open.acc.bin','wt+')

            start=True
            end=False
            spos_i=0
            sequence=''
        elif len(sequence)>=100000:
            namename='_'.join([seqchr,str(spos_i),str(spos_i+len(sequence))])
            sys.stderr.write(namename+'\n')
            scores=saveplfoldscores(sequence,namename,ws,maxl,ulength,strand)
            spos_i=(spos_i+len(sequence))-2*ws
            if not start:
                scores=scores[ws*ulength:]
            start=False
            for i in range(len(scores)-(ws*ulength)):
                number=int(round(scores[i]*10.0))
                if number>255:
                    number=255
                elif number<0:
                    number=0
                outf.write(chr(number))
            if strand != '+':
                sequence=complement(line.rstrip())+sequence[:(2*ws)]
            else:
                sequence=sequence[-(2*ws):]+line.rstrip()
        else:
            if strand != '+':
                sequence=complement(line.rstrip())+sequence
            else:
                sequence+=line.rstrip()

    if outf!=None and len(sequence)>0:
        namename='_'.join([seqchr,str(spos_i),str(spos_i+len(sequence))])
        sys.stderr.write(namename+'\n')
        scores = saveplfoldscores(sequence,namename,ws,maxl,ulength,strand)
        if not start:
            scores=scores[ws*ulength:]
        for i in range(len(scores)):
            number=int(round(scores[i]*10.0))
            if number>255:
                number=255
            elif number<0:
                number=0
            outf.write(chr(number))
        outf.close()
    in_f.close()

# Main function
if __name__ == '__main__':
    if len(sys.argv)==4:
        runWithParts(sys.argv[1])
    elif len(sys.argv)==5:
        runWithParts(sys.argv[1],strand=sys.argv[4])



