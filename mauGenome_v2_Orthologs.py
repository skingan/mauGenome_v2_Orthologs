#!/usr/bin/env python

# make dictionary of XLOCv2 - FBgn orthologs
# XLOCv2FBgn = {XLOCv2: {Transcript_Orth: FBgn, CDS_Orth: FBgn}}
import csv
file = '/home/LCPG/skingan/MauGenome/OrthologCalls/XLOCv2_FBgnMap.txt'
XLOCv2FBgn = {}
with open(file) as f:
    for line in csv.DictReader(f, delimiter='\t'):
        subdict = {'Transcript_Orth': line['Transcript_Orth'], 'CDS_Orth': line['CDS_Orth']}
        XLOCv2FBgn[line['XLOC']] = subdict


# make dictionary of TCONSv2 - TCONSv1 pairs
# TCONSv2TCONSv1 = {TCONSv2: TCONSv1}
file = '/home/LCPG/skingan/MauGenome/OrthologCalls/mauGenomeVersionOrthsReformated.txt'
TCONSv2TCONSv1 = {}
with open(file) as f:
    for line in csv.DictReader(f, delimiter='\t'):
        TCONSv2TCONSv1[line['#version1.1']] = line['version1.0']

# function to parse line of a gff file
# creates dictionary for 9 fields
def parse_gff(line):
    import csv
    dict = {}
    keys = ['chrom', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    for k,v in zip(keys,line):
        dict[k] = v
    return dict

# function to create dictionary of TCONS - XLOC pairs from gtf file
# all pairs in file are saved in memory
def TCONS_XLOC_dict(gtf):
    from re import search
    dict={}
    with open(gtf) as f:
        for line in csv.reader(f, delimiter='\t'):
            gff_dict = parse_gff(line)
            XLOC_pattern = r'(XLOC_[0-9]{6})'
            TCONS_pattern = r'(TCONS_[0-9]{8})'
            XLOC = ''
            TCONS = ''
            s = search(XLOC_pattern,gff_dict['attributes'])
            if s:
                XLOC = s.group(1)
            s = search(TCONS_pattern,gff_dict['attributes'])
            if s:
                TCONS = s.group(1)
            dict[TCONS] = XLOC
    return dict

# make dictionaries for v1 and v2 of gtf files
file = '/home/LCPG/skingan/MauGenome/OrthologCalls/v1.gtf'
TCONSv1XLOC = TCONS_XLOC_dict(file)
file = '/home/LCPG/skingan/MauGenome/OrthologCalls/v2.gtf'
TCONSv2XLOC = TCONS_XLOC_dict(file)


# function to create dictionary of XLOC - FBgn orthologs from gff file
# all pairs in file are saved in memory
# separate keys for transcript and CDS orthologs
def XLOC_FBgn(gff):
    from re import sub
    dict={}
    with open(gff) as f:
        for line in csv.reader(f, delimiter='\t'):
            gff_dict = parse_gff(line)
            if gff_dict['type'] == 'gene':
                XLOC = ''
                FBgn = ''
                transcript_FBgn = ''
                CDS_FBgn = '' 
                for a in gff_dict['attributes'].split(';'):         
                    if 'XLOC' in a:
                        XLOC = sub('ID=','',a)
                    if 'Dmel_ortholog' in a:
                        FBgn = sub('Dmel_ortholog=','',a)
                    if 'Dmel_transcript_ortholog' in a:
                        transcript_FBgn = sub('Dmel_transcript_ortholog=','',a)
                    if 'Dmel_CDS_ortholog' in a:
                        CDS_FBgn = sub('Dmel_CDS_ortholog=','',a)
                if FBgn != '':
                    dict[XLOC] =  {'transcript_FBgn': FBgn, 'CDS_FBgn': FBgn}
                else:
                    dict[XLOC] =  {'transcript_FBgn': 'NA', 'CDS_FBgn': 'NA'}
                if transcript_FBgn != '':
                    dict[XLOC]['transcript_FBgn'] = transcript_FBgn
                if CDS_FBgn != '':
                    dict[XLOC]['CDS_FBgn'] = CDS_FBgn
    return dict

file = '/home/LCPG/skingan/MauGenome/OrthologCalls/mau12v1.0.gff'            
XLOCv1FBgn = XLOC_FBgn(file)


# make dictionary of FBgn and best blast hit to TCONSv2
#file = '/home/LCPG/skingan/MauGenome/OrthologCalls/dmelTranscript2TCONSv2.blast'
#BestBlast = {}
#with open(file) as f:
#    for line in csv.reader(f, delimiter='\t'):
# add only first line with FBtr into dictionary
#        if line[0] in BestBlast:
#            continue
#        else:
#            BestBlast[line[0]] = {'TCONS': line[1], 'FBtr_length': line[2], 'TCONS_length': line[3], 'aln_length': line[4]}

# make dictionary of CDS and best blast hit to TCONSv2
file = '/home/LCPG/skingan/MauGenome/OrthologCalls/dmelCDS_vs_mau12v2ORF.blast'
BestBlast = {}
with open(file) as f:
    for line in csv.reader(f, delimiter='\t'):
# add only first line with FBtr into dictionary
        if line[0].replace('-','_') in BestBlast:
            continue
        else:
            BestBlast[line[0].replace('-','_')] = {'TCONS': line[1].split('|')[0], 'CDS_length': line[2], 'TCONS_length': line[3], 'aln_length': line[4]}


# make dictionary of FBtr and FBgn
#file = '/home/LCPG/skingan/MauGenome/OrthologCalls/dmel-all-transcript-r5.55.fasta'
#FBtrFBgn = {}
#from re import search
#with open(file) as f:
#    for line in f:
#        if search(r'^>',line):
#            FBtr = ''
#            FBgn = ''
#            s1 = search(r'(FBtr[0-9]{7})',line.split(';')[0])
#            if s1:
#                FBtr = s1.group(1)
#            s2 = search(r'parent=(FBgn[0-9]{7})',line)
#            if s2:
#                FBgn = s2.group(1)
#            FBtrFBgn[FBtr] = FBgn
           
# make dictionary of CDS and FBgn
file = '/home/LCPG/skingan/MauGenome/OrthologCalls/dmel-all-CDS-r5.55.fasta'
CDSFBgn = {}
from re import search
with open(file) as f:
    for line in f:
        if search(r'^>',line):
            CDS = ''
            FBgn = ''
            tmp = line.split(';')[0]
            tmp1 = tmp.replace(' type=CDS','')
            tmp2 = tmp1.replace('>','')
            CDS = tmp2.replace('-','_')
            s2 = search(r'parent=(FBgn[0-9]{7})',line)
            if s2:
                FBgn = s2.group(1)
            CDSFBgn[CDS] = FBgn


# start with version2 TCONS
for tcons2 in sorted(TCONSv2TCONSv1.keys()):
# convert to version 2 XLOC
    xloc2 = TCONSv2XLOC[tcons2]
# get FBgn orthologs
# route #1, from ortholog calls of TCONS new in version 2
    if xloc2:
        if xloc2 in XLOCv2FBgn:
            CDS_fbgn = XLOCv2FBgn[xloc2]['CDS_Orth']
# route #2, from ortholog calls in version 1
        else:
            tcons1 = TCONSv2TCONSv1[tcons2] # translate to TCONSv1
            if tcons1 in TCONSv1XLOC:
                if tcons1 != 'NA':
                    xloc1 = TCONSv1XLOC[tcons1] # get XLOC
                    CDS_fbgn = XLOCv1FBgn[xloc1]['CDS_FBgn'] # get transcript ortholog(s)
                else:
                    CDS_fbgn = 'NA'
            
    else:
        print "error! no xloc2!"
# now I have tcons2, xloc2, CDS_fbgn!
    output = [tcons2, xloc2, CDS_fbgn]
# find fbgn in FBtrFBgn dictionary by searching through values
    for k in CDSFBgn.keys():
        if CDSFBgn[k] == CDS_fbgn:
# set alnL as variable to choose best transcript to print
            alnL = 0
            add = ['NA', 'NA', 'NA', 'NA']
# find transcript in Blast dictionary by searching through keys
            if k in BestBlast:
                if BestBlast[k]['aln_length'] > alnL:
# from BestBlast dictionary get lengths
# convert TCONS to XLOC to check reciprocal blast
                    aln = BestBlast[k]['aln_length']
                    xloc_check = TCONSv2XLOC[BestBlast[k]['TCONS']]
                    add = [xloc_check,BestBlast[k]['CDS_length'],BestBlast[k]['TCONS_length'],BestBlast[k]['aln_length']]
    output.extend(add)
    if xloc_check == xloc2:
        output.append('match')
    else:
        output.append('error')
    print "\t".join(output)    
