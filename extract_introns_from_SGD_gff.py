#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 11:26:31 2020

@author: kevinroy
"""
import pandas as pd
GFF_infilename = '/u/home/k/kevinh97/project-guillom/COMPASS/S288C_reference_genome_R64-2-1_20150113/saccharomyces_cerevisiae_R64-2-1_20150113_no_seq.gff'
intron_outfilename = '/u/home/k/kevinh97/project-guillom/COMPASS/S288C_reference_genome_R64-2-1_20150113/saccharomyces_cerevisiae_R64-2-1_20150113_introns.tsv'
  
def get_attribute_field(attribute, field):
    '''
    input: annotation in gff format
    output: parent_ID
    
    chrII	SGD	gene	2907	5009	.	-	.	ID=YBL111C;Name=YBL111C;Ontology_term=GO:0003674,GO:0005739,GO:0008150;Note=Helicase-like%20protein%20encoded%20within%20the%20telomeric%20Y%27%20element%3B%20relocalizes%20from%20mitochondrion%20to%20cytoplasm%20upon%20DNA%20replication%20stress;display=Helicase-like%20protein%20encoded%20within%20the%20telomeric%20Y%27%20element;dbxref=SGD:S000002151;orf_classification=Verified
    chrII	SGD	CDS	2907	4116	.	-	1	Parent=YBL111C_mRNA;Name=YBL111C_CDS;orf_classification=Verified
    chrII	SGD	CDS	4216	5009	.	-	0	Parent=YBL111C_mRNA;Name=YBL111C_CDS;orf_classification=Verified
    chrII	SGD	intron	4117	4215	.	-	.	Parent=YBL111C_mRNA;Name=YBL111C_intron;orf_classification=Verified
    chrII	SGD	mRNA	2907	5009	.	-	.	ID=YBL111C_mRNA;Name=YBL111C_mRNA;Parent=YBL111C
    '''
    dist = len(field)
    if field not in attribute:
        return None
    for idx in range(len(attribute)):
        if attribute[idx:idx+dist] == field:
            id_start = idx+dist
            break
    for idx in range(id_start, len(attribute)):
        if attribute[idx] == ';':
            id_end = idx
            break
    return attribute[id_start:id_end]

def load_systematic_to_common_gene_name_dict(GFF_infilename):
    '''
    input: gff file with intron annotations
    output: dictionary of chromosomes of start and stop coordinate tuples of strand and gene name tuple
    '''
    systematic_to_common_gene_name_dict = {}
    with open(GFF_infilename,'r') as GFF:
        for line in GFF:
            if not line[0] == '#':
                break 
        for line in GFF:
            if line.startswith('###'):
                break
            info = line.strip().split('\t')
            chromosome, region, sequence_type, start, end, period, strand, frame, annotation = info
            if chromosome != 'chrMito' and chromosome != 'chrmt':  ## exlude mitochondrial introns
                start, end = int(start), int(end)
                if sequence_type == 'gene':
                    systematic_gene_name = get_attribute_field(annotation, 'ID=')
                    common_gene_name = get_attribute_field(annotation, 'gene=')
                    if common_gene_name == None:
                        common_gene_name = systematic_gene_name
                    systematic_to_common_gene_name_dict[systematic_gene_name] = common_gene_name
    GFF.close()
    return systematic_to_common_gene_name_dict

systematic_to_common_gene_name_dict = load_systematic_to_common_gene_name_dict(GFF_infilename)

intron_dict = {}
intron_colnames = 'chrom, start, end, strand, type, Name, gene, intron_type'.split(', ')
for colname in intron_colnames:
    intron_dict[colname] = []
  
with open(GFF_infilename, 'r') as GFF:
    for line in GFF:
        if not line[0] == '#':
            break 
    for line in GFF:
        if line.startswith('###'):
            break
        seqname, source, feature, start, end, score, strand, frame, attribute = line.strip().split()
        if 'intron' in feature:
            Name = get_attribute_field(attribute, 'Name=')
            Name = Name.replace('%28', '(')
            Name = Name.replace('%29', ')')
            gene = systematic_to_common_gene_name_dict.get( Name.replace('_intron', ''), Name)
            if seqname == 'chrmt':
                intron_type = 'mitochondrial_mRNA_intron'
            else:
                intron_type = 'spliceosomal_intron'
            intron_dict['chrom'] += [seqname]
            intron_dict['start'] += [start]
            intron_dict['end'] += [end]
            intron_dict['strand'] += [strand]
            intron_dict['type'] += [feature]
            intron_dict['Name'] += [Name]
            intron_dict['gene'] += [gene]
            intron_dict['intron_type'] += [intron_type]


    GFF.close()
        
intron_df = pd.DataFrame.from_dict(intron_dict)
intron_df.columns 

intron_df.to_csv(intron_outfilename, index = False, sep = '\t' )
        
        