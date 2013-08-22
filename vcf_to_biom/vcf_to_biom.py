#!/usr/bin/env python
# File created on 11 Jul 2013
from __future__ import division

__author__ = "John Chase"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["John Chase"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "John Chase"
__email__ = "chasejohnh@gmail.com"
__status__ = "Development"

from biom.table import table_factory, SparseOTUTable
from biom.parse import parse_mapping, generatedby
from numpy import array
import os

def process_data_entry_line (line, ids):
    fields = line.split('\t') 
    chr_n = fields[0] 
    pos_n = fields[1]
    rs_n = fields[2]
    ref = fields[3] 
    alt = fields[4]
    indiv_ids = zip(ids, indiv_snp_variation(line))
    return (chr_n, pos_n, rs_n, ref, alt, indiv_ids)

#process the information line. This line contains all of the information about the individuals. 
#The function returns a list of the individuals 
def process_header_entry(line):
    fields = line.split() 
    return fields [9:]

#Creates a list of allelic variation for each SNP. The list itself contains one entry of two numbers for each individual. 
def indiv_snp_variation(line): 
    genotypes = line.split('\t')[9:]
    formatted_genotypes = [] 
    for i in genotypes: 
        formatted_genotypes.append(map(int,i.split(':')[0].split('|')))
    return formatted_genotypes
 
def create_biom_table(f):
    sample_ids = None
    observation_ids = []
    observation_md = []
    data = []
    for line in f:
        if line.startswith('##'):
            pass
        elif line.startswith('#CHROM'): 
            sample_ids = process_header_entry(line)
            sample_md = []
            for i in sample_ids:
                sample_md.append({})
        else:
            if sample_ids == None:
                raise ValueError, "Didn't find '#CHROM' line before data lines. Can't continue."
            else: 
                chr_n, pos_n, rs_n, ref, alt, indiv_ids = \
                process_data_entry_line(line, sample_ids)
                for i in ref.split(','):
                    if len(i) > 1: 
                        check = True
                    else: 
                        check = False
                if check == True:
                    pass
                else: 
                    observation_id = "%s:%s" %(chr_n, pos_n)
                    if observation_id in observation_ids:
                        pass
                    else:    
                        observation_ids.append(observation_id)
                        meta_dic = {"alleles":(ref, alt),"rs":rs_n}
                        observation_md.append(meta_dic)
                        data_row = []
                        for indiv, variation in indiv_ids:
                            if variation == [0, 0]:
                                data_row.append(0)
                            elif variation == [0, 1]:
                                data_row.append(1)
                            elif variation == [1, 0]: 
                                data_row.append(1)                        
                            elif variation == [1, 1]:
                                data_row.append(2)
                        data.append(data_row)
    if len(data) == 0:
        raise ValueError, """No valid SNP data was present in the file. Indels will be 
ignored"""
    else: 
        data = array(data)
    return data, sample_ids, observation_ids, sample_md, observation_md
    
def create_biom_file(vcf_fp, output_fp, mapping_fp=None):
    vcf_f = open(vcf_fp, 'U')
    data, sample_ids, observation_ids, sample_md, observation_md =\
    create_biom_table(vcf_f)
    biom_table = table_factory(data, 
                              sample_ids, 
                              observation_ids,
                              sample_md, 
                              observation_md,
                              constructor=SparseOTUTable)
    if mapping_fp != None:
        mapping_f = parse_mapping(open(mapping_fp,'U'))
        biom_table.addSampleMetadata(mapping_f)
    output_f = open(output_fp, 'w')
    output_f.write(biom_table.getBiomFormatJsonString(generatedby()))
    output_f.close()
    
# biom biom.util.biom_open 
# and look for python module that can take zipped files.
# parrallel merge otu tables
# bl2seek on the command 


                        