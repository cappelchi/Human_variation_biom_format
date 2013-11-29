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


from shutil import rmtree
from os.path import exists, join
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files, create_dir
from qiime.util import (get_qiime_temp_dir, 
                        get_tmp_filename)
from qiime.test import initiate_timeout, disable_timeout
from numpy import array


from vcf_to_biom import (process_data_entry_line, process_header_entry, 
indiv_snp_variation, create_biom_table)

class ExampleTests(TestCase):
    """ Tests the generation of lists from vcf files. """
    
    def setUp(self):
        """ Initialize variables to be used by the tests """
        self.example_file1 = example_file1 
        self.example_file2 = example_file2
        self.example_file3 = example_file3
        self.example_file4 = example_file4
        self.example_file5 = example_file5
        
##Tests for process_data_entry_line
    def test_process_data_entry_line(self):
        """Does this work when given correct input?"""
        ids = ['HG00096', 'HG00097', 'HG00099', 'HG00100', 'HG00101']
        line = '10	89623323	rs1044322	G	A	100	PASS	.	GT:AP	0|0:0.015,0.000	0|0:0.000,0.000	0|0:0.002,0.000	0|0:0.000,0.052	0|0:0.000,0.000'
        expected = ('10', '89623323', 'rs1044322', 'G', 'A', [('HG00096',[0,0]), ('HG00097', [0,0]), ('HG00099', [0,0]), ('HG00100', [0,0]), ('HG00101', [0,0])])
        self.assertEqual(process_data_entry_line(line, ids), expected)


    def test_process_data_entry_line_2(self):
        """Does this work when given correct input?"""
        ids = ['HG00096', 'HG00097', 'HG00099', 'HG00100', 'HG00101', 'HG00102', 'HG00103', 'HG00104', 'HG00106', 'HG00108', 'HG00109', 'HG00110']
        line = '10	89674997	rs116819638	A	G	100	PASS	ERATE=0.0004;AN=2184;AC=13;VT=SNP;RSQ=0.9179;AA=A;AVGPOST=0.9989;LDAF=0.0063;SNPSOURCE=LOWCOV;THETA=0.0013;AF=0.01;AFR_AF=0.03	GT:DS:GL	0|0:0.000:-0.03,-1.23,-5.00	0|0:0.000:-0.02,-1.33,-5.00	0|0:0.000:-0.02,-1.33,-5.00	0|0:0.000:-0.03,-1.21,-5.00	0|0:0.000:-0.06,-0.87,-5.00	0|0:0.000:-0.48,-0.48,-0.48	0|0:0.000:-0.10,-0.69,-4.70	0|0:0.000:-0.03,-1.22,-5.00	0|0:0.000:-0.01,-1.49,-5.00	0|0:0.000:-0.00,-2.05,-5.00	0|0:0.000:-0.06,-0.87,-5.00	0|0:0.000:-0.06,-0.90,-5.00'
        expected = ('10', '89674997', 'rs116819638', 'A', 'G', [('HG00096', [0, 0]), ('HG00097', [0, 0]), ('HG00099', [0, 0]), ('HG00100', [0, 0]), ('HG00101', [0, 0]), ('HG00102', [0, 0]), ('HG00103', [0, 0]), ('HG00104', [0, 0]), ('HG00106', [0, 0]), ('HG00108', [0, 0]), ('HG00109', [0, 0]), ('HG00110', [0, 0])])
        self.assertEqual(process_data_entry_line(line, ids), expected)

    def test_process_data_entry_line_3(self):
        """Does this work when given correct input?"""
        ids = ['HG00096', 'HG00097', 'HG00099', 'HG00100', 'HG00101', 'HG00102', 'HG00103', 'HG00104', 'HG00106', 'HG00108', 'HG00109', 'HG00110']
        line = '10	89674997	rs116819638	A	TGC,GA	100	PASS	ERATE=0.0004;AN=2184;AC=13;VT=SNP;RSQ=0.9179;AA=A;AVGPOST=0.9989;LDAF=0.0063;SNPSOURCE=LOWCOV;THETA=0.0013;AF=0.01;AFR_AF=0.03	GT:DS:GL	0|0:0.000:-0.03,-1.23,-5.00	0|0:0.000:-0.02,-1.33,-5.00	0|0:0.000:-0.02,-1.33,-5.00	0|0:0.000:-0.03,-1.21,-5.00	0|0:0.000:-0.06,-0.87,-5.00	0|0:0.000:-0.48,-0.48,-0.48	0|0:0.000:-0.10,-0.69,-4.70	0|0:0.000:-0.03,-1.22,-5.00	0|0:0.000:-0.01,-1.49,-5.00	0|0:0.000:-0.00,-2.05,-5.00	0|0:0.000:-0.06,-0.87,-5.00	0|0:0.000:-0.06,-0.90,-5.00'
        expected = ('10', '89674997', 'rs116819638', 'A', 'TGC,GA', [('HG00096', [0, 0]), ('HG00097', [0, 0]), ('HG00099', [0, 0]), ('HG00100', [0, 0]), ('HG00101', [0, 0]), ('HG00102', [0, 0]), ('HG00103', [0, 0]), ('HG00104', [0, 0]), ('HG00106', [0, 0]), ('HG00108', [0, 0]), ('HG00109', [0, 0]), ('HG00110', [0, 0])])
        self.assertEqual(process_data_entry_line(line, ids), expected)
        
#Tests for process_header_entry##
    def test_process_header_entry(self):
        """ Does this return a list of the fields when given the correct input?"""
        input = '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG00096	HG00097	HG00099	HG00100	HG00101'
        expected = ['HG00096', 'HG00097', 'HG00099', 'HG00100', 'HG00101']
        self.assertEqual(process_header_entry(input), expected)

##Test that indiv_snp_variation returns a list of the allele variation information. 
    def test_indiv_snp_variation(self): 
        """Does the function return correct output when given correct output"""
        input  = '10	89674917	rs182708158	T	G	100	PASS	ERATE=0.0004;LDAF=0.0024;AA=T;AN=2184;THETA=0.0008;VT=SNP;AVGPOST=0.9993;RSQ=0.9008;AC=5;SNPSOURCE=LOWCOV;AF=0.0023;ASN_AF=0.01	GT:DS:GL	0|0:0.000:-0.03,-1.22,-5.00	0|0:0.000:-0.01,-1.62,-5.00	0|0:0.000:-0.01,-1.87,-5.00	0|0:0.000:-0.18,-0.48,-2.19	0|0:0.000:-0.12,-0.62,-3.85	0|0:0.000:-0.14,-0.56,-3.22	0|0:0.000:-0.03,-1.15,-5.00	0|0:0.000:-0.01,-1.74,-5.00	0|0:0.000:-0.03,-1.19,-5.00	0|0:0.000:-0.05,-0.99,-5.00	0|0:0.000:-0.08,-0.75,-5.00	0|0:0.000:-0.02,-1.28,-5.00	0|0:0.000:-0.00,-2.07,-5.00	0|0:0.000:-0.01,-1.54,-5.00	0|0:0.000:-0.03,-1.15,-5.00	0|0:0.000:-0.03,-1.20,-5.00	0|0:0.000:-0.06,-0.91,-5.00	0|0:0.000:-0.00,-1.97,-5.00	0|0:0.000:-0.01,-1.51,-5.00'
        expected = [[0, 0], [0, 0], [0, 0], [0, 0], [0, 0],[0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0],[0, 0], [0, 0], [0, 0], [0, 0], [0, 0]]
        self.assertEqual(indiv_snp_variation(input), expected)
       
    def test_indiv_snp_variation_2(self): 
        """Does the function return correct output when given correct output"""
        input  = '10	89623323	rs1044322	G	A	100	PASS	.	GT:AP	0|0:0.015,0.000	0|0:0.000,0.000	0|0:0.002,0.000	0|0:0.000,0.052	0|0:0.000,0.000'
        expected = [[0, 0], [0, 0], [0, 0], [0, 0], [0, 0]]
        self.assertEqual(indiv_snp_variation(input), expected)
        
    def test_indiv_snp_variation_3(self): 
        """Does the function return correct output when given correct output"""
        input  = '10	89623323	rs1044322	G	A	100	PASS	.	GT:AP	1|0:0.015,0.000	0|0:0.000,0.000	0|0:0.002,0.000	1|1:0.000,0.052	0|0:0.000,0.000'
        expected = [[1, 0], [0, 0], [0, 0], [1, 1], [0, 0]]
        self.assertEqual(indiv_snp_variation(input), expected)
        
#Test for create_biom_table
    def test_create_biom_table(self):
        """Does the function return the correct output when given the correct input?"""
        vcf = list(create_biom_table(self.example_file1))
        expected = [array([[0, 0, 0, 0, 0]]), ['HG00096', 'HG00097', 'HG00099', 'HG00100', 'HG00101'], ['10:89623323'], [{}, {}, {}, {}, {}], [{'alleles': ('G', 'A'), 'rs': 'rs1044322'}]]
        self.assertEqual(vcf, expected)
 
    def test_create_biom_table_2(self):
        """Does the function return the correct output when given the correct input?"""
        vcf = list(create_biom_table(self.example_file2))
        expected = [array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
       [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]), ['HG00096', 'HG00097', 'HG00099', 'HG00100', 'HG00101', 'HG00102', 'HG00103', 'HG00104', 'HG00106', 'HG00108', 'HG00109', 'HG00110'], ['10:89674917', '10:89674997'], [{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}], [{'alleles': ('T', 'G'), 'rs': 'rs182708158'}, {'alleles': ('A', 'G'), 'rs': 'rs116819638'}]]
        self.assertEqual(vcf, expected)
        
    def test_create_biom_table_3(self):
        """Does the function return the correct output when given the correct input?"""
        self.assertRaises(ValueError, create_biom_table, self.example_file3)
        
    def test_create_biom_table_4(self):
        """Does the function return the correct output when given the correct input?"""
        vcf = list(create_biom_table(self.example_file4))
        expected = [array([[0, 0,],[0, 0,],[0, 0,],[0, 1]]), ['HG00096', 'HG00097'], ['10:89674917', '10:89674997', '10:89675036', '10:89675296' ], [{}, {}], [{'alleles': ('T', 'G'), 'rs': 'rs182708158'}, {'alleles': ('A', 'G'), 'rs': 'rs116819638'}, {'alleles': ('C', 'T'), 'rs': 'rs111627758'}, {'alleles': ('A', 'G'), 'rs': 'rs1234224'}]]
        self.assertEqual(vcf, expected)
        
    def test_create_biom_table_5(self):
        """Does the function return the correct output when given the correct input?"""
        vcf = list(create_biom_table(self.example_file5))
        expected = [array([[2],[1],[1],[1],[1],[1],[1]]), ['TNT028'], ['1:25611035', '1:25627613', '1:25627627', '1:25627628', '1:25656531', '1:25656673', '1:25688901'], [{}], [{'alleles': ('G', 'C'), 'rs': 'rs2301153'}, {'alleles': ('C', 'A'), 'rs': '.'}, {'alleles': ('C', 'CG'), 'rs': '.'}, {'alleles': ('A', 'C'), 'rs': '.'}, {'alleles': ('T', 'C'), 'rs': '.'}, {'alleles': ('C', 'T'), 'rs': 'rs667179'}, {'alleles': ('T', 'G'), 'rs': '.'}]]
        self.assertEqual(vcf, expected)
    
example_file1 = """##fileformat=VCFv4.0													
##source=BCM:SNPTools:hapfuse													
##reference=1000Genomes-NCBI37													
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">													
##FORMAT=<ID=AP,Number=2,Type=Float,Description="Allelic Probability, P(Allele=1|Haplotype)">													
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG00096	HG00097	HG00099	HG00100	HG00101
10	89623323	rs1044322	G	A	100	PASS	.	GT:AP	0|0:0.015,0.000	0|0:0.000,0.000	0|0:0.002,0.000	0|0:0.000,0.052	0|0:0.000,0.000""".split('\n')
        
        
example_file2 = """##reference=GRCh37																				
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG00096	HG00097	HG00099	HG00100	HG00101	HG00102	HG00103	HG00104	HG00106	HG00108	HG00109	HG00110
10	89674917	rs182708158	T	G	100	PASS	ERATE=0.0004;LDAF=0.0024;AA=T;AN=2184;THETA=0.0008;VT=SNP;AVGPOST=0.9993;RSQ=0.9008;AC=5;SNPSOURCE=LOWCOV;AF=0.0023;ASN_AF=0.01	GT:DS:GL	0|0:0.000:-0.03,-1.22,-5.00	0|0:0.000:-0.01,-1.62,-5.00	0|0:0.000:-0.01,-1.87,-5.00	0|0:0.000:-0.18,-0.48,-2.19	0|0:0.000:-0.12,-0.62,-3.85	0|0:0.000:-0.14,-0.56,-3.22	0|0:0.000:-0.03,-1.15,-5.00	0|0:0.000:-0.01,-1.74,-5.00	0|0:0.000:-0.03,-1.19,-5.00	0|0:0.000:-0.05,-0.99,-5.00	0|0:0.000:-0.08,-0.75,-5.00	0|0:0.000:-0.02,-1.28,-5.00
10	89674997	rs116819638	A	G	100	PASS	ERATE=0.0004;AN=2184;AC=13;VT=SNP;RSQ=0.9179;AA=A;AVGPOST=0.9989;LDAF=0.0063;SNPSOURCE=LOWCOV;THETA=0.0013;AF=0.01;AFR_AF=0.03	GT:DS:GL	0|0:0.000:-0.03,-1.23,-5.00	0|0:0.000:-0.02,-1.33,-5.00	0|0:0.000:-0.02,-1.33,-5.00	0|0:0.000:-0.03,-1.21,-5.00	0|0:0.000:-0.06,-0.87,-5.00	0|0:0.000:-0.48,-0.48,-0.48	0|0:0.000:-0.10,-0.69,-4.70	0|0:0.000:-0.03,-1.22,-5.00	0|0:0.000:-0.01,-1.49,-5.00	0|0:0.000:-0.00,-2.05,-5.00	0|0:0.000:-0.06,-0.87,-5.00	0|0:0.000:-0.06,-0.90,-5.00""".split('\n')

example_file3 = """##fileformat=VCFv4.0													
##source=BCM:SNPTools:hapfuse													
##reference=1000Genomes-NCBI37													
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">													
##FORMAT=<ID=AP,Number=2,Type=Float,Description="Allelic Probability, P(Allele=1|Haplotype)">													
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG00096	HG00097	HG00099	HG00100	HG00101
10	89623323	rs1044322	GT	A,TC	100	PASS	.	GT:AP	2|0:0.015,0.000	1|2:0.000,0.000	0|0:0.002,0.000	0|0:0.000,0.052	1|1:0.000,0.000""".split('\n')

example_file4 = """##fileformat=VCFv4.1										
##INFO=<ID=LDAF,Number=1,Type=Float,Description=""MLE Allele Frequency Accounting for LD"">										
##INFO=<ID=AVGPOST,Number=1,Type=Float,Description=""Average posterior probability from MaCH/Thunder"">										
##INFO=<ID=RSQ,Number=1,Type=Float,Description=""Genotype imputation quality from MaCH/Thunder"">								
##INFO=<ID=ERATE,Number=1,Type=Float,Description=""Per-marker Mutation rate from MaCH/Thunder"">										
##INFO=<ID=THETA,Number=1,Type=Float,Description=""Per-marker Transition rate from MaCH/Thunder"">										
##INFO=<ID=CIEND,Number=2,Type=Integer,Description=""Confidence interval around END for imprecise variants"">										
##INFO=<ID=SNPSOURCE,Number=.,Type=String,Description=""indicates if a snp was called when analysing the low coverage or exome alignment data"">									
##reference=GRCh37										
##reference=GRCh37										
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG00096	HG00097
10	89674917	rs182708158	T	G	100	PASS	ERATE=0.0004;LDAF=0.0024;AA=T;AN=2184;THETA=0.0008;VT=SNP;AVGPOST=0.9993;RSQ=0.9008;AC=5;SNPSOURCE=LOWCOV;AF=0.0023;ASN_AF=0.01	GT:DS:GL	0|0:0.000:-0.03,-1.22,-5.00	0|0:0.000:-0.01,-1.62,-5.00
10	89674997	rs116819638	A	G	100	PASS	ERATE=0.0004;AN=2184;AC=13;VT=SNP;RSQ=0.9179;AA=A;AVGPOST=0.9989;LDAF=0.0063;SNPSOURCE=LOWCOV;THETA=0.0013;AF=0.01;AFR_AF=0.03	GT:DS:GL	0|0:0.000:-0.03,-1.23,-5.00	0|0:0.000:-0.02,-1.33,-5.00
10	89674997	rs151009112	T	C	100	PASS	ERATE=0.0005;AA=T;AN=2184;LDAF=0.0018;THETA=0.0005;VT=SNP;RSQ=0.4609;SNPSOURCE=LOWCOV;AC=2;AVGPOST=0.9974;AF=0.0009;ASN_AF=0.0035	GT:DS:GL	0|0:0.000:-0.03,-1.23,-5.00	0|0:0.000:-0.04,-1.05,-5.00
10	89675036	rs111627758	C	T	100	PASS	ERATE=0.0004;LDAF=0.0047;AA=T;AN=2184;RSQ=0.8612;THETA=0.0005;VT=SNP;AC=9;SNPSOURCE=LOWCOV;AVGPOST=0.9984;AF=0.0041;AMR_AF=0.01;EUR_AF=0.01	GT:DS:GL	0|0:0.000:-0.03,-1.16,-5.00	0|0:0.000:-0.02,-1.28,-5.00
10	89675296	rs1234224	A	G	100	PASS	ERATE=0.0004;AC=915;RSQ=0.9931;THETA=0.0004;AA=G;AN=2184;VT=SNP;LDAF=0.4188;SNPSOURCE=LOWCOV;AVGPOST=0.9961;AF=0.42;ASN_AF=0.48;AMR_AF=0.43;AFR_AF=0.42;EUR_AF=0.37	GT:DS:GL	0|0:0.000:-0.01,-1.50,-5.0	1|0:1.000:-5.00,-0.00,-4.40""".split('\n')

example_file5 = """##reference=file://E:\Genomes\Homo_sapiens\UCSC\hg19\Sequence\WholeGenomeFASTA\genome.fa									
##source=GATK 1.6									
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TNT028
1	25611035	rs2301153	G	C	10465.01	PASS	AC=2;AF=1.0;AN=2;DP=470;QD=22.27;TI=NM_001127691,NM_016124;GI=RHD,RHD;FC=Silent,Silent	GT:AD:DP:GQ:PL:VF:GQX	1/1:0,470:470:99:10465,713,0:1.000:99
1	25627613	.	C	A	4033.01	PASS	AC=1;AF=0.5;AN=2;DP=167;QD=24.15;TI=NM_001127691,NM_016124;GI=RHD,RHD;FC=Silent,Silent	GT:AD:DP:GQ:PL:VF:GQX	0/1:55,112:167:99:4033,0,891:0.671:99
1	25627627	.	C	CG	3582	PASS	AC=1;AF=0.5;AN=2;DP=163;QD=21.98;TI=NM_001127691,NM_016124;GI=RHD,RHD;FC=Noncoding,Noncoding	GT:AD:DP:GQ:PL:VF:GQX	0/1:57,106:167:99:3582,0,1837:0.650:99
1	25627628	.	A	C	3555.01	PASS	AC=1;AF=0.5;AN=2;DP=163;QD=21.81;TI=NM_001127691,NM_016124;GI=RHD,RHD;FC=Silent,Silent	GT:AD:DP:GQ:PL:VF:GQX	0/1:57,106:163:99:3555,0,982:0.650:99
1	25656531	.	T	C	2414.01	PASS	AC=1;AF=0.5;AN=2;DP=532;QD=4.54;TI=NM_001127691,NM_016124;GI=RHD,RHD;FC=Silent,Silent;EXON	GT:AD:DP:GQ:PL:VF:GQX	0/1:387,145:532:99:2444,0,7968:0.273:99
1	25656673	rs667179	C	T	2233.01	PASS	AC=1;AF=0.5;AN=2;DP=626;QD=3.57;TI=NM_001127691,NM_016124;GI=RHD,RHD;FC=Silent,Silent;EXON	GT:AD:DP:GQ:PL:VF:GQX	0/1:483,143:629:99:2263,0,9941:0.228:99
1	25688901	.	T	G	11.34	LowVariantFreq	AC=1;AF=0.5;AN=2;DP=22;QD=0.52;TI=NM_138617,NM_020485,NM_138618,NM_138616;GI=RHCE,RHCE,RHCE,RHCE;FC=Silent,Silent,Silent,Silent;EXON	GT:AD:DP:GQ:PL:VF:GQX	0/1:20,2:22:40.83:41,0,442:0.091:11""".split('\n')



if __name__ == "__main__":
    main()
    
