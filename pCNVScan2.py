"""
=========================================================================================================================
 CNVScan: Algorithm for searching Copy Number Variations in a sample file (*.pileup) among a list of genes in a ref file (*.gbk or *.gff).
 
 Beghain J et al. Plasmodium copy number variation scan: gene copy numbers evaluation in haploid genomes.
 Malar J. 2016 Apr 12;15:206. doi: 10.1186/s12936-016-1258-x.
 PMID: 27066902
 
 Required: 	- Python 3.5.3
			- BioPython [1]
			
 REQUIRED OPTIONS:
 Pattern size. Size of the patterns for algorithms. Recommended number: 6. Too High number will decrease program's speed and increase memory. Too low number reduce the sensitivity.
 Pileup sample file. Path to the input Pileup file.
 GeneBank reference file. Path to the input reference file.
 
 Version: 2.0.1 (sept. 2017)
 Author: Johann Beghain
 Contact: johann.beghain@inserm.fr
=========================================================================================================================
"""
#!/usr/bin/env python
import argparse
import ast
import logging
import multiprocessing as mp
import subprocess
import os
import sys
import builtins
#import numpy
from Bio import SeqIO
import pprint
from BCBio import GFF
from BCBio.GFF import GFFExaminer


###############################################################################
#############                CLASS                            #################
###############################################################################

#Definition of Chromosomes
class Chromosome:
    def __init__(self, name, length, seq):
        self.name = name
        self.seq = seq
        self.coverage = [0] * (length+1)
        self.length = length
        self.genes = {}
    def set_genes(self, genes):
        self.genes = genes
    def set_coverage(self, coverage):
        self.coverage = coverage

###############################################################################
#############                FUNCTIONS                            #############
###############################################################################

# Argument selection and help
def ArgumentsParser():
  parser=argparse.ArgumentParser(description = '''CNVScan: Algorithm for searching Copy Number Variations in a sample file (*.pileup) among a list of genes in a ref file (A gff+ fasta file OR a gbk file)''',
                                 epilog = '''The reads must be aligned on the same version of the genome as the gff or gbk+fasta provided.''')
  parser.add_argument('-p', '--pileup', type = str, required = True,
                      help = 'Path to a valid .pileup file. The alignment must have been done on the same reference as the reference file.')
  parser.add_argument('-o', '--output', type = str, required = True,
                      help = 'Path to the output file. Output a list of copy numbers by genes.')
  parser.add_argument('-i', '--genbank', type = str, required = False,
                      help = 'Path to a valid input .gbk reference file. No fasta file required.')
  parser.add_argument('-g', '--gff', type = str, required = False,
                      help = 'Path to a valid .gff reference file. If set you have to give also a fasta file in -f option.')
  parser.add_argument('-f', '--fasta', type = str, required = False,
                      help = 'Path to a valid .fasta file. Should be the same reference as the .gff file.')
  parser.add_argument('-k', '--kmer', type = str, required = False,
                      help = 'k-mer size for iterative search. Default: 6.')
  logger_args = parser.add_argument_group('Logger configuration')
  logger_args.add_argument( '--log_level',
                         type = str,
                         nargs = '?',
                         default = 'INFO',
                         help = 'log level',
                         choices = ['ERROR', 'error', 'WARNING', 'warning', 'INFO', 'info', 'DEBUG', 'debug'],
                         required = False )
  logger_args.add_argument( '--log_file',
                         type = str,
                         nargs = '?',
                         help = 'log file (use the stderr by default)',
                         required = False )
  return(parser.parse_args())

#Read a fasta and gff file for building a dict of Chromosome class
def extract_all_from_gff(gff_file_name, fasta_file_name):
    chromosomes = dict()
    #Record the sequence of chromosome from fasta
    for seq_record in SeqIO.parse(fasta_file_name, "fasta"):
        chromosomes[seq_record.id] = Chromosome(seq_record.id, len(seq_record), seq_record.seq)
    #Record the genes positions from gff
    gff_handle = open(gff_file_name)
    for rec in GFF.parse(gff_handle):
        print(rec.id)
        if rec.id not in chromosomes:
            logger.error(rec.id+" not in fasta. Please check the names of chromosomes: they must be the same in fasta and gff.")
            continue
        genes = dict()
        for gene in rec.features:
            #print(gene.id)
            for mRNA in gene.sub_features:
                if mRNA.type != "mRNA":
                    continue
                for exon in mRNA.sub_features:
                    #print(exon.type)
                    if exon.type == "exon":
                        #if rec.id not in genes:
                        #    genes[rec.id] = {}
                        if gene.id in genes:
                            genes[gene.id] = genes[gene.id]+'@'+str(exon.location.nofuzzy_start)+'_'+str(exon.location.nofuzzy_end)
                        else:
                            genes[gene.id] = str(exon.location.nofuzzy_start)+'_'+str(exon.location.nofuzzy_end)
        chromosomes[rec.id].set_genes(genes)
    gff_handle.close()
    return(chromosomes)

#Read a gbk file for building a dict of Chromosome class
def extract_all_from_gbk(gbk_file_name):
    record = SeqIO.read(gbk_file_name, "genbank")
    logger.info("Not implementend yet.")

#This function need a dict of Chromosome class.
#It is used after extract_all_from_gbk or either extract_all_from_gff
#Record coverage in all the genome
def extract_coverage_from_pileup(pileup_file_name, chromosomes):
    coverage = {}
    with open(pileup_file_name, 'r') as pileup:
        for line in pileup:
            #Record coverage
            record = line.split("\t")
            chromosomes[record[0]].coverage[int(record[1])]=int(record[3])
    return(chromosomes)

#This function first browse all the genes and compute a theorical ratio for all k-mers
def compute_ratios(chromosomes):
    kmers = dict()
    kmers_number = dict()
    for chromosome in chromosomes:
        #Get seq /!\ seq begin at 0 /!\
        sequence = list(chromosomes[chromosome].seq)
        coverage = chromosomes[chromosome].coverage
        #Get genes
        nb_genes = 0
        #print(chromosomes[chromosome].name)
        genes = chromosomes[chromosome].genes
        #print(chromosome)
        for gene in genes:
            #print(genes[gene])
            exons = genes[gene].split('@')
            nb_genes=nb_genes+1
            for exon in exons:
                list_exon = exon.split("_")
                if len(list_exon) > 2:
                    logger.warning("Error: end, start and another number in exon from "+gene+" Keeping just the 2 firsts.")
                    list_exon = list_exon[0:1]
                #For each nt in the exons, compute the mean coverage of kmers
                for i in range(int(min(list_exon)),int(max(list_exon))):
                    s=""
                    kmer_sequence = s.join(sequence[i:i+kmer_size])
                    kmer_coverage =  sum(float(embedding) for embedding in coverage[i:i+kmer_size]) / len(coverage[i:i+kmer_size])
                    #kmer_coverage = numpy.mean(coverage[i:i+kmer_size])
                    if kmer_sequence not in kmers:
                        kmers[kmer_sequence] = kmer_coverage
                        kmers_number[kmer_sequence] = 1
                    else:
                        # Average formula : _Xn= _Xn-1+(Xn-_Xn-1)/n
                        kmers_number[kmer_sequence] = kmers_number[kmer_sequence]+1
                        kmers[kmer_sequence]=kmers[kmer_sequence]+(kmer_coverage-kmers[kmer_sequence])/kmers_number[kmer_sequence]
    return(kmers)

def cnv_scan(chromosomes, theorical_ratios, output_name):
    output = open(output_name, 'w')
    for chromosome in chromosomes:
        sequence = list(chromosomes[chromosome].seq)
        coverage = chromosomes[chromosome].coverage
        genes = chromosomes[chromosome].genes
        for gene in genes:
            ratio = 0
            gene_size = 1
            exons = genes[gene].split('@')
            for exon in exons:
                list_exon = exon.split("_")
                if len(list_exon) > 2:
                    logger.warning("Error: end, start and another number in exon from "+gene+" Keeping just the 2 firsts.")
                    list_exon = list_exon[0:1]
                #For each nt in the exons, compute the coverage of kmers
                for i in range(int(min(list_exon)),int(max(list_exon))):
                    s=""
                    kmer_sequence = s.join(sequence[i:i+kmer_size])
                    kmer_coverage =  sum(float(embedding) for embedding in coverage[i:i+kmer_size]) / len(coverage[i:i+kmer_size])
                    #kmer_coverage = numpy.mean(coverage[i:i+kmer_size])
                    if kmer_sequence not in theorical_ratios:
                        logger.warning(kmer_sequence+" not in precomputed theorical ratios!")
                        continue
                    else:
                        #Ratio is the observed under theorical one
                        ratio = ratio+kmer_coverage/theorical_ratios[kmer_sequence]
                        gene_size=gene_size+1
            #Compute the mean ratio
            final_ratio=ratio/gene_size
            #print(final_ratio)
            output.write(gene+"\t"+str(round(final_ratio, 2))+"\n")

###############################################################################
#############                MAIN                                 #############
###############################################################################
if __name__ == "__main__":
    args=ArgumentsParser()
    # Logger config
    logging_std_format = '[%(levelname)s] %(message)s'
    logging_debug_format = '%(asctime)s [%(levelname)s] [%(threadName)s - %(name)s] %(message)s'
    log_level = args.log_level.upper()
    if ( log_level == 'DEBUG' ):
        logging_std_format = logging_debug_format
    logging_datefmt = '%d/%m/%Y - %H:%M:%S'
    if ( args.log_file != None ):
        logging.basicConfig( format = logging_std_format,
                             datefmt = logging_datefmt,
                             filename = args.log_file,
                             filemode = 'w',
                             level = log_level )
    else:
        logging.basicConfig( format = logging_std_format,
                             datefmt = logging_datefmt,
                             level = log_level )
    logger = logging.getLogger( 'main' )
    
    #Default number of k-mer
    kmer_size = 6
    if args.kmer:
        kmer_size = args.kmer
    chromosomes = dict()
    #Check input format for reference
    if args.genbank:
        logger.info("Input GBK file: "+args.genbank)
        extract_all_from_gbk(args.genbank)
    else:
        if args.gff:
            if args.fasta:
                logger.info("Input GFF file: "+args.gff+". Input FASTA file: "+args.fasta)
                chromosomes = extract_all_from_gff(args.gff, args.fasta)
            else:
                logger.error("Error: missing fasta reference input. With -g option you should have a valid -f fasta file.")
                sys.exit(1)
        else:
            logger.error("Error: missing reference input. Use either option -g + -f (gff + fasta file) or -i (gbk file)")
            sys.exit(1)
    #Parse pileup
    chromosomes = extract_coverage_from_pileup(args.pileup, chromosomes)
    #Compute ratios on all genome
    theorical_ratios = compute_ratios(chromosomes)
    #Compute ratios for each genes
    cnv_scan(chromosomes, theorical_ratios, args.output)
    logger.info("End Program")
