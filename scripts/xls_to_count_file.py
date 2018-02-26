#!/usr/bin/env python

#Load required modules
import sys, os
import pandas as pd
import common
import logging

logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
#handler.setFormatter(NiceFormatter())
root = logging.getLogger()
root.addHandler(handler)
root.setLevel(logging.INFO)

logger.info('Reading file `%s`...', snakemake.input[0])
df = pd.read_excel(snakemake.input[0], header=3)

col_names = list(df.columns.values)


# ----------------- untreated vs. stimulated -----------------
logger.info('Writing file `count_files/untreated_vs_stimulated.txt`...')
untreated_vs_stimulated = open('count_files/untreated_vs_stimulated.txt', 'w')

untreated_vs_stimulated.write("GeneID\tuntreated_1\tuntreated_2\tuntreated_3\tstimulated_1\tstimulated_2\tstimulated_3\n")
for index, row in df.iterrows():
    untreated_vs_stimulated.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (row['Geneid'], row['U373minusdox1_ATCACG_L001_R1_001.sorted.bam'], row['U373minusdox2_TTAGGC_L001_R1_001.sorted.bam'], row['U373minusdox3_ACTTGA_L001_R1_001.sorted.bam'],
	row['U373plusdox1_GATCAG_L001_R1_001.sorted.bam'], row['U373plusdox2_TAGCTT_L001_R1_001.sorted.bam'], row['U373plusdox3_GGCTAC_L001_R1_001.sorted.bam']))

# ------------------ stimulated versus US28 ------------------
logger.info('Writing file `count_files/stimulated_vs_US28.txt`...')
stimulated_vs_US28 = open('count_files/stimulated_vs_US28.txt', 'w')

stimulated_vs_US28.write("GeneID\tstimulated_1\tstimulated_2\tstimulated_3\tUS28_1\tUS28_2\tUS28_3\n")
for index, row in df.iterrows():
    stimulated_vs_US28.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (row['Geneid'], row['U373plusdox1_GATCAG_L001_R1_001.sorted.bam'], row['U373plusdox2_TAGCTT_L001_R1_001.sorted.bam'], row['U373plusdox3_GGCTAC_L001_R1_001.sorted.bam'], row['U373iUS28plusdox1_GTGGCC_L001_R1_001.sorted.bam'], row['U373iUS28plusdox2_GTTTCG_L001_R1_001.sorted.bam'], row['U373iUS28plusdox3_CGTACG_L001_R1_001.sorted.bam']))

# ------------------ stimulated versus UL33 ------------------
logger.info('Writing file `count_files/stimulated_vs_UL33.txt`...')
stimulated_vs_UL33 = open('count_files/stimulated_vs_UL33.txt', 'w')

stimulated_vs_UL33.write("GeneID\tstimulated_1\tstimulated_2\tstimulated_3\tUL33_1\tUL33_2\tUL33_3\n")
for index, row in df.iterrows():
    stimulated_vs_UL33.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (row['Geneid'], row['U373plusdox1_GATCAG_L001_R1_001.sorted.bam'], row['U373plusdox2_TAGCTT_L001_R1_001.sorted.bam'], row['U373plusdox3_GGCTAC_L001_R1_001.sorted.bam'], row['U373iUL33plusdox1_GAGTGG_L001_R1_001.sorted.bam'], row['U373iUL33plusdox2_ACTGAT_L001_R1_001.sorted.bam'], row['U373iUL33plusdox3_ATTCCT_L001_R1_001.sorted.bam']))
