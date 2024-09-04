# Merge sample-level GTFs into one GTF

import pandas as pd
import os
import sys
import pyranges as pr

p = os.path.dirname(os.getcwd())+'/scripts/'
sys.path.append(p)

from sm_utils import *

include: '../snakemake/cerberus.smk'

configfile: 'config.yml'
config_tsv = 'config.tsv'
df = parse_config(config_tsv)
analysis = ['lyric']
# analysis = ['flair', 'iq', 'espresso']

# df = df.loc[df.tech_rep.isin(['GM10493_1',
#                               'GM12878_1',
#                               'GM22300_1',
#                               'GM19117_1'])]


wildcard_constraints:
    tech_rep='|'.join([re.escape(x) for x in df.tech_rep.tolist()]),
    lab_rep='|'.join([re.escape(x) for x in df.lab_rep.tolist()]),
    end_mode='|'.join([re.escape(x) for x in end_modes]),

rule all:
    input:
        expand(config['lyric']['gtf'],
             tech_rep=df.tech_rep.tolist())

rule gff_to_gff:
    input:
        gff = config['lyric']['gff']
    output:
        gtf = config['lyric']['gtf']
    run:
        df = pr.read_gff(input.gff)
        df.to_gtf(output.gtf)
