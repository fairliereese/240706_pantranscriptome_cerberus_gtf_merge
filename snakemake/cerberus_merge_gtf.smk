# Merge sample-level GTFs into one GTF

import pandas as pd

p = os.getcwd()
sys.path.append(p)

from sm_utils import *

include: 'cerberus.smk'

configfile: 'snakemake/config.yml'
config_tsv = config['meta']['library']

df = parse_config(config_tsv)
df['analysis'] = 'cerb_masked_genomic'

rule all:
    input:
        rules.all_cerberus_merge_gtf_all.input

def get_df_val(df, col, col2_val, col2='sample'):
    temp = df.loc[df[col2]==col2_val]
    assert len(temp.index) == 1
    return temp[col].values[0]

# actual rule calls
use rule gtf_to_ic as ref_gtf_to_ic with:
    input:
        gtf = config['ref']['gtf']
    output:
        tsv = config['ref']['cerberus']['ics']

use rule gtf_to_ends as ref_gtf_to_ends with:
    input:
        gtf = rules.ref_gtf_to_ic.input.gtf
    params:
        extend = lambda wc: config['params']['cerberus'][wc.end_mode]['extend'],
        slack = lambda wc: config['params']['cerberus'][wc.end_mode]['slack']
    output:
        bed = config['cerberus']['ends']

use rule gtf_to_ic as cerb_gtf_to_ic with:
    input:
        gtf = lambda wc: expand(config['iq']['gtf'],
                                lab_rep=get_df_val(df,
                                         'lab_rep',
                                         wc.tech_rep,
                                         'tech_rep'))
    output:
        tsv = config['cerberus']['ics']

use rule agg_ics_cfg as cerb_agg_ic_cfg with:
    input:
        ref_tsv = rules.ref_gtf_to_ic.output.tsv,
        tsvs = expand(rules.cerb_gtf_to_ic.output.tsv,
                     tech_rep=df.tech_rep.tolist())
    params:
        ref_source = config['ref']['gtf_ver'],
        sources = df.tech_rep.tolist()
    output:
        cfg = config['cerberus']['agg']['ic_cfg']

use rule agg_ics as cerb_agg_ic with:
    input:
        cfg = rules.cerb_agg_ic_cfg.output.cfg
    output:
        tsv = config['cerberus']['agg']['ics']

use rule gtf_to_ends as cerb_gtf_to_ends with:
    input:
        gtf = rules.cerb_gtf_to_ic.input.gtf
    params:
        extend = lambda wc: config['params']['cerberus'][wc.end_mode]['extend'],
        slack = lambda wc: config['params']['cerberus'][wc.end_mode]['slack']
    output:
        bed = config['cerberus']['ends']

rule all_cerberus_merge_gtf:
    expand(config['cerberus']['ics'],
           tech_rep=tech_rep=df.tech_rep.tolist()),
    expand(config['cerberus']['ends'],
           tech_rep=tech_rep=df.tech_rep.tolist(),
           end_mode=end_modes)
