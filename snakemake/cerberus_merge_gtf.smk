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

wildcard_constraints:
    tech_rep='|'.join([re.escape(x) for x in df.tech_rep.tolist()]),
    lab_rep='|'.join([re.escape(x) for x in df.lab_rep.tolist()]),
    end_mode='|'.join([re.escape(x) for x in end_modes]),


rule all:
    input:
        expand(config['cerberus']['agg']['ics']),
        expand(config['cerberus']['agg']['ends'],
               tech_rep=df.tech_rep.tolist(),
               end_mode=end_modes)

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
    resources:
        threads = 1,
        nodes = 2
    params:
        dist =  lambda wc: config['params']['cerberus'][wc.end_mode]['dist'],
        slack = lambda wc: config['params']['cerberus'][wc.end_mode]['slack']
    output:
        bed = config['ref']['cerberus']['ends']

use rule gtf_to_ic as cerb_gtf_to_ic with:
    input:
        gtf = lambda wc: expand(config['iq']['gtf'],
                                lab_rep=get_df_val(df,
                                         'lab_rep',
                                         wc.tech_rep,
                                         'tech_rep'))
    output:
        tsv = config['cerberus']['ics']

use rule agg_ics_cfg as cerb_agg_ics_cfg with:
    input:
        ref_tsv = rules.ref_gtf_to_ic.output.tsv,
        tsvs = expand(rules.cerb_gtf_to_ic.output.tsv,
                     tech_rep=df.tech_rep.tolist())
    params:
        ref_source = config['ref']['gtf_ver'],
        sources = df.tech_rep.tolist()
    output:
        cfg = config['cerberus']['agg']['ics_cfg']

use rule agg_ics as cerb_agg_ic with:
    input:
        cfg = rules.cerb_agg_ics_cfg.output.cfg,
        sample_tsvs = rules.cerb_agg_ics_cfg.input.tsvs,
        ref_tsv = rules.cerb_agg_ics_cfg.input.ref_tsv
    output:
        tsv = config['cerberus']['agg']['ics']

use rule gtf_to_ends as cerb_gtf_to_ends with:
    input:
        gtf = rules.cerb_gtf_to_ic.input.gtf
    params:
        dist =  lambda wc: config['params']['cerberus'][wc.end_mode]['dist'],
        slack = lambda wc: config['params']['cerberus'][wc.end_mode]['slack']
    output:
        bed = config['cerberus']['ends']

use rule agg_ends_cfg as cerb_agg_config with:
    input:
        ref_ends = rules.ref_gtf_to_ends.output.bed,
        sample_ends = lambda wc: expand(rules.cerb_gtf_to_ends.output.bed,
                    tech_rep=df.tech_rep.tolist(),
                    end_mode=wc.end_mode)
    params:
        add_ends = [True]+[True for i in range(len(df.tech_rep.tolist()))],
        refs = [True]+[False for i in range(len(df.tech_rep.tolist()))],
        ref_source = config['ref']['gtf_ver'],
        sample_sources = df.tech_rep.tolist()
    output:
        cfg = config['cerberus']['agg']['ends_cfg']

use rule agg_ends as cerb_agg_ends with:
    input:
        ref_ends = rules.cerb_agg_config.input.ref_ends,
        sample_ends = rules.cerb_agg_config.input.sample_ends,
        cfg = rules.cerb_agg_config.output.cfg
    params:
        slack = lambda wc: config['params']['cerberus'][wc.end_mode]['agg_slack']
    output:
        agg_ends = config['cerberus']['agg']['ends']


# use rule agg_ends_cfg as cerb_agg_tss_config with:
#     input:
#         ref_ends = lambda wc: expand(rules.ref_gtf_to_ends.output.bed,
#                     end_mode='tss'),
#         sample_ends = lambda wc: expand(rules.cerb_gtf_to_ends.output.bed,
#                     tech_rep=df.tech_rep.tolist(),
#                     end_mode='tss'),
#     params:
#         add_ends = [True]+[True for i in range(len(df.tech_rep.tolist()))],
#         refs = [True]+[False for i in range(len(df.tech_rep.tolist()))],
#         ref_source = config['ref']['gtf_ver'],
#         sample_sources = df.tech_rep.tolist()
#     output:
#         cfg = expand(config['cerberus']['agg']['ends_cfg'],
#                      end_mode='tss')[0]
#
# use rule agg_ends_cfg as cerb_agg_tes_config with:
#     input:
#         ref_ends = lambda wc: expand(rules.ref_gtf_to_ends.output.bed,
#                     end_mode='tes'),
#         sample_ends = lambda wc: expand(rules.cerb_gtf_to_ends.output.bed,
#                     tech_rep=df.tech_rep.tolist(),
#                     end_mode='tes')
#     params:
#         add_ends = [True]+[True for i in range(len(df.tech_rep.tolist()))],
#         refs = [True]+[False for i in range(len(df.tech_rep.tolist()))],
#         ref_source = config['ref']['gtf_ver'],
#         sample_sources = df.tech_rep.tolist()
#     output:
#         cfg = expand(config['cerberus']['agg']['ends_cfg'],
#                      end_mode='tes')[0]
#
# use rule agg_ends as cerb_agg_tss with:
#     input:
#         ref_ends = rules.cerb_agg_tss_config.input.ref_ends,
#         sample_ends = rules.cerb_agg_tss_config.input.sample_ends,
#         cfg = rules.cerb_agg_tss_config.output.cfg
#     params:
#         slack = lambda wc: config['params']['cerberus']['tss']['agg_slack']
#     output:
#         agg_ends = expand(config['cerberus']['agg']['ends'], end_mode='tss')[0]
#
# use rule agg_ends as cerb_agg_tes with:
#     input:
#         ref_ends = rules.cerb_agg_tes_config.input.ref_ends,
#         sample_ends = rules.cerb_agg_tes_config.input.sample_ends,
#         cfg = rules.cerb_agg_tes_config.output.cfg
#     params:
#         slack = config['params']['cerberus']['tes']['agg_slack']
#     output:
#         agg_ends = expand(config['cerberus']['agg']['ends'], end_mode='tes')[0]
