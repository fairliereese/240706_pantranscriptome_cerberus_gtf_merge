# Merge sample-level GTFs into one GTF

import pandas as pd

p = os.getcwd()
sys.path.append(p)

from sm_utils import *

include: 'cerberus.smk'

configfile: 'snakemake/config.yml'
# config_tsv = config['meta']['library']

# # settings for running  w/ different sets of data
# analysis = 'espresso_pseudomasked_genomic'
# tool = 'espresso'
# config_tsv = f'snakemake/config_{analysis}_expression.tsv'
# df = parse_config(config_tsv)
# df['analysis'] = analysis
# input_gtf = config[analysis]['gtf']

# analysis = 'pseudomasked_genomic_isoquant_guided'
# tool = 'iq'
# config_tsv = f'snakemake/config_{analysis}.tsv'
# df = parse_config(config_tsv)
# df['analysis'] = analysis
# input_gtf = config[analysis]['gtf']

analysis = 'pseudomasked_genomic_flair_guided'
tool = 'flair'
config_tsv = f'snakemake/config_{analysis}.tsv'
df = parse_config(config_tsv)
df['analysis'] = analysis
input_gtf = config[analysis]['gtf']

# df = df.loc[df.tech_rep.isin(['GM10493_1',
#                               'GM12878_1',
#                               'GM22300_1'])]

# df = df.loc[df.tech_rep=='HG03732_1']

wildcard_constraints:
    tech_rep='|'.join([re.escape(x) for x in df.tech_rep.tolist()]),
    lab_rep='|'.join([re.escape(x) for x in df.lab_rep.tolist()]),
    end_mode='|'.join([re.escape(x) for x in end_modes]),


rule all:
    input:
        expand(config['cerberus']['update']['gtf'],
               tech_rep='HG03732_1',
               analysis=analysis)

        # expand(config['cerberus']['ends'],
        #        tech_rep=df.tech_rep.tolist(),
        #        end_mode=end_modes,
        #        analysis=analysis)

        # expand(config['cerberus']['merge']['h5'],
        #  analysis=analysis),
        # expand(config['cerberus']['merge']['gtf'],
        #        analysis=analysis)




        # expand(config['cerberus']['update']['gtf'],
        #        tech_rep=df.tech_rep.tolist(),
        #        analysis=analysis)

        # config['cerberus']['ref']['h5']
        # expand(config['cerberus']['ics']),
        # expand(config['cerberus']['ends'],
        #        tech_rep=df.tech_rep.tolist(),
        #        end_mode=end_modes)

def get_df_val(df, col, col2_val, col2='sample'):
    temp = df.loc[df[col2]==col2_val]
    assert len(temp.index) == 1
    return temp[col].values[0]

##########################
####### Novel gene stuff
##########################
# give novel genes the same IDs across the different samples
rule get_novel_gene_gtf:
    input:
        gtf = lambda wc: expand(input_gtf,
                                lab_rep=get_df_val(df,
                                'lab_rep',
                                wc.tech_rep,
                                'tech_rep'))[0]
    params:
        tool = tool
    resources:
        threads = 1,
        nodes = 1
    output:
        # gtf = temporary(config['fmt']['novel_gene_gtf'])
        gtf = config['fmt']['novel_gene_gtf']
    run:
        get_novel_gene_gtf(input.gtf, output.gtf, params.tool)

# give novel genes the same IDs across the different samples
rule get_novel_gene_bed:
    input:
        gtf = lambda wc: expand(input_gtf,
                                lab_rep=get_df_val(df,
                                'lab_rep',
                                wc.tech_rep,
                                'tech_rep'))[0]
    params:
        tool = tool
    resources:
        threads = 1,
        nodes = 1
    output:
        bed = temporary(config['fmt']['novel_gene_bed'])
    run:
        get_novel_gene_bed(input.gtf, output.bed, params.tool)

rule cat_novel_gene_gtfs:
    input:
        gtfs = lambda wc: expand(rules.get_novel_gene_gtf.output.gtf,
                                 analysis=wc.analysis,
                                 tech_rep=df.tech_rep.tolist())
    params:
        cli_gtf = lambda wc: fmt_list_for_cli(expand(rules.get_novel_gene_gtf.output.gtf,
                                 analysis=wc.analysis,
                                 tech_rep=df.tech_rep.tolist()))
    resources:
        threads = 1,
        nodes = 1
    output:
        gtf = temporary(config['fmt']['novel_gene_merge_gtf'])
    shell:
        """
        cat {params.cli_gtf} > {output.gtf}
        """

# build loci using Julien's tool
rule buildloci:
    input:
        gtf = rules.cat_novel_gene_gtfs.output.gtf
    resources:
        threads = 1,
        nodes = 1
    conda:
        'base'
    params:
        buildLoci = config['software']['buildLoci'],
        prefix = 'novel_gene'
    output:
        gtf = temporary(config['fmt']['novel_gene_merge_build_loci_gtf'])
    shell:
        """
        module load bedtools
        bedtools intersect \
            -s \
            -wao \
            -a {input.gtf} \
            -b {input.gtf} | \
            {params.buildLoci} - \
                locPrefix {params.prefix} \
                > {output.gtf}
        """

rule get_novel_tid_gid:
    input:
        gtf = rules.buildloci.output.gtf
    resources:
        threads = 1,
        nodes = 1
    output:
        tsv = config['fmt']['novel_gene_tid_to_gid']
    run:
        df = pr.read_gtf(input.gtf, as_df=True)
        df = df[['transcript_id', 'gene_id']].drop_duplicates()
        assert len(df.loc[df.transcript_id.duplicated(keep=False)].index) == 0
        df.to_csv(output.tsv, sep='\t', index=False)

rule fmt_novel_gene_rename:
    input:
        tsv = rules.get_novel_tid_gid.output.tsv,
        gtf = lambda wc: expand(input_gtf,
                                lab_rep=get_df_val(df,
                                'lab_rep',
                                wc.tech_rep,
                                'tech_rep'))[0]
    resources:
        threads = 1,
        nodes = 1
    params:
        tool = tool
    output:
        # gtf = temporary(config['fmt']['novel_gene_rename_gtf'])
        gtf = config['fmt']['novel_gene_rename_gtf']

    run:
        rename_novel_genes(input.gtf,
                           input.tsv,
                           output.gtf,
                           params.tool)

# # merge novel gene intervals across samples
# # and give them a unique number
# rule get_merged_novel_gene_bed:
#     input:
#         beds = expand(rules.get_novel_gene_bed.output.bed,
#                      analysis=analysis,
#                      tech_rep=df.tech_rep.tolist())
#     resources:
#         threads = 1,
#         nodes = 2
#     output:
#         bed = temporary(config['fmt']['novel_gene_merge_bed'])
#     run:
#         merge_beds(list(input.beds), output.bed)
#
# rule fmt_novel_gene_rename:
#     input:
#         bed = config['fmt']['novel_gene_merge_bed'],
#         gtf = lambda wc: expand(input_gtf,
#                                 lab_rep=get_df_val(df,
#                                 'lab_rep',
#                                 wc.tech_rep,
#                                 'tech_rep'))[0]
#     resources:
#         threads = 1,
#         nodes = 1
#     params:
#         tool = tool
#     output:
#         gtf = temporary(config['fmt']['novel_gene_rename_gtf'])
#     run:
#         rename_novel_genes(input.gtf,
#                            input.bed,
#                            output.gtf,
#                            params.tool)

# format the gtf corrrectly first
rule fmt_gtf:
    input:
        gtf = config['fmt']['novel_gene_rename_gtf'],
        ref_gtf = config['ref']['gtf']
    resources:
        threads = 1,
        nodes = 2
    params:
        tool = tool
    output:
        gtf = config['fmt']['gtf']
    conda:
        'cerberus'
    shell:
        """
        python snakemake/refmt_gtf.py \
            {input.gtf} \
            {input.ref_gtf} \
            {output.gtf} \
            {params.tool}
        """


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
        bed = temporary(config['ref']['cerberus']['ends'])

use rule gtf_to_ic as cerb_gtf_to_ic with:
    input:
        gtf = rules.fmt_gtf.output.gtf
        # gtf = lambda wc: expand(input_gtf,
        #                         lab_rep=get_df_val(df,
        #                                  'lab_rep',
        #                                  wc.tech_rep,
        #                                  'tech_rep'))
    output:
        tsv = temporary(config['cerberus']['ics'])

use rule agg_ics_cfg as cerb_agg_ics_cfg with:
    input:
        ref_tsv = rules.ref_gtf_to_ic.output.tsv,
        tsvs = expand(rules.cerb_gtf_to_ic.output.tsv,
                     analysis=analysis,
                     tech_rep=df.tech_rep.tolist())
    params:
        ref_source = config['ref']['gtf_ver'],
        sources = df.tech_rep.tolist()
    output:
        cfg = temporary(config['cerberus']['agg']['ics_cfg'])

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
        # bed = temporary(config['cerberus']['ends'])
        bed = config['cerberus']['ends']


use rule agg_ends_cfg as cerb_agg_config with:
    input:
        ref_ends = rules.ref_gtf_to_ends.output.bed,
        sample_ends = lambda wc: expand(rules.cerb_gtf_to_ends.output.bed,
                    analysis=wc.analysis,
                    tech_rep=df.tech_rep.tolist(),
                    end_mode=wc.end_mode)
    params:
        add_ends = [True]+[True for i in range(len(df.tech_rep.tolist()))],
        refs = [True]+[False for i in range(len(df.tech_rep.tolist()))],
        ref_source = config['ref']['gtf_ver'],
        sample_sources = df.tech_rep.tolist()
    output:
        cfg = temporary(config['cerberus']['agg']['ends_cfg'])

use rule agg_ends as cerb_agg_ends with:
    input:
        ref_ends = rules.cerb_agg_config.input.ref_ends,
        sample_ends = rules.cerb_agg_config.input.sample_ends,
        cfg = rules.cerb_agg_config.output.cfg
    params:
        slack = lambda wc: config['params']['cerberus'][wc.end_mode]['agg_slack']
    output:
        agg_ends = temporary(config['cerberus']['agg']['ends'])

use rule write_ref as cerb_write_ref with:
    input:
        tss = expand(rules.cerb_agg_ends.output.agg_ends,
                     analysis=analysis,
                     end_mode='tss'),
        tes = expand(rules.cerb_agg_ends.output.agg_ends,
                     analysis=analysis,
                     end_mode='tes'),
        ics = rules.cerb_agg_ic.output.tsv
    output:
        h5 = config['cerberus']['ref']['h5']

use rule annot_transcriptome as cerb_annot_sample with:
    input:
        gtf = rules.cerb_gtf_to_ic.input.gtf,
        h5 = rules.cerb_write_ref.output.h5
    params:
        source = lambda wc: wc.tech_rep
    output:
        h5 = config['cerberus']['ref']['annot_h5']

use rule update_gtf as cerb_update_gtf_sample with:
    input:
        gtf = rules.cerb_annot_sample.input.gtf,
        h5 = rules.cerb_annot_sample.output.h5
    params:
        source = lambda wc: wc.tech_rep
    output:
        gtf = config['cerberus']['update']['gtf']

use rule annot_transcriptome as cerb_annot_ref with:
    input:
        gtf = rules.ref_gtf_to_ic.input.gtf,
        h5 = rules.cerb_write_ref.output.h5
    params:
        source = config['ref']['gtf_ver']
    output:
        h5 = config['ref']['cerberus']['ref']['annot_h5']

use rule update_gtf as cerb_update_gtf_ref with:
    input:
        gtf = rules.ref_gtf_to_ic.input.gtf,
        h5 = rules.cerb_annot_ref.output.h5
    params:
        source = config['ref']['gtf_ver']
    output:
        gtf = config['ref']['cerberus']['update']['gtf']


################################
############ merge cerb gtfs
################################
rule make_cerb_agg_gtf_cfg:
    input:
        gtfs = lambda wc: expand(rules.cerb_update_gtf_sample.output.gtf,
                                 analysis=wc.analysis,
                                 tech_rep=df.tech_rep.tolist()),
    resources:
        threads = 1,
        nodes = 1
    output:
        tsv = temporary(config['cerberus']['merge']['cfg'])
    run:
        make_cerb_agg_gtf_cfg(input.gtfs, output.tsv)

rule cerb_merge_gtf:
    input:
        gtfs = rules.make_cerb_agg_gtf_cfg.input.gtfs,
        cfg = rules.make_cerb_agg_gtf_cfg.output.tsv,
        h5 = rules.cerb_write_ref.output.h5
    output:
        gtf = config['cerberus']['merge']['gtf']
    conda:
        'cerberus'
    resources:
        nodes = 4,
        threads = 1
    shell:
        """
        python snakemake/cerb_agg_gtfs.py \
            {input.cfg} \
            {input.h5} \
            {output.gtf}
        """


##
rule make_cerb_agg_t_map_cfg:
    input:
        h5s = lambda wc: expand(rules.cerb_annot_sample.output.h5,
                                 analysis=wc.analysis,
                                 tech_rep=df.tech_rep.tolist()),
    resources:
        threads = 1,
        nodes = 1
    output:
        tsv = temporary(config['cerberus']['merge']['cfg_h5'])
    run:
        make_cerb_agg_gtf_cfg(input.h5s, output.tsv)

rule cerb_merge_t_map:
    input:
        gtfs = rules.make_cerb_agg_t_map_cfg.input.h5s,
        cfg = rules.make_cerb_agg_t_map_cfg.output.tsv
    output:
        h5 = config['cerberus']['merge']['h5']
    conda:
        'cerberus'
    resources:
        nodes = 4,
        threads = 1
    shell:
        """
        python snakemake/cerb_agg_t_maps.py \
            {input.cfg} \
            {output.h5}
        """


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
