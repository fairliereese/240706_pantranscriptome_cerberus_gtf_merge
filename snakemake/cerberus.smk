end_modes = ['tss', 'tes']

# reusable rules
rule agg_ics_cfg:
    resources:
        threads = 1,
        nodes = 1
    run:
        temp = pd.DataFrame()
        temp['fname'] = list([input.ref_tsv])+list(input.tsvs)
        temp['ref'] = [True]+[False for i in range(len(list(input.tsvs)))]
        temp['source'] = [params.ref_source]+params.sources
        temp.to_csv(output.cfg, header=None, index=None)

rule gtf_to_ic:
    resources:
        threads = 1,
        nodes = 2
    conda:
        'cerberus'
    shell:
        """
        cerberus gtf_to_ics \
            --gtf {input.gtf} \
            -o {output.tsv}
        """

rule agg_ics:
    resources:
        threads = 1,
        nodes = 1
    conda:
        'cerberus'
    shell:
        """
        cerberus agg_ics \
            --input {input.cfg} \
            -o {output.tsv}
        """

rule gtf_to_ends:
    resources:
        threads = 1,
        nodes = 1
    conda:
        'cerberus'
    shell:
        """
        cerberus gtf_to_bed \
            --gtf {input.gtf} \
            --mode {wildcards.end_mode} \
            --dist {params.dist} \
            --slack {params.slack} \
            -o {output.bed}
        """

rule agg_ends_cfg:
    resources:
        threads = 1,
        nodes = 1
    run:
        temp = pd.DataFrame()
        temp['fname'] = list([input.ref_ends])+list(input.sample_ends)
        temp['ref'] = params.refs
        temp['add_ends'] = params.add_ends
        temp['source'] = [params.ref_source]+params.sample_sources
        temp = temp[['fname', 'add_ends', 'ref', 'source']]
        temp.to_csv(output.cfg, header=None, index=None)

rule agg_ends:
    resources:
        threads = 1,
        nodes = 2
    conda:
        'cerberus'
    shell:
        """
        cerberus agg_ends \
            --input {input.cfg} \
            --mode {wildcards.end_mode} \
            --slack {params.slack} \
            -o {output.agg_ends}
        """

rule write_ref:
    resources:
        nodes = 1,
        threads = 1
    conda:
        'cerberus'
    shell:
        """
        cerberus write_reference \
            --tss {input.tss} \
            --tes {input.tes} \
            --ics {input.ics} \
            -o {output.h5}
        """

rule annot_transcriptome:
    resources:
        threads = 1,
        nodes = 2
    conda:
        'cerberus'
    shell:
        """
        cerberus annotate_transcriptome \
            --gtf {input.gtf} \
            --h5 {input.h5} \
            --source {params.source} \
            -o {output.h5}
        """

rule update_gtf:
    resources:
        threads = 1,
        nodes = 2
    conda:
        'cerberus'
    shell:
        """
        cerberus replace_gtf_ids \
            --h5 {input.h5} \
            --gtf {input.gtf} \
            --source {params.source} \
            --update_ends \
            --collapse \
            -o {output.gtf}
        """
