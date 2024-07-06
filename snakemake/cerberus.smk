end_modes = ['tss', 'tes']

# reusable rules
rule agg_ics_cfg:
    resources:
        threads = 1,
        nodes = 1
    run:
        df = pd.DataFrame()
        df['fname'] = [input.ref_tsv]+list(input.tsvs)
        df['ref'] = [True]+[False for i in range(len(list(input.tsvs)))]
        df['source'] = [params.ref_source]+params.sources
        df.to_csv(output.cfg, header=None, index=None)

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
    run:
        gtf_to_bed(input.gtf,
                   wildcards.end_mode,
                   output.bed,
                   dist=params.extend,
                   slack=params.slack)
