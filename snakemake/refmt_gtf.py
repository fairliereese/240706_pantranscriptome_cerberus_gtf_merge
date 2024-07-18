import sys
import pyranges as pr

def rename_novel_genes(gtf_df, bed, tool='iq'):
    """
    Rename novel genes using a unified set of novel gene
    bed coordinates across the different samples
    """
    if tool == 'iq':
        nov_df = gtf_df.loc[(gtf_df.Feature=='gene')&\
                    (gtf_df.gene_id.str.contains('novel_gene'))]
    l1 = len(nov_df.index)

    # merge w/ bed
    cols = nov_df.columns.tolist()
    nov_df = pr.PyRanges(nov_df)
    nov_df = nov_df.join(df,
                          how='left',
                          strandedness='same')
    nov_df = nov_df.df
    l2 = len(nov_df.index)
    assert l1 == l2
    assert len(nov_df.loc[nov_df.Name.isnull()].index) == 0

    # get a dictionary to map sample-level novel gene
    # ids to cross-sample novel gene ids
    gid_dict = dict([(o,n) for o,n in zip(nov_df.gene_id.tolist(),
                                          nov_df.Name.tolist())])

    # now update full gtf
    novel_gids = nov_df.gene_id.unique().tolist()
    inds = gtf_df.loc[gtf_df.gene_id.isin(novel_gids)].index
    gtf_df.loc[inds, 'gene_id'] = gtf_df.loc[inds, 'gene_id'].map(gid_dict)

    return gtf_df

def fmt_iq_gtf(ifile, bed, ofile, tool):
    """
    Add gene name and transcript name
    columns in a more sensible way.
    """
    df = pr.read_gtf(ifile).as_df()

    # just grab all the gene ids to be the gene names
    df['gene_name'] = df['gene_id']

    # for transcripts / other transcript-level features do the same but restrict to those feats
    inds = df.loc[df.Feature != 'gene'].index
    df.loc[inds, 'transcript_name'] = df.loc[inds, 'transcript_id']

    # fix up novel genes by orienting them using the same naming system
    df = rename_novel_genes(df, bed, tool)

    # convert + save
    df = pr.PyRanges(df)
    df.to_gtf(ofile)

def main():
    if len(sys.argv) != 4:
        print("Usage: python script_name.py <input_gtf_file> <input_bed_file> <output_gtf_file> <tool>")
        sys.exit(1)

    input_file = sys.argv[1]
    bed_file = sys.argv[2]
    output_file = sys.argv[3]
    tool = sys.argv[4]

    fmt_iq_gtf(input_file, bed_file, output_file, tool,)

if __name__ == "__main__":
    main()
