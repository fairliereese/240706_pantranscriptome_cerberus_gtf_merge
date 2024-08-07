import sys
import pyranges as pr
import cerberus
import pandas as pd

def make_hier_entry(df, how='t'):
    """
    kind {'g','t'}
    """
    agg_dict = {'min_coord': 'min', 'max_coord': 'max'}
    t_df = df.copy(deep=True)
    t_df['min_coord'] = t_df[['Start', 'End']].min(axis=1)
    t_df['max_coord'] = t_df[['Start', 'End']].max(axis=1)
    if how == 't':
        gb_cols = ['Chromosome', 'Strand', 'gene_name',
                   'gene_id', 'transcript_id', 'transcript_name',
                   'tss_id', 'tes_id',
                   'new_transcript_id', 'original_transcript_id',
                   'original_transcript_name', 'ag1', 'ag2']
    elif how == 'g':
        gb_cols = ['Chromosome', 'Strand', 'gene_name',
                   'gene_id']
    gb_cols = list(set(gb_cols)&(set(t_df.columns)))

    cols = gb_cols + ['min_coord', 'max_coord']
    t_df = t_df[cols]
    t_df = t_df.groupby(gb_cols, observed=True).agg(agg_dict).reset_index()
    t_df.rename({'min_coord': 'Start', 'max_coord': 'End'}, axis=1, inplace=True)
    if how == 't':
        t_df['Feature'] = 'transcript'
    elif how == 'g':
        t_df['Feature'] = 'gene'

    return t_df

def fmt_gtf(ifile, ref_file, ofile, tool):
    """
    Add gene name and transcript name
    columns in a more sensible way.
    """

    df = pr.read_gtf(ifile).df

    # limit to just transcripts and exons b/c we'll just
    # reconstruct the genes
    df = df.loc[df.Feature.isin(['transcript', 'exon'])]
    df['gene_name'] = df['gene_id']

    gtf_df = df.copy(deep=True)

    # flair -- rename antisense genes
    if tool == 'flair':
        ref_df = pr.read_gtf(ref_file).df

        # get strand, gene id, and transcript id of flair transcripts
        df = df[['gene_id', 'transcript_id', 'Strand']].drop_duplicates()
        df['gid_stable'] = cerberus.get_stable_gid(df, 'gene_id')

        #  error check to prevent buildLoci errors
        temp = df[['gene_id', 'Strand']].drop_duplicates()
        temp = temp.loc[temp.gene_id.duplicated(keep=False)]
        temp = temp.loc[temp.gene_id.str.startswith('LOC')]
        if len(temp.index) != 0:
            raise ValueError('Found buildLoci locations w/ + and - strand')

        # get strand and gene id for reference
        ref_df = ref_df[['gene_id', 'Strand']].drop_duplicates()
        ref_df['gid_stable'] = cerberus.get_stable_gid(ref_df, 'gene_id')

        # merge w/ reference to figure out which transcripts are antisense
        df = df.merge(ref_df,
                      how='left',
                      on='gid_stable',
                      suffixes=('', '_ref'))

        # limit to things on the opposite strand as ref
        # and known things
        # and replace gene id for those
        # import pdb; pdb.set_trace()
        temp = df.loc[(df.Strand!=df.Strand_ref)&\
                      (df.Strand_ref.notnull())]
        tids = temp.transcript_id.unique().tolist()
        inds = gtf_df.loc[gtf_df.transcript_id.isin(tids)].index

        gtf_df.loc[inds, 'gene_id'] = gtf_df.loc[inds, 'gene_id']+'_antisense'

        l1 = len(gtf_df.gene_id.unique())
        l2 = len(gtf_df[['Strand', 'gene_id']].drop_duplicates())

        # make sure we have the same number of gene ids at the end (l2)
        # as we did unique gene id + strand combos at the beginning (l1)
        assert l1 == l2

        # definitionally we are making novel loci by merging
        # on strands. So make sure we don't have any novel genes
        # w/ antisense designation
        assert len(df.loc[(df.gene_id.str.contains('LOC'))&\
                     (df.gene_id.str.contains('antisense'))].index) == 0

    # make new gene entries for everything
    df = gtf_df.copy(deep=True)
    l1 = len(df.gene_id.unique().tolist())

    # make gene entry
    g_df = make_hier_entry(df, how='g')

    if tool == 'espresso': source = 'Espresso'
    elif tool == 'flair': source = 'FLAIR'
    elif tool == 'iq': source = 'IsoQuant'
    g_df['Source'] = source
    g_df['Frame'] = '.'
    g_df['Score'] = '.'
    l2 = len(g_df.loc[g_df.Feature=='gene'].index)
    assert l1 == l2

    # concat them and then sort gtf
    df = pd.concat([df, g_df], axis=0)
    df = cerberus.sort_gtf(df)

    # convert + save
    df = pr.PyRanges(df)
    df.to_gtf(ofile)

def main():
    if len(sys.argv) != 5:
        print("Usage: python script_name.py <input_gtf_file> <ref_gtf_file> <output_gtf_file> <tool>")
        sys.exit(1)

    input_file = sys.argv[1]
    ref_file = sys.argv[2]
    output_file = sys.argv[3]
    tool = sys.argv[4]

    fmt_gtf(input_file, ref_file, output_file, tool)

if __name__ == "__main__":
    main()
