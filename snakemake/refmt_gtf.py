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

def fmt_gtf(ifile, ofile, tool):
    """
    Add gene name and transcript name
    columns in a more sensible way.
    """

    df = pr.read_gtf(ifile).df

    # espresso -- add gene entries
    if tool == 'espresso' or tool == 'flair':
        l1 = len(df.gene_id.unique().tolist())

        # make gene entry
        g_df = make_hier_entry(df, how='g')
        g_df['Source'] = 'Espresso'
        g_df['Frame'] = '.'
        g_df['Score'] = '.'
        l2 = len(g_df.loc[g_df.Feature=='gene'].index)
        assert l1 == l2

        # concat them and then sort gtf
        df = pd.concat([df, g_df], axis=0)
        df = cerberus.sort_gtf(df)

    # iq -- add gene names and transcript names
    elif tool == 'iq':

        # just grab all the gene ids to be the gene names
        df['gene_name'] = df['gene_id']

        # for transcripts / other transcript-level features do the
        # same but restrict to those feats
        inds = df.loc[df.Feature != 'gene'].index
        df.loc[inds, 'transcript_name'] = df.loc[inds, 'transcript_id']

    # convert + save
    df = pr.PyRanges(df)
    df.to_gtf(ofile)

def main():
    if len(sys.argv) != 4:
        print("Usage: python script_name.py <input_gtf_file> <output_gtf_file> <tool>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    tool = sys.argv[3]

    fmt_gtf(input_file, output_file, tool)

if __name__ == "__main__":
    main()
