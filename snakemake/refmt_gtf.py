import sys
import pyranges as pr

def fmt_iq_gtf(ifile, ofile):
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

    # convert + save
    df = pr.PyRanges(df)
    df.to_gtf(ofile)

def main():
    if len(sys.argv) != 3:
        print("Usage: python script_name.py <input_gtf_file> <output_gtf_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    fmt_iq_gtf(input_file, output_file)

if __name__ == "__main__":
    main()
