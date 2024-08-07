import cerberus
import pyranges as pr
import pandas as pd
import argparse

def make_ic_gtf(gtf_files, ofile)
    ic_df = pd.DataFrame()
    for f in gtf_files:
        gtf_df = pr.read_gtf(f, duplicate_attr=True)
        gtf_df = gtf_df.df
        gtf_df = pr.PyRanges(gtf_df)
        df = cerberus.get_ic(gtf_df)

        # remove monoxonic
        df = df.loc[df.ic!='-']

        # merge to get starts for each sample-level thing
        tss_df = gtf_df.features.tss().df
        tes_df = gtf_df.features.tes().df

        tss_df = tss_df[['transcript_id', 'Start']].rename({'Start':'tss_start'}, axis=1)
        tes_df = tes_df[['transcript_id', 'Start']].rename({'Start':'tes_start'}, axis=1)

        df = df[['Chromosome', 'Strand',
                 'transcript_id', 'ic']]
        df = df.merge(tss_df, how='left', on='transcript_id')
        df = df.merge(tes_df, how='left', on='transcript_id')
        df = df.drop('transcript_id', axis=1)

        # concat w/ original ic df
        ic_df = pd.concat([df, ic_df], axis=0)
        ic_df.drop_duplicates(inplace=True)

        # keep the longest for each
        fwd, rev = cerberus.get_stranded_gtf_dfs(ic_df)

        fwd = fwd.groupby(['Chromosome', 'Strand', 'ic']).agg(tss_start=("tss_start", "min"),
                                                              tes_start=("tes_start", "max")).reset_index()

        rev = rev.groupby(['Chromosome', 'Strand', 'ic']).agg(tss_start=("tss_start", "max"),
                                                              tes_start=("tes_start", "min")).reset_index()
        ic_df = pd.concat([fwd, rev], axis=0)

    ic_df.to_csv(ofile, sep='\t', index=False)

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Process a GTF file to get unique intron chains.")
    parser.add_argument("input_gtfs", type=str, help="CFG file with paths to the input GTF file.")
    parser.add_argument("output_ics", type=str, help="Path to the output CSV file for intron chains.")

    # Parse the arguments
    args = parser.parse_args()
    df = pd.read_csv(args.input_gtfs)
    gtfs = df['gtf'].tolist()


    # Process the GTF file
    make_ic_gtf(gtfs, args.output_ics)

if __name__ == "__main__":
    main()
