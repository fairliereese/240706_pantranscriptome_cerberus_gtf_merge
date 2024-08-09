import cerberus
import pyranges as pr
import pandas as pd
import argparse

def process_table_to_gtf_df(df):
    gtf_entries = []

    for index, row in df.iterrows():
        chrom = row['Chromosome']
        strand = row['Strand']
        ic_list = list(map(int, row['ic'].split('-')))
        tss = int(row['tss'])
        tes = int(row['tes'])
        source = row['source']
        gene_id = f"gene_{index}"
        transcript_id = f"transcript_{index}"

        # Define exons based on strand
        if strand == '+':
            exon_coords = [(tss, ic_list[0])]  # First exon
            for i in range(1, len(ic_list) - 1, 2):
                exon_coords.append((ic_list[i], ic_list[i + 1]))  # Middle exons
            exon_coords.append((ic_list[-1], tes))  # Last exon
        else:  # Reverse strand '-'
            exon_coords = [(tes, ic_list[-1])]  # First exon (on reverse strand)
            for i in range(len(ic_list) - 2, 0, -2):
                exon_coords.append((ic_list[i + 1], ic_list[i]))  # Middle exons
            exon_coords.append((ic_list[0], tss))  # Last exon (on reverse strand)

        # Ensure each exon has start < end and create GTF entries
        for i, (start, end) in enumerate(exon_coords):
            start, end = min(start, end), max(start, end)
            attributes = f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; exon_number "{i+1}";'
            gtf_entries.append([chrom, source, 'exon', start, end, '.', strand, '.', attributes])

    # Convert the list of GTF entries to a DataFrame
    gtf_df = pd.DataFrame(gtf_entries, columns=['Chromosome', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Attributes'])

    return gtf_df


def make_ic(gtf_files):

    gb_cols = ['Chromosome', 'Strand', 'ic']

    ic_df = pd.DataFrame()
    source_ic_df = pd.DataFrame()
    for f in gtf_files:

        # get info about which each ic was detected in
        analysis = f.split('data/')[1].split('/')[0]

        tech_rep = f.rsplit('/', maxsplit=1)[1].split('.')[0]
        source = f'{analysis}_{tech_rep}'

        gtf_df = pr.read_gtf(f, duplicate_attr=True)
        gtf_df = gtf_df.df
        gtf_df = pr.PyRanges(gtf_df)
        df = cerberus.get_ic(gtf_df)
        df['source'] = source

        # remove monoxonic
        df = df.loc[df.ic!='-']

        # agg. sources; groupby and add commas
        source_ic_df = pd.concat([source_ic_df, df[gb_cols+['source']]],
                                 axis=0)
        source_ic_df = source_ic_df.groupby(gb_cols, observed=True).agg({'source': ','.join}).reset_index()

        # merge to get starts for each sample-level thing
        tss_df = gtf_df.features.tss().df
        tes_df = gtf_df.features.tes().df

        tss_df = tss_df[['transcript_id', 'Start']].rename({'Start':'tss'}, axis=1)
        tes_df = tes_df[['transcript_id', 'Start']].rename({'Start':'tes'}, axis=1)

        df = df[gb_cols+['transcript_id']]
        df = df.merge(tss_df, how='left', on='transcript_id')
        df = df.merge(tes_df, how='left', on='transcript_id')
        df = df.drop('transcript_id', axis=1)

        # concat w/ original ic df
        ic_df = pd.concat([df, ic_df], axis=0)
        ic_df.drop_duplicates(inplace=True)

        # keep the longest for each
        fwd, rev = cerberus.get_stranded_gtf_dfs(ic_df)

        fwd = fwd.groupby(gb_cols, observed=True).agg(tss=("tss", "min"),
                                                      tes=("tes", "max")).reset_index()

        rev = rev.groupby(gb_cols, observed=True).agg(tss=("tss", "max"),
                                                      tes=("tes", "min")).reset_index()
        ic_df = pd.concat([fwd, rev], axis=0)

    # merge in sources
    ic_df = ic_df.merge(source_ic_df, on=gb_cols, how='left')
    # ic_df.to_csv(ofile, sep='\t', index=False)
    return ic_df

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
    df = make_ic(gtfs)
    df = process_table_to_gtf_df(df)

    df = pr.PyRanges(df)
    df.to_gtf(args.output_ics)

if __name__ == "__main__":
    main()
