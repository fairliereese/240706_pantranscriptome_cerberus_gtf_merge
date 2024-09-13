import cerberus
import pyranges as pr
import pandas as pd
import argparse

def df_to_gtf(df):
    gtf_entries = []

    for index, row in df.iterrows():
        chromosome = row['Chromosome']
        strand = row['Strand']
        ic_coords = list(map(int, row['ic'].split('-')))
        tss = row['tss']
        tes = row['tes']
        source = row['source']
        samples = source.split(',')

        # change tss and tes depending on strand to combat
        # off by one errors
        if strand == '+':
            tes += 1
        elif strand == '-':
            tss += 1

        # Ensure Start < End
        if tss > tes:
            tss, tes = tes, tss

        # Create transcript entry
        transcript_id = f"transcript_{index}"
        gtf_entries.append({
            'Chromosome': chromosome,
            'Source': 'ChatGPT',
            'Feature': 'transcript',
            'Start': tss,
            'End': tes,
            'Score': '.',
            'Strand': strand,
            'Frame': '.',
            'transcript_id': transcript_id,
            'samples': ','.join(samples)
        })

        # Generate exon entries
        exons = []
        if strand == '+':
            # For forward strand genes
            coords = [tss] + ic_coords + [tes]
            exons = [(coords[i], coords[i + 1]) for i in range(0,len(coords) - 1,2)]
        else:
            # For reverse strand genes
            coords = [tss] + ic_coords[::-1] + [tes]
            exons = [[coords[i], coords[i + 1]] for i in range(0,len(coords) - 1,2)][::-1]

        # Add exon entries
        for i, (start, end) in enumerate(exons):
            if start > end:
                start, end = end, start
            gtf_entries.append({
                'Chromosome': chromosome,
                'Source': 'ChatGPT',
                'Feature': 'exon',
                'Start': start,
                'End': end,
                'Score': '.',
                'Strand': strand,
                'Frame': '.',
                'transcript_id': transcript_id,
                'exon_number': i + 1,
                'samples': ','.join(samples)
            })

    # Convert the list of entries to a DataFrame
    gtf_df = pd.DataFrame(gtf_entries, columns=[
        'Chromosome', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand', 'Frame', 'transcript_id', 'exon_number', 'samples'
    ])

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
        # remove dupes
        source_ic_df.drop_duplicates(inplace=True, keep='first')

        source_ic_df = source_ic_df.groupby(gb_cols, observed=True).agg({'source': ','.join}).reset_index()

        # if there are no transcripts, have to get tss / tes differently
        # looking at you lyric
        if 'transcript' not in gtf_df.df.Feature.unique().tolist():
            fwd, rev = cerberus.get_stranded_gtf_dfs(gtf_df.df)
            tss_df = pd.DataFrame()
            tes_df = pd.DataFrame()
            for strand, strand_df in zip(['+', '-'], [fwd,rev]):
                strand_df['max_coord'] = strand_df[['Start', 'End']].max(axis=1)
                strand_df['min_coord'] = strand_df[['Start', 'End']].min(axis=1)
                strand_df = strand_df[['transcript_id',
                                       'min_coord',
                                       'max_coord']].groupby('transcript_id').agg(
                                           {'min_coord': min,
                                            'max_coord': max}).reset_index()

                if strand == '-':
                    strand_df['tss_end'] = strand_df['max_coord']
                    strand_df['tss_start'] = strand_df['max_coord']-1
                    strand_df['tes_start'] = strand_df['min_coord']
                    strand_df['tes_end'] = strand_df['tes_start']+1
                elif strand == '+':
                    strand_df['tss_start'] = strand_df['min_coord']
                    strand_df['tss_end'] = strand_df['min_coord'] + 1
                    strand_df['tes_start'] = strand_df['max_coord'] - 1
                    strand_df['tes_end'] = strand_df['max_coord']
                tss_df = pd.concat([tss_df,
                            strand_df[['transcript_id', 'tss_start', 'tss_end']]],
                            axis=0)
                tes_df = pd.concat([tes_df,
                            strand_df[['transcript_id', 'tes_start', 'tes_end']]],
                            axis=0)
            tss_df.rename({'tss_start': 'Start',
                           'tss_end': 'End'},
                          axis=1, inplace=True)
            tes_df.rename({'tes_start': 'Start',
                           'tes_end': 'End'},
                          axis=1, inplace=True)
        else:
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
    df = df_to_gtf(df)

    df = pr.PyRanges(df)
    df.to_gtf(args.output_ics)

if __name__ == "__main__":
    main()
