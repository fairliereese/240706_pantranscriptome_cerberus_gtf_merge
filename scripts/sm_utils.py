import pandas as pd
import pyranges as pr

def parse_config(fname):
    df = pd.read_csv(fname, sep='\t')
    df['sample'] = df.lab_rep.str.rsplit('_', n=1, expand=True)[1]

    # get tech rep id -- each lab replicate w/ the same sample is a tech rep
    temp = df.drop_duplicates()
    temp['tech_rep_num'] = temp.sort_values(['sample',
                             'lab_rep'],
                              ascending=[True, True])\
                              .groupby(['sample']) \
                              .cumcount() + 1
    temp.sort_values(by=['sample', 'tech_rep_num'])
    df = df.merge(temp, how='left', on=['sample', 'lab_rep'])

    df['tech_rep'] = df['sample']+'_'+df['tech_rep_num'].astype(str)

    return df

def make_cerb_agg_gtf_cfg(gtfs, ofile):
    df = pd.DataFrame()
    df['gtf'] = gtfs
    df.to_csv(ofile, index=False, header=None, sep='\t')

def get_novel_gene_bed(gtf, bed, how='iq'):
    df = pr.read_gtf(gtf).df

    # for isoquant
    if how == 'iq':
        df = df.loc[(df.Feature=='gene')&\
                    (df.gene_id.str.contains('novel_gene'))]

    # get just bed fields
    df = df[['Chromosome', 'Source', 'Start', 'End', 'Strand', 'gene_id']]
    df.rename({'gene_id': 'Name'}, axis=1, inplace=True)

    # merge overlapping guys
    df = pr.PyRanges(df)
    df = df.merge(strand=True,
                  slack=0)

    # save to bed
    df.to_bed(bed)
