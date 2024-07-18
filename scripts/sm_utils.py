import pandas as pd
import pyranges as pr
import os
import yaml

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
    # TODO possibly make a get novel genes function
    # should ask for what level ie transcript gene exon
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

def load_config(config_file=None):
    """
    Load snakemake config file as a dictionary
    """
    if not config_file:
        d = os.path.dirname(__file__)
        od = f'{d}/../snakemake/'
        config_file = f'{od}/config.yml'

    with open(config_file) as f:
        config = yaml.safe_load(f)

    return config

def merge_beds(beds, bed):
    """
    Merge beds into non-overlapping invervals
    """
    df = pd.DataFrame()
    for b in beds:
        temp = pr.read_bed(b).df
        df = pd.concat([df, temp], axis=0)
    df = pr.PyRanges(df)
    df = df.merge(strand=True,
              slack=0)

    # assign each a random number
    df = df.df
    df['num'] = [i for i in range(len(df.index))]
    df['Name'] = 'novel_gene_'+df['num'].astype(str)
    df.drop('num', axis=1, inplace=True)

    df = pr.PyRanges(df)
    df.to_bed(bed)
