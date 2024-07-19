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

def get_gtf_novel_genes(gtf, tool):
    """
    Get a list of gene IDs from novel genes
    from a GTF
    """
    if tool == 'espresso':
        nov_gids = ['NA']
    elif tool == 'iq':
        nov_gids = gtf.loc[(gtf.Feature=='gene')&\
                         (gtf.gene_id.str.contains('novel_gene')&\
                         (gtf.Source=='IsoQuant'))].gene_id.unique().tolist()
    return nov_gids

def get_novel_gene_bed(gtf, bed):
    df = pr.read_gtf(gtf).df

    # get all entries belonging to novel genes
    novel_gids = get_gtf_novel_genes(df, tool)
    df = df.loc[df.gene_id.isin(novel_gids)]

    # get just bed fields
    df = df[['Chromosome', 'Source', 'Start', 'End', 'Strand', 'gene_id']]
    df.rename({'gene_id': 'Name'}, axis=1, inplace=True)

    # merge overlapping guys
    df = pr.PyRanges(df)
    df = df.merge(strand=True,
                  slack=0)

    # save to bed
    df.to_bed(bed)

def rename_novel_genes(ifile, bed, ofile, tool):
    """
    Rename novel genes based on the name of overlapping bed regions
    in a GTF file
    """

    # get entries w/ novel GIDs
    gtf_df = pr.read_gtf(ifile, as_df=True)
    nov_gids = get_gtf_novel_genes(gtf_df, tool)
    nov_df = gtf_df.loc[gtf_df.gene_id.isin(nov_gids)]
    l1 = len(nov_df.index)

    # merge w/ bed
    df = pr.read_bed(bed)
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

    gtf_df = pr.PyRanges(gtf_df)
    gtf_df.to_gtf(ofile)

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
