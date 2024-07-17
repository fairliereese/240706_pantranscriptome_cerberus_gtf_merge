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
