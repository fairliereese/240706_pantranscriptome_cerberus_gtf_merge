import sys
import pyranges as pr
import pandas as pd
import cerberus

def merge_cerb_gtfs(gtfs, h5, ofile):

    # get tss / tes stuff
    ca = cerberus.read(h5)
    tss = ca.tss
    tss.rename({'Name': 'tss_id'}, axis=1, inplace=True)
    tes = ca.tes
    tes.rename({'Name': 'tes_id'}, axis=1, inplace=True)
    tss = pr.PyRanges(tss)
    tes = pr.PyRanges(tes)

    # for each gtf, concat, sort, update ends, and agg
    df = pd.DataFrame()
    for f in gtfs:
        df = pd.concat([pr.read_gtf(f, as_df=True), df], axis=0)
        df = cerberus.sort_gtf(df)
        df = cerberus.update_gtf_ends(df, tss, tes)
        df = cerberus.agg_gtf(df)

    df = pr.PyRanges(df)
    df.to_gtf(ofile)


def main():
    if len(sys.argv) != 3:
        print("Usage: python cerb_agg_gtfs.py <input txt file w/ list of gtfs> <output_gtf_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    df = pd.read_csv(input_file, sep='\t', header=None)
    gtfs = df[0].tolist()
    h5 = sys.argv[2]
    ofile = sys.argv[3]

    merge_cerb_gtfs(gtfs, h5, ofile)

if __name__ == "__main__":
    main()
