import sys
import pyranges as pr
import pandas as pd
import cerberus

def merge_cerb_h5s(h5s, ofile):

    df = pd.DataFrame()
    for f in h5s:
        ca = cerberus.read(f)
        t_map = ca.t_map.copy(deep=True)
        df = pd.concat([df, t_map], axis=0)

    ca.t_map = t_map
    ca.write(ofile)


def main():
    if len(sys.argv) != 3:
        print("Usage: python cerb_agg_t_maps.py <input txt file w/ list of h5s> <output_h5_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    df = pd.read_csv(input_file, sep='\t', header=None)
    gtfs = df[0].tolist()
    ofile = sys.argv[2]

    merge_cerb_h5s(h5s, ofile)

if __name__ == "__main__":
    main()
