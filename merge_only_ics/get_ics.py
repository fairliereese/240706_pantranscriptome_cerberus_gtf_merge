import argparse
import cerberus
import pyranges as pr

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Process a GTF file to get unique intron chains.")
    parser.add_argument("input_gtf", type=str, help="Path to the input GTF file.")
    parser.add_argument("output_ics", type=str, help="Path to the output CSV file for intron chains.")

    # Parse the arguments
    args = parser.parse_args()

    # Read the GTF file and process it
    df = pr.read_gtf(args.input_gtf, rename_attr=True, duplicate_attr=True)
    df = df.df
    df = pr.PyRanges(df)
    df = cerberus.get_ic(df)

    # Write the resulting intron chains to the output CSV file
    df.to_csv(args.output_ics, sep='\t')

if __name__ == "__main__":
    main()
