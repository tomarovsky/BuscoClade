#!/usr/bin/env python3
__author__ = 'tomarovsky'
import pandas as pd
import argparse
import re


def main():
    data_list = []
    for file_name in args.input:
        print(file_name)
        with open(file_name, 'r') as file:
            lines = file.readlines()
            target_line = lines[8]
            matches = re.findall(r'\d+\.\d+|\d+', target_line)
            numbers = [float(match) if '.' in match else int(match) for match in matches]
            species = file_name.split('/')[-1][14:-4].replace('_', ' ')
            line = [species, numbers[1], numbers[2], numbers[3], numbers[4], numbers[5]]
            data_list.append(line)
    df = pd.DataFrame(data_list, columns=['Species', 'S', 'D', 'F', 'M', 'N'])
    print(df)
    print(data_list)
    df.to_csv(args.output, sep='\t', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="BUSCO summaries to TSV file")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', nargs="+", type=str, help="short_summary_{species}.txt files")
    group_required.add_argument('-o', '--output', type=str, help="output TSV file name")
    args = parser.parse_args()
    main()