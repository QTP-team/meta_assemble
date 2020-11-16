### tianliu@genomics.cn
### 2020-08-20

import pandas as pd
import argparse

parser = argparse.ArgumentParser(
        usage = 'python MAGs_summary.py *.depth -o merge.depth',
        description = "merge jgi_summarize_bam_contig_depths results"
        )
parser.add_argument("infiles", nargs = "+", help = "One or more files to join")
parser.add_argument("-o", help = "outfile")

args = parser.parse_args()
infiles = args.infiles
outfile = args.o

def load_df(depth_file):
    tmp = pd.read_table(depth_file).drop(columns='totalAvgDepth')
    tmp.iloc[:,2] = tmp.iloc[:,2].round(decimals = 6)
    tmp.iloc[:,3] = tmp.iloc[:,3].round(decimals = 6)
    return tmp

merge_res = load_df(infiles[0])
for infile in infiles[1:]:
    tmp_df = load_df(infile)
    merge_res = merge_res.merge(tmp_df, on = ['contigName', 'contigLen'])

col_list = [x for x in list(merge_res) if x.endswith('sorted.bam')]
merge_res['totalAvgDepth'] = merge_res[col_list].sum(axis=1)

totalAvgDepth = merge_res.pop('totalAvgDepth')
merge_res.insert(2,'totalAvgDepth',totalAvgDepth)

merge_res.to_csv(outfile, index=None, sep='\t')

print("merge %s depth files to %s" %(len(args.infiles), outfile))
