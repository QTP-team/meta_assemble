### tianliu@genomics.cn
### 2020-05-11

import sys
import pandas as pd

trim_file = sys.argv[1]
rmhost_file = sys.argv[2]
merge_file = sys.argv[3]

trim = pd.read_table(trim_file)
rmhost = pd.read_table(rmhost_file)

df = trim.merge(rmhost, on = ["Library_ID"])
df['Nohomo_reads_ratio'] = df['Nohomo_reads_count'] / df['Clean_reads_count']
df['Nohomo_bases_ratio'] = df['Nohomo_bases_count'] / df['Clean_bases_count']
df['Nohomo_reads_ratio'] = df['Nohomo_reads_ratio'].apply(lambda x: format(x, '.2%'))
df['Nohomo_bases_ratio'] = df['Nohomo_bases_ratio'].apply(lambda x: format(x, '.2%'))

df.to_csv(merge_file, sep="\t", index = False)