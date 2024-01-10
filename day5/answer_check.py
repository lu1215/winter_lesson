import pandas as pd
df1 = pd.read_csv("miRNA_domain_map_id.csv")
df2 = pd.read_csv("gene_domain_map_id.csv")

total_df = pd.concat([df1, df2])
print(total_df)
differ_df = total_df.drop_duplicates(keep=False)
print(differ_df)