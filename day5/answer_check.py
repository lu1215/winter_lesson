import pandas as pd
df1 = pd.read_csv("miRNA_domain_map_id.csv")
df2 = pd.read_csv("miRNA_domain_map_id_v2.csv")

total_df = pd.concat([df1, df2])
# print(total_df)
differ_df = total_df.drop_duplicates(keep=False)
print(differ_df)