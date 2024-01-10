import pandas as pd

def create_domain_table():
    df_ori = pd.read_csv("../day4/Homo_sapiens_miRNA.csv")[['mirna_name', 'gene_name']]
    # df_ori = pd.read_csv("Homo_sapiens_TarBase-v9.tsv", sep='\t')[['mirna_name', 'gene_name']]
    df_ori.drop_duplicates(inplace=True)
    df_ori.dropna(inplace=True)
    # print(len(df_ori[["gene_name"]].drop_duplicates()))
    df_ori['gene_name'] = df_ori['gene_name'].astype(str)
    df = df_ori.groupby('mirna_name').agg({'gene_name': ','.join})
    df["count"] = df_ori.groupby(by=['mirna_name']).agg({'gene_name': len})['gene_name']
    df = df.sort_values(by="mirna_name").reset_index()
    print(max(df["count"]))
    new_col = ["mirna_name", "count",  "gene_name"]
    df = df.reindex(columns=new_col)
    print(f"最多基因的:\n {df.iloc[df['count'].idxmax()]}")
    print(f"最少基因的:\n {df.iloc[df['count'].idxmin()]}")
    df.to_csv("miRNA_domain_map_id_v2.csv", index=False)

# def count_number():
#     df = pd.read_csv("protein_domain_map_id.csv")
#     print(len(df[["mirna_name"]].drop_duplicates()))
#     print(df["count"].sum())

if __name__ == "__main__":
    create_domain_table()
    # count_number()

