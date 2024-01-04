import pandas as pd

def create_domain_table():
    df_ori = pd.read_csv("../day4/Homo_sapiens_miRNA.csv")[['mirna_name', 'gene_name']]
    df_ori.drop_duplicates(inplace=True)
    print(len(df_ori[["mirna_name"]].drop_duplicates()))
    df = df_ori.groupby('gene_name').agg({'mirna_name': ','.join})
    df["count"] = df_ori.groupby(by=['gene_name']).agg({'mirna_name': len})['mirna_name']
    df.to_csv("gene_domain_map_id.csv")

def count_number():
    df = pd.read_csv("protein_domain_map_id.csv")
    print(len(df[["mirna_name"]].drop_duplicates()))
    print(df["count"].sum())

if __name__ == "__main__":
    create_domain_table()
    # count_number()

