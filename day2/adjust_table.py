import pandas as pd
import os

current_path = os.path.dirname(os.path.abspath(__file__))
df = pd.read_csv(f"{current_path}/TCGA-LIHC_genes_FPKM_Cufflinks.csv")
df = df[["gene_name", "stage_1", "stage_2", "stage_3", "stage_4"]]
df.to_csv(f"{current_path}/TCGA-LIHC_genes_FPKM_Cufflinks_v2.csv", index=False)