import pandas as pd

df = pd.read_csv('hg38_gene_transcripts_20180130.csv')
df = df[['gene']]
df.to_csv('all_genes.csv', index=False)