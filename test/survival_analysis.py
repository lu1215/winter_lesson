from tqdm import tqdm
import pandas as pd
import numpy as np
import multiprocessing as mp
import os
import pandas as pd
from lifelines.statistics import logrank_test
import sqlite3
from statsmodels.stats.multitest import multipletests
current_path = os.path.dirname(__file__)


def Survival_analysis(input_data):
    gene = input_data[0]
    pat_data = input_data[1]
    high_percentile = input_data[2]
    low_percentile = input_data[3]
    patient_df = pd.DataFrame(
        columns=['Expression','Patient', 'Days', 'Status'])
    data_list = pat_data.split(',')
    for index, data in enumerate(data_list):
        data_split = data.split('|')
        patient_df.loc[index, ['Expression','Patient', 'Days', 'Status']] = [
            data_split[0], data_split[1], data_split[2], data_split[3]]
    patient_df['Days'] = pd.to_numeric(patient_df['Days'], errors='coerce')
    patient_df['Expression'] = pd.to_numeric(patient_df['Expression'], errors='coerce')
    patient_df.drop_duplicates(subset='Patient',keep='last',inplace=True)
    patient_df.drop(patient_df[patient_df['Days']== 0.0].index, inplace=True)
    patient_df.dropna(inplace=True)
    patient_df.sort_values(by='Expression',inplace=True)
    FPKM_list = np.array(patient_df['Expression'].tolist())
    low_quartile = float(np.percentile(FPKM_list, low_percentile))
    high_quartile = float(np.percentile(FPKM_list, 100-high_percentile))
    status_dict = {'alive': False, 'dead': True}
    high_patient_df = patient_df[patient_df['Expression'] >= high_quartile]
    low_patient_df = patient_df[patient_df['Expression'] <= low_quartile]
    T1 = high_patient_df['Days'].tolist()
    E1 = high_patient_df['Status'].replace(status_dict).tolist()
    T2 = low_patient_df['Days'].tolist()
    E2 = low_patient_df['Status'].replace(status_dict).tolist()
    if (T2 != [] and E2 != []) and (T1 != [] and E1 != []):
        logrank_result = logrank_test(T1, T2, E1, E2)
        logrank_p_value = logrank_result.p_value
    else:
        logrank_p_value = 1
    result_dict = {'name': gene, 'logrank_p_value': logrank_p_value}
    return (result_dict)

def cal_pvalue_main(input_type, cancer, input_stage, high_percent, low_percent, input_pvalue, cor_method):
    print(input_type, cancer, input_stage, high_percent,
          low_percent, input_pvalue, cor_method)
    high_percent = float(high_percent)
    low_percent = float(low_percent)
    primary_site, project = cancer.split("|")
    table_name = f'{project}_{input_type}_FPKM_Cufflinks'
    conn = sqlite3.connect(
        f"{current_path}/../db.sqlite3", check_same_thread=False)
    table_name = f"{project}_genes_FPKM_Cufflinks"
    stage_dict = {'stage i': 'stage_1', 'stage ii': 'stage_2',
                  'stage iii': 'stage_3', 'stage iv': 'stage_4'}
    if input_stage == 'all stage':
        query = f"SELECT * FROM '{table_name}'"
        gene_df = pd.read_sql_query(query, conn)
        # gene_list = gene_df['gene_name']
        gene_df.set_index('gene_name', inplace=True)
        new_column = gene_df.apply(lambda row: ",".join(row), axis=1)
        gene_df = gene_df.assign(patients=new_column)
        gene_df.drop(columns=["stage_1", "stage_2",
                     "stage_3", "stage_4"], inplace=True)
        gene_df.reset_index(inplace=True)
    else:
        stage = stage_dict[input_stage]
        query = f"SELECT gene_name,{stage} FROM '{table_name}'"
        gene_df = pd.read_sql_query(query, conn)
        gene_df.rename(columns={stage: 'patients'}, inplace=True)
    print(gene_df)
    df = pd.DataFrame()
    result_list = []
    processes = []
    gene_list = gene_df['gene_name'].tolist()
    data_list = gene_df['patients'].tolist()
    with mp.Pool(processes=os.cpu_count()) as pool:
        for gene, data in tqdm(zip(gene_list, data_list), desc="Assigning workload", unit="iterations"):
            input_value = [gene, data, high_percent, low_percent]
            processes.append(pool.apply_async(
                Survival_analysis, args=(input_value,)))
        result_list = [p.get() for p in tqdm(
            processes, desc="Waiting for Complete", unit="iterations")]
        df = pd.DataFrame.from_dict(x for x in tqdm(
            result_list, desc="Converting data to dataframe", unit="iterations"))
    cut_off = 0.01
    df['logrank_p_value'] = pd.to_numeric(df['logrank_p_value'], errors='coerce')
    p_value_list = df['logrank_p_value'].tolist()
    P_value_corr_FDR = multipletests(
        p_value_list, alpha=cut_off, method="fdr_bh")
    P_value_corr_Bon = multipletests(
        p_value_list, alpha=cut_off, method="bonferroni")
    df['FDR'] = P_value_corr_FDR[1]
    df['Bonferroni'] = P_value_corr_Bon[1]

    if cor_method == "None":
        df = df[df['logrank_p_value'] <= input_pvalue].reset_index(drop=True)
    else:
        df = df[df[cor_method] <= input_pvalue].reset_index(drop=True)
    return df

if __name__ == "__main__":
    result = cal_pvalue_main(
        'genes', 'Liver|TCGA-LIHC', 'all stage', 50, 50, 0.05, 'FDR')
    print(result)
