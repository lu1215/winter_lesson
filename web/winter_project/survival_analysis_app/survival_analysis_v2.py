import os
import sys
import sqlite3
import numpy as np
import pandas as pd
from scipy import stats
from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt
import matplotlib
from lifelines import KaplanMeierFitter
import base64
import io
from tqdm import tqdm, trange
import multiprocessing
import logging
matplotlib.use('Agg')
logging.getLogger('matplotlib.font_manager').disabled = True
current_path = os.path.dirname(__file__)


def cal_pvalue_main(input_type, cancer, stage, high_percent, low_percent, input_pvalue):
    cancer = "Liver_cancer|TCGA-LIHC"
    primary_site, project = cancer.split("|")
    table_name = '%s_%s_FPKM_Cufflinks'%(project,input_type)
    column_table = "%s|%s"%(stage,table_name)
    primary_key = 'gene_name' if input_type == 'genes' else 'isoform_name'
    stage_list = ['stage_1','stage_2','stage_3','stage_4']
    all_cancer_data = get_allcancer_data(column_table, input_type)
    result_list = []
    # multiprocessing
    with multiprocessing.Manager() as manager:
        result_list_m = manager.list()
        # 創建一個 multiprocessing.Pool
        pool = multiprocessing.Pool()
        # 使用 map 函數來平行處理數據
        # 注意: 如果 all_cancer_data 很大，可以考慮使用 imap 或者 imap_unordered 來節省記憶體
        pool.starmap(process_data, [(data, low_percent, high_percent, primary_site, stage, input_pvalue, result_list_m) for data in tqdm(all_cancer_data)])
        # 關閉 pool
        pool.close()
        pool.join()
        # 在這裡，result_list_m 包含所有符合條件的結果
        # print(result_list_m)
        result_list = list(result_list_m)
    str_colmns = ','.join(stage_list)

    # 2d result list to df
    df = pd.DataFrame.from_dict(result_list)
    return df, result_list

def organize_and_cal_pvalue(survival_data: list, Low_Percentile:float, High_Percentile:float, primary_site:str, survival_select:str) -> float:
    survival_str = ""
    case_id_list = []
    # print(f"before: {survival_data}")
    GT_input = survival_data[0]
    survival_data = survival_data[1:]
    # print(f"len: {len(survival_data)}")
    FPKM_list = [float(y.split("|")[0]) for x in survival_data for y in x.split(',')]
    low_quartile = np.percentile(FPKM_list, float(Low_Percentile))
    high_quartile = np.percentile(FPKM_list, 100-float(High_Percentile))
    T1 = [] #high 存活天數
    E1 = [] #high 是否死亡
    T2 = []
    E2 = []
    high_case = []
    low_case = []
    high_FPKM = []
    low_FPKM = []
    survival_days = []
    ele = "".join(survival_data).split(",")
    survival_days = [float(x.split("|")[2]) for x in ele if x.split("|")[2] != 'None']
    # print(survival_days)
    max_survival_days = max(survival_days)
    # survival_days = max(T1+T2)
    for stage in survival_data:
        for info in stage.split(','):
            FPKM = float(info.split('|')[0])
            case_id = info.split('|')[1]
            survival_times = float(info.split('|')[2]) if info.split('|')[2] != 'None' else info.split('|')[2] #存活天數
            # print(case_id,survival_times)
            survival_events = False if info.split('|')[3] == 'alive' else True #是否死亡
            if FPKM > high_quartile and (survival_times != 0 and survival_times != 'None') and survival_times <= float(max_survival_days):
                T1 += [survival_times]
                E1 += [survival_events]
                case_id_list += [case_id]
                high_case += [case_id]
                high_FPKM += [FPKM]
            elif FPKM < low_quartile and (survival_times != 0 and survival_times != 'None') and survival_times <= float(max_survival_days):
                T2 += [survival_times]
                E2 += [survival_events]
                case_id_list += [case_id]
                low_case += [case_id]
                low_FPKM += [FPKM]
    if (T2 != [] and E2 != []) and (T1 != [] and E1 != []):
        logrank_result = logrank_test(T1, T2, E1, E2)
        logrank_p_value = logrank_result.p_value
        # draw and get img info
        _, img_str = survival_plot(T1,E1,T2,E2,GT_input,primary_site,0,Low_Percentile,High_Percentile,max(T1+T2)+1,survival_select)
    else:
        # if error condition occur make pvalue very big 
        logrank_p_value = 1
        img_str = ""
    # print(logrank_p_value)
    return logrank_p_value, max_survival_days, img_str

def get_allcancer_data(column_table,search_by):
    stage_dict = {
        'stage i' : 'stage_1',
        'stage ii' : 'stage_2',
        'stage iii' : 'stage_3',
        'stage iv' : 'stage_4',
    }
    db_path = f"{current_path}/../database/db"
    with sqlite3.connect(db_path, check_same_thread=False) as db_conn:
        cursor = db_conn.cursor()
        primary_key = 'gene_name' if search_by == 'genes' else 'isoform_name'
        column = f"{primary_key}, {stage_dict[column_table.split('|')[0]] if column_table.split('|')[0] != 'all stage' else ','.join(stage_dict.values())}"
        print(column)
        table_name = column_table.split('|')[1]
        cursor.execute("SELECT %s FROM `%s`"%(column,table_name))
        result = list(cursor.fetchall())
    return result

def survival_plot(T1,E1,T2,E2,GT_input,primary_site,random_id,Low_Percentile,High_Percentile,survival_days,survival_select):
    kmf = KaplanMeierFitter()
    dpi = 100
    logrank_result = logrank_test(T1, T2, E1, E2)
    logrank_p_value = logrank_result.p_value
    logrank_test_statistic = logrank_result.test_statistic
    if T2 != [] and E2 != []:
        kmf.fit(T2, event_observed=E2, label='Low Expression (n={}, {}%)'.format(len(T2),Low_Percentile))
        ax = kmf.plot(ci_show=False,color='green',show_censors=True, figsize=(1200/dpi, 800/dpi))
    if T1 != [] and E1 != []:
        kmf.fit(T1, event_observed=E1, label='High Expression (n={}, {}%)'.format(len(T1),High_Percentile))
        kmf.plot(ax=ax,ci_show=False,color='red',show_censors=True) 
    font = {'family' : 'verdana'}
    matplotlib.rc('font', **font)
    plt.subplots_adjust(left=0.06, right=0.94, top=0.94, bottom=0.06)
    plt.title("%s"%GT_input)
    ax.text(0.1, 0.2, 'Low Percentile: {}%\nHigh Percentile: {}%\nDays: {}\nCondition: {}'.format(Low_Percentile,High_Percentile,survival_days,survival_select), horizontalalignment='left',verticalalignment='center', transform=ax.transAxes, bbox={'facecolor':'#F5F5F5', 'pad':5},fontsize=10)
    ax.text(0.1, 0.1, 'logrank pvalue: {}'.format(round(logrank_p_value,3)), horizontalalignment='left',verticalalignment='center', transform=ax.transAxes, bbox={'facecolor':'#F5F5F5', 'pad':5},fontsize=10)
    ax.set_ylim(ymin=0)
    ax.set_xlabel('Days',fontsize=12)
    ax.set_ylabel('Survival probability',fontsize=12)
    ax.grid(True)
    gridlines = ax.get_xgridlines() + ax.get_ygridlines()
    for line in gridlines:
        line.set_linestyle('-.')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    img_str = io.BytesIO()
    plt.savefig(img_str,dpi=300,format = 'png')
    plt.clf()
    plt.close()
    img_str = base64.b64encode(img_str.getvalue()).decode("utf-8").replace("\n", "")
    plt.close()
    return logrank_p_value, img_str

## function for multiprocessing star map
def process_data(data, low_percent, high_percent, primary_site, survival_select, input_pvalue, result_list):
    p_value, max_time, img_str = organize_and_cal_pvalue(data, low_percent, high_percent,  primary_site, survival_select)
    if p_value <= float(input_pvalue):
        result_list.append({"name": data[0], "logrank_p_value": p_value, "max_time": max_time, "img_str": img_str})

if __name__ == "__main__":
    print(cal_pvalue_main("genes", "Liver_cancer|TCGA-LIHC", "all stage", 50, 50, 0.05))