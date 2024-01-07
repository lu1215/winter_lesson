import os
import sys
import sqlite3
import numpy as np
import pandas as pd
from scipy import stats
from lifelines.statistics import logrank_test
current_path = os.path.dirname(__file__)


def cal_pvalue_main(input_type, cancer, stage, high_percent, low_percent, input_pvalue):
    # input_type = request.POST["type"]
    # cancer = request.POST["cancer"]
    # stage = request.POST["stage"]
    # high_percent = request.POST["high_percent"]
    # low_percent = request.POST["low_percent"]
    # input_pvalue = request.POST["pvalue"]
    cancer = "liver_cancer|TCGA-LIHC"
    primary_site, project = cancer.split("|")
    table_name = '%s_%s_FPKM_Cufflinks'%(project,input_type)
    column_table = "%s|%s"%(stage,table_name)
    primary_key = 'gene_name' if input_type == 'genes' else 'isoform_name'
    stage_list = ['stage_1','stage_2','stage_3','stage_4']
    all_cancer_data = get_allcancer_data(column_table, input_type)
    result_list = []

def organize_and_cal_pvalue(survival_data: list, Low_Percentile:float, High_Percentile:float) -> float:
    survival_str = ""
    case_id_list = []
    # print(f"before: {survival_data}")
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
    else:
        # if error condition occur make pvalue very big 
        logrank_p_value = 1
    # print(logrank_p_value)
    return logrank_p_value, max_survival_days

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
        # cursor = connections['edward_Cufflinks'].cursor()
        primary_key = 'gene_name' if search_by == 'genes' else 'isoform_name'
        # column = primary_key
        # column = stage_dict[column_table.split('|')[0]] if column_table.split('|')[0] != 'all stage' else ','.join(stage_dict.values())
        column = f"{primary_key}, {stage_dict[column_table.split('|')[0]] if column_table.split('|')[0] != 'all stage' else ','.join(stage_dict.values())}"
        print(column)
        table_name = column_table.split('|')[1]
        cursor.execute("SELECT %s FROM `%s`"%(column,table_name))
        result = list(cursor.fetchall())
    return result