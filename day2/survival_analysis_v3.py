from tqdm import tqdm
import pandas as pd
import numpy as np
import multiprocessing
import logging
import os
import pandas as pd
from lifelines.statistics import logrank_test
import sqlite3
import io
import base64
import matplotlib
from lifelines import KaplanMeierFitter
from matplotlib import pyplot as plt

current_path = os.path.dirname(__file__)
stage_list = ['stage_1','stage_2','stage_3','stage_4']
stage_dict = {
    'stage i' : 'stage_1',
    'stage ii' : 'stage_2',
    'stage iii' : 'stage_3',
    'stage iv' : 'stage_4',
}

def cal_pvalue_main(input_type, cancer, stage, high_percent, low_percent, input_pvalue):
    # input_type = request.POST["type"]
    # cancer = request.POST["cancer"]
    # stage = request.POST["stage"]
    # high_percent = request.POST["high_percent"]
    # low_percent = request.POST["low_percent"]
    # input_pvalue = request.POST["pvalue"]
    primary_site, project = cancer.split("|")
    table_name = '%s_%s_FPKM_Cufflinks'%(project,input_type)
    column_table = "%s|%s"%(stage,table_name)
    primary_key = 'gene_name' if input_type == 'genes' else 'isoform_name'
    stage_list = ['stage_1','stage_2','stage_3','stage_4']
    all_cancer_data = get_allcancer_data(column_table, input_type)
    result_list = []

    ## multiprocessing
    # with multiprocessing.Manager() as manager:
    #     result_list_m = manager.list()
    #     # 創建一個 multiprocessing.Pool
    #     pool = multiprocessing.Pool()
    #     # 使用 map 函數來平行處理數據
    #     # 注意: 如果 all_cancer_data 很大，可以考慮使用 imap 或者 imap_unordered 來節省記憶體
    #     pool.starmap(process_data, [(data, low_percent, high_percent, input_pvalue, result_list_m) for data in tqdm(all_cancer_data)])
    #     # 關閉 pool
    #     pool.close()
    #     pool.join()
    #     # 在這裡，result_list_m 包含所有符合條件的結果
    #     # print(result_list_m)
    #     result_list = list(result_list_m)

    ## original
    # for i in tqdm(range(len(all_cancer_data))):
    for i in tqdm(range(1000)):
        p_value, max_time = organize_and_cal_pvalue(all_cancer_data[i], low_percent, high_percent)
        if p_value <= float(input_pvalue):
            result_list.append({"name":all_cancer_data[i][0], "logrank_p_value":p_value, "max_time": max_time})
    # print(f"result_list: {result_list}")
    str_colmns = ','.join(stage_list)

    # 2d result list to df
    df = pd.DataFrame.from_dict(result_list)
    return df, result_list

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
        logrank_p_value = 1
    return logrank_p_value, max_survival_days

def process_data(data, low_percent, high_percent, input_pvalue, result_list):
    p_value, max_time = organize_and_cal_pvalue(data, low_percent, high_percent)
    if p_value <= float(input_pvalue):
        result_list.append({"name": data[0], "logrank_p_value": p_value, "max_time": max_time})

def survival_plot_realtime(project, primary_site, search_by, GT_input,random_id, Low_Percentile, High_Percentile, survival_select, survival_days:int =None):
    if survival_days==None:
        survival_days = survival_max_days(project, GT_input, search_by, survival_select)["max_survival_days"]
    # project = request.POST['project']
    # primary_site = request.POST['primary_site']
    # search_by = request.POST['search_by']
    # GT_input = request.POST['GT_input']
    # random_id = request.POST['random_id']
    # Low_Percentile = request.POST['Low_Percentile']
    # High_Percentile = request.POST['High_Percentile']
    # survival_days = request.POST['survival_days']
    # survival_select = request.POST['survival_select']
    table_name = '%s_%s_FPKM_Cufflinks'%(project,search_by)
    column_table = "%s|%s"%(survival_select,table_name)
 #### patched by t50504
    survival_data = Survival_plot.survival_data_realtime(column_table,search_by,GT_input)
    survival_str = ""
    case_id_list = []
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
    for stage in survival_data:
        for info in stage.split(','):
            FPKM = float(info.split('|')[0])
            case_id = info.split('|')[1]
            survival_times = float(info.split('|')[2]) if info.split('|')[2] != 'None' else info.split('|')[2] #存活天數
            # print(case_id,survival_times)
            survival_events = False if info.split('|')[3] == 'alive' else True #是否死亡
            if FPKM > high_quartile and (survival_times != 0 and survival_times != 'None') and survival_times <= float(survival_days):
                T1 += [survival_times]
                E1 += [survival_events]
                case_id_list += [case_id]
                high_case += [case_id]
                high_FPKM += [FPKM]
            elif FPKM < low_quartile and (survival_times != 0 and survival_times != 'None') and survival_times <= float(survival_days):
                T2 += [survival_times]
                E2 += [survival_events]
                case_id_list += [case_id]
                low_case += [case_id]
                low_FPKM += [FPKM]
    if (T2 != [] and E2 != []) and (T1 != [] and E1 != []):
        _, img_str = Survival_plot.survival_plot(T1,E1,T2,E2,GT_input,primary_site,random_id,Low_Percentile,High_Percentile,max(T1+T2),survival_select)
        # survival_download = Survival_plot.survival_download(T1,E1,T2,E2,high_case,low_case,high_FPKM,low_FPKM,GT_input,primary_site,Low_Percentile,High_Percentile,'all stage')
        # survival_csv = survival_download(T1,E1,T2,E2,high_case,low_case,high_FPKM,low_FPKM,GT_input,primary_site, random_id, Low_Percentile,High_Percentile,survival_select)
    else:
        survival_str = ' Survival analysis is not available for '+GT_input+' since more than half of the samples have zero expression.'
        img_str = ''
    return img_str

def survival_max_days(project: str, GT_input: str, search_by: str, survival_select: str) -> dict:
    table_name = '%s_%s_FPKM_Cufflinks'%(project,search_by)
    column_table = "%s|%s"%(survival_select,table_name)
    survival_data = "".join(Survival_plot.survival_data_realtime(column_table,search_by,GT_input)).split(",")
    # print(f"survival_data: {survival_data}")
    survival_days = [float(x.split("|")[2]) for x in survival_data if x.split("|")[2] != 'None']
    max_survival_days = max(survival_days)
    return {"max_survival_days": max_survival_days}

class Survival_plot():
    def survival_data_realtime(column_table,search_by,GT_input):
        table_name = column_table.split('|')[1]
        db_path = f"{current_path}/../database/db"
        with sqlite3.connect(db_path, check_same_thread=False) as db_conn:
            cursor = db_conn.cursor() 
            # cursor = connections['edward_Cufflinks'].cursor()
            primary_key = 'gene_name' if search_by == 'genes' else 'isoform_name'
            if column_table.split('|')[0] != 'all stage':
                column = stage_dict[column_table.split('|')[0]] if column_table.split('|')[0] != 'all stage' else ','.join(stage_dict.values())
            else:
                cursor.execute("SHOW columns FROM `{}`".format(table_name))
                data_columns = [column[0] for column in cursor.fetchall()]
                data_columns = data_columns[1:]
                column=[]
                for i in data_columns:
                    if i != 'normal' and i !='all_stage' and i != 'pvalue(normal_allstage)' and i != 'pvalue(allstage_normal)':
                        column.append('`'+i+'`')
                column = ','.join(column)
            print("SELECT %s FROM `%s` WHERE `%s` = '%s'"%(column,table_name,primary_key,GT_input))
            cursor.execute("SELECT %s FROM `%s` WHERE `%s` = '%s'"%(column,table_name,primary_key,GT_input))
            print("SELECT %s FROM `%s` WHERE `%s` = '%s'"%(column,table_name,primary_key,GT_input))
            result = list(cursor.fetchall()[0])
        return result
        
    def survival_plot(T1,E1,T2,E2,GT_input,primary_site,random_id,Low_Percentile,High_Percentile,survival_days,survival_select):
        plot_path = f"{current_path}/../static/Survival_plot/survival_plot_%s_%s.png"%(primary_site.replace(' ','_').replace("(",'').replace(")",''),random_id)
        try:
            os.remove(plot_path)
        except OSError:
            pass
        kmf = KaplanMeierFitter()
        dpi = 100
        # T1 = [232,565,205,4967,370,276,59,1649,819,590,82,522,3364,640,211,223,783,1556,3981,1005,638,1949,615,823,484,1884,899,588,539,562,344,88,413,789,273,311,1639,510,278,37,20,122,851,118,2423,1048,540,547,685,610,237,89,406,603,572,272,415,303,377,599,328,248,272,508,539,62,579,1108,2109,251,974,384,154,633,173,2027,467,1110,382,5041,1582,1830,491,128,90,372,474,522,200,4343,253,272,480,646,758,408,400,251,2828,1947]
        # E1 = [True,True,True,False,True,False,False,False,True,True,False,False,False,False,True,True,False,True,False,True,False,False,True,True,False,False,False,False,False,False,True,True,True,False,True,True,False,True,True,False,True,True,False,True,False,False,False,True,True,False,True,False,True,False,False,True,True,True,False,True,True,True,True,True,True,True,True,False,False,True,True,False,True,False,True,False,True,False,False,False,False,False,False,True,True,False,True,True,True,False,True,True,False,False,False,True,True,False,True,False] #是否死亡
        # T2 = [13,1326,906,474,1090,76,469,2177,95,55,680,544,149,544,144,2008,542,129,2954,734,415,460,182,385,467,1845,324,578,949,240,3183,399,20,262,612,2380,361,453,428,370,168,1029,1912,366,369,105,1561,642,2703,859,333,410,1348,3011,1350,455,700,921,1792,690,641,832,577,798,324,512,117,565,691,57,1869,1064,389,415,384,1454,28,1072,1621,368,385,2139,864,56,2020,345,1604,1804,365,388,508,364,224,2312,142,997,1522,67,294,1529]
        # E2 = [False,False,False,False,False,True,False,False,False,False,True,True,True,True,True,False,False,False,True,True,False,True,True,True,True,False,True,False,True,False,True,False,False,True,True,False,False,True,False,False,True,False,False,False,False,False,False,False,False,True,False,False,True,False,False,False,False,False,False,True,False,False,True,False,True,False,False,True,False,True,True,True,False,False,False,False,False,False,False,False,True,False,False,True,False,False,False,True,False,True,False,False,False,False,True,False,False,False,True,False]
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
        plt.show()
        plt.savefig(plot_path,dpi=dpi)
        img_str = io.BytesIO()
        plt.savefig(img_str,dpi=300,format = 'png')
        plt.clf()
        plt.close()
        img_str = base64.b64encode(img_str.getvalue()).decode("utf-8").replace("\n", "")
        return logrank_p_value, img_str

    def survival_download(T1,E1,T2,E2,high_case,low_case,high_FPKM,low_FPKM,GT_input,primary_site,random_id,Low_Percentile,High_Percentile,survival_select):
        output_data = [
                        ['Query: %s'%(GT_input)],
                        ['Primary site: %s'%(primary_site)],
                        ['Low Percentile: %s'%(Low_Percentile)],
                        ['High Percentile: %s'%(High_Percentile)],
                        ['Condition: %s'%(survival_select)],
                        [],
                        ['Patient','Days','Status','Expression','Group']
        ]
        for idx,row in enumerate(low_case):
            Status = 'Alive' if E2[idx] == False else 'Dead'
            output_data += [[row,T2[idx],Status,low_FPKM[idx],'Low']]
        for idx,row in enumerate(high_case):
            Status = 'Alive' if E1[idx] == False else 'Dead'
            output_data += [[row,T1[idx],Status,high_FPKM[idx],'High']]
        print('survival_download')
        return output_data

if __name__ == "__main__":
    df_pvalue, _ = cal_pvalue_main()
