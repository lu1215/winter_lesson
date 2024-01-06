#-*- coding:utf-8 -*-
#####################################################################################
# python ~/pipeline/prepare_data/cuffdiff_FPKM_gather/cuffdiff_FPKM_TEST.py \
#        ~/DIFF_output/TCGA-COAD/ \         (diff資料夾)
#        genes.read_group_tracking \        (檔案名稱)(or isoforms.read_group_tracking)
#        TCGA-COAD                          (project名稱)

# 輸出: TCGA-COAD_isoforms_FPKM(mysql), ~/EXP_output/TCGA-COAD/TCGA-COAD_isoforms_FPKM(text)
#####################################################################################
import sys
import os
from collections import OrderedDict
from scipy import stats
import os
import numpy as np
from statsmodels.stats.multitest import multipletests
import pandas as pd
import math
from tqdm import tqdm
from itertools import islice
current_path = os.path.dirname(os.path.abspath(__file__))


def KS_test(x,y):
    if x.size > 0 and y.size > 0:
        less = stats.mstats.ks_2samp(x, y,alternative = 'less')[1] 
        greater = stats.mstats.ks_2samp(x, y,alternative = 'greater')[1]
        if (x.size == y.size) or (sum(x) == 0 and sum(y) == 0): #樣本相同 或 兩個樣本皆為0
            two_sided = 1.0
        else:
            two_sided = stats.mstats.ks_2samp(x, y,alternative = 'two-sided')[1]
    else:
        less = greater = two_sided = 0

    return {'two_sided':two_sided, 'greater':less, 'less':greater}
    
def T_test(x,y):
    if x.size > 0 and y.size > 0:
        d, two_sided = stats.ttest_ind(x, y, equal_var=False)
        if d < 0:
            greater = 1 - two_sided/2 #"greater" is the alternative that x has a larger mean than y
            less = two_sided/2 #"less" is the alternative that x has a larger mean than y
        elif d >0:
            greater = two_sided/2
            less = 1 - two_sided/2
        else:
            greater = less = two_sided/2
    else:
        less = greater = two_sided = 0
    return {'two_sided':two_sided, 'greater':greater, 'less':less}

def U_test(x,y):
    if x.size > 0 and y.size > 0:
        d, two_sided = stats.ranksums(x, y)
        if d < 0:
            greater = 1 - two_sided/2
            less = two_sided/2
        elif d >0:
            greater = two_sided/2
            less = 1 - two_sided/2
        else:
            greater = less = two_sided/2
    else:
        less = greater = two_sided = 0
        
    return {'two_sided':two_sided, 'greater':greater, 'less':less}

def main(folder_path, file_name, project_name, stage):
    # folder_path = sys.argv[1]
    # file_name = sys.argv[2]
    # project_name = sys.argv[3]
    # stage = sys.argv[4]
    total_GI = []
    total_KS_test_two_sided = []
    total_KS_test_greater = []
    total_KS_test_less = []
    total_T_test_two_sided = []
    total_T_test_greater = []
    total_T_test_less = []
    total_U_test_two_sided = []
    total_U_test_greater = []
    total_U_test_less = []
    level = file_name.split('.')[0]
    answer_dict = OrderedDict()
    # Permutations = list(itertools.permutations(['1','2'],2))#for debug n =1
    stage1, stage2 = stage.split('_')
    Permutations = [(f'{stage1}',f'{stage2}')]
    test_output_path = f"{current_path}/test_output/TCGA-LIHC/{stage}test_result.csv"
    # 檢查文件是否存在
    if os.path.exists(test_output_path):
        # 如果存在，刪除文件
        os.remove(test_output_path)
    f_test_output = open(f"{current_path}/test_output/TCGA-LIHC/{stage}test_result.csv",'w', encoding='utf-8')
    f_test_output.write("gene_name, KS_test_two_sided, KS_test_greater, KS_test_less, T_test_two_sided, T_test_greater, T_test_less, U_test_two_sided, U_test_greater, U_test_less \n")
    for combination in Permutations:
        c1 = 'normal' if combination[0] == 'n' else 'stage_%s'%combination[0]
        c2 = 'normal' if combination[1] == 'n' else 'stage_%s'%combination[1]
        condition_pair = "%s_%s"%(combination[0],combination[1])
        condition_pair_upper = "%s_%s"%(combination[0].upper(),combination[1].upper())
        file_path = "%s%s/%s"%(folder_path,condition_pair,file_name)
        if os.path.isfile(file_path) == True:
            TEST_dict = OrderedDict()
            f = open(file_path,'r')
            for row in tqdm(f):
                row = row.strip().split('\t')
                GI = row[0]
                tracking_c = row[1]
                FPKM = row[6]
                if tracking_c == 'q1':
                    if GI not in answer_dict:
                        answer_dict[GI] = OrderedDict()
                        answer_dict[GI]["%s(%s)"%(c1,condition_pair_upper)] = [FPKM]
                    else:
                        if "%s(%s)"%(c1,condition_pair_upper) not in answer_dict[GI]:
                            answer_dict[GI]["%s(%s)"%(c1,condition_pair_upper)] = [FPKM]
                        else:
                            answer_dict[GI]["%s(%s)"%(c1,condition_pair_upper)] += [FPKM]
                    if GI not in TEST_dict:
                        TEST_dict[GI] = OrderedDict()
                        TEST_dict[GI]["%s(%s)"%(c1,condition_pair_upper)] = [FPKM]
                    else:
                        if "%s(%s)"%(c1,condition_pair_upper) not in TEST_dict[GI]:
                            TEST_dict[GI]["%s(%s)"%(c1,condition_pair_upper)] = [FPKM]
                        else:
                            TEST_dict[GI]["%s(%s)"%(c1,condition_pair_upper)] += [FPKM]
                if tracking_c == 'q2':
                    if GI not in answer_dict:
                        answer_dict[GI] = OrderedDict()
                        answer_dict[GI]["%s(%s)"%(c2,condition_pair_upper)] = [FPKM]
                    else:
                        if "%s(%s)"%(c2,condition_pair_upper) not in answer_dict[GI]:
                            answer_dict[GI]["%s(%s)"%(c2,condition_pair_upper)] = [FPKM]
                        else:
                            answer_dict[GI]["%s(%s)"%(c2,condition_pair_upper)] += [FPKM]
                    if GI not in TEST_dict:
                        TEST_dict[GI] = OrderedDict()
                        TEST_dict[GI]["%s(%s)"%(c2,condition_pair_upper)] = [FPKM]
                    else:
                        if "%s(%s)"%(c2,condition_pair_upper) not in TEST_dict[GI]:
                            TEST_dict[GI]["%s(%s)"%(c2,condition_pair_upper)] = [FPKM]
                        else:
                            TEST_dict[GI]["%s(%s)"%(c2,condition_pair_upper)] += [FPKM]
            diff_table_name = "%s_%s_%s"%(project_name,condition_pair_upper,level)

            for GI,stage_FPKM in tqdm(TEST_dict.items()):
            # for GI,stage_FPKM in tqdm(islice(TEST_dict.items(), 100)):
                
                TEST_arg_FPKM = []
                for stage,FPKM_list in stage_FPKM.items():
                    TEST_arg_FPKM += [map(float,FPKM_list)]
                # print (condition_pair_upper,GI)
                parameter1 = np.fromiter(TEST_arg_FPKM[1], dtype=float)
                parameter2 = np.fromiter(TEST_arg_FPKM[0], dtype=float)
                KS_test_result = KS_test(parameter1,parameter2)
                T_test_result = T_test(parameter1,parameter2)
                U_test_result = U_test(parameter1,parameter2)
                TEST_result = []
                TEST_result.extend(KS_test_result.values())
                TEST_result.extend(T_test_result.values())
                TEST_result.extend(U_test_result.values())
                
                total_GI.append(GI)
                total_KS_test_two_sided.append(TEST_result[0])
                total_KS_test_greater.append(TEST_result[1])
                total_KS_test_less.append(TEST_result[2])
                total_T_test_two_sided.append(TEST_result[3])
                total_T_test_greater.append(TEST_result[4])
                total_T_test_less.append(TEST_result[5])
                total_U_test_two_sided.append(TEST_result[6])
                total_U_test_greater.append(TEST_result[7])
                total_U_test_less.append(TEST_result[8])
                f_test_output.write(f"{GI}, ")
                f_test_output.write(f"{','.join(map(str, TEST_result))} \n")
    f_test_output.close()
    ## correction
    cut_off = 1
    cor_method = "fdr_bh"
    # cor_method = "bonferroni"
    total_KS_test_two_sided = [2 if math.isnan(x) else x for x in total_KS_test_two_sided]
    total_KS_test_greater = [2 if math.isnan(x) else x for x in total_KS_test_greater]
    total_KS_test_less = [2 if math.isnan(x) else x for x in total_KS_test_less]
    total_T_test_two_sided = [2 if math.isnan(x) else x for x in total_T_test_two_sided]
    total_T_test_greater = [2 if math.isnan(x) else x for x in total_T_test_greater]
    total_T_test_less = [2 if math.isnan(x) else x for x in total_T_test_less]
    total_U_test_two_sided = [2 if math.isnan(x) else x for x in total_U_test_two_sided]
    total_U_test_greater = [2 if math.isnan(x) else x for x in total_U_test_greater]
    total_U_test_less = [2 if math.isnan(x) else x for x in total_U_test_less]
    df = pd.DataFrame()
    # print(total_T_test_two_sided)
    df['gene_name'] = total_GI
    df['KS_test_two_sided'] = multipletests(total_KS_test_two_sided,alpha=cut_off, method= cor_method)[1]
    df['KS_test_greater'] = multipletests(total_KS_test_greater,alpha=cut_off, method= cor_method)[1]
    df['KS_test_less'] = multipletests(total_KS_test_less,alpha=cut_off, method= cor_method)[1]
    df['T_test_two_sided'] = multipletests(total_T_test_two_sided,alpha=cut_off, method= cor_method)[1]
    df['T_test_greater'] = multipletests(total_T_test_greater,alpha=cut_off, method= cor_method)[1]
    df['T_test_less'] = multipletests(total_T_test_less,alpha=cut_off, method= cor_method)[1]
    df['U_test_two_sided'] = multipletests(total_U_test_two_sided,alpha=cut_off, method= cor_method)[1]
    df['U_test_greater'] = multipletests(total_U_test_greater,alpha=cut_off, method= cor_method)[1]
    df['U_test_less'] = multipletests(total_U_test_less,alpha=cut_off, method= cor_method)[1]
    # df["first_FPKM"]
    # df["second_FPKM"]

    # P_value_corr_FDR = multipletests(total_KS_test_two_sided,alpha=cut_off, method= "fdr_bh")
    # P_value_corr_Bon = multipletests(total_KS_test_two_sided,alpha=cut_off, method= cor_method)
    
    # P_value_corr_FDR = multipletests(test,alpha=cut_off, method= "fdr_bh")
    # P_value_corr_Bon = multipletests(test,alpha=cut_off, method= "bonferroni")
    table_name = "%s_%s_FPKM_Cuffdiff"%(project_name,level)
    f_FPKM = open(f"{current_path}/EXP_output/%s/%s"%(project_name,table_name),'w')
    for GI,stage_FPKM in answer_dict.items():
        insert = [GI]
        for stage,FPKM_list in stage_FPKM.items():
            FPKM_sting =  ','.join(FPKM_list)
            insert += [FPKM_sting]
        f_FPKM.write('\t'.join(insert)+"\r\n")
    f_FPKM.close()
    
    df_FPKM = pd.read_csv(f"{current_path}/EXP_output/%s/%s"%(project_name,table_name), sep='\t' , header=None)
    df_FPKM = df_FPKM.rename(columns={0: 'gene', 1: 'first_FPKM', 2: 'second_FPKM'})
    print(df_FPKM.columns)
    df["avg_f_FPKM"] = df_FPKM["first_FPKM"].apply(lambda x: np.mean(list(map(float, x.split(',')))))
    df["avg_s_FPKM"] = df_FPKM["second_FPKM"].apply(lambda x: np.mean(list(map(float, x.split(',')))))
    df["foldchange"] = df["avg_s_FPKM"] / df["avg_f_FPKM"]
    df.to_csv(f"{current_path}/test_output/TCGA-LIHC/corr_test_result.csv", index=False)

if __name__ == '__main__':
    folder_path = sys.argv[1]
    file_name = sys.argv[2]
    project_name = sys.argv[3]
    stage = sys.argv[4]
    main(folder_path, file_name, project_name, stage)