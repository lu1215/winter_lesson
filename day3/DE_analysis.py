#-*- coding:utf-8 -*-
#####################################################################################
# python ~/pipeline/prepare_data/cuffdiff_FPKM_gather/cuffdiff_FPKM_TEST.py \
#        ~/DIFF_output/TCGA-COAD/ \         (diff資料夾)
#        genes.read_group_tracking \        (檔案名稱)(or isoforms.read_group_tracking)
#        TCGA-COAD                          (project名稱)

# 輸出: TCGA-COAD_isoforms_FPKM(mysql), ~/EXP_output/TCGA-COAD/TCGA-COAD_isoforms_FPKM(text)
#####################################################################################
import sys
# import MySQLdb
import sqlite3
import os
import glob
from collections import OrderedDict
import itertools
import collections
import subprocess
from scipy import stats
import os
import numpy as np
import copy
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
        print('x or y is empty')
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
def main():
    folder_path = sys.argv[1]
    file_name = sys.argv[2]
    project_name = sys.argv[3]
    stage = sys.argv[4]
    level = file_name.split('.')[0]
    answer_dict = OrderedDict()
    # Permutations = list(itertools.permutations(['1','2'],2))#for debug n =1
    stage1, stage2 = stage.split(',')
    Permutations = [(f'{stage1}',f'{stage2}')]
    test_output_path = f"{current_path}/test_output/TCGA-LIHC/test_result.csv"
    # 檢查文件是否存在
    if os.path.exists(test_output_path):
        # 如果存在，刪除文件
        os.remove(test_output_path)
    f_test_output = open(f"{current_path}/test_output/TCGA-LIHC/test_result.csv",'w', encoding='utf-8')
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
            for row in f:
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

            for GI,stage_FPKM in TEST_dict.items():
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
                f_test_output.write(f"{condition_pair_upper}_{GI}, ")
                f_test_output.write(f"{','.join(map(str, TEST_result))} \n")
                
    f_test_output.close()
    table_name = "%s_%s_FPKM_Cuffdiff"%(project_name,level)
    f_FPKM = open(f"{current_path}/EXP_output/%s/%s"%(project_name,table_name),'w')
    for GI,stage_FPKM in answer_dict.items():
        insert = [GI]
        for stage,FPKM_list in stage_FPKM.items():
            FPKM_sting =  ','.join(FPKM_list)
            insert += [FPKM_sting]
        f_FPKM.write('\t'.join(insert)+"\r\n")
    f_FPKM.close()