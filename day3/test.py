# -*- coding:utf-8 -*-
#####################################################################################
# python ~/pipeline/prepare_data/cuffdiff_FPKM_gather/cuffdiff_FPKM_TEST.py \
#        ~/DIFF_output/TCGA-COAD/ \         (diff資料夾)
#        genes.read_group_tracking \        (檔案名稱)(or isoforms.read_group_tracking)
#        TCGA-COAD                          (project名稱)

# 輸出: TCGA-COAD_isoforms_FPKM(mysql), ~/EXP_output/TCGA-COAD/TCGA-COAD_isoforms_FPKM(text)
#####################################################################################
import sys
# import MySQLdb
import os
import glob
from collections import OrderedDict
import itertools
import collections
import subprocess
from scipy import stats
import os

current_path = os.path.dirname(os.path.abspath(__file__))
def KS_test(x, y):
    if x and y:
        less = stats.ks_2samp(x, y, alternative='less')[1]
        greater = stats.ks_2samp(x, y, alternative='greater')[1]
        if (x == y) or (sum(x) == 0 and sum(y) == 0):  # 樣本相同 或 兩個樣本皆為0
            two_sided = 1.0
        else:
            two_sided = stats.ks_2samp(x, y, alternative='two-sided')[1]
    else:
        less = greater = two_sided = 0

    return {'two_sided': two_sided, 'greater': less, 'less': greater}

def T_test(x, y):
    if x and y:
        d, two_sided = stats.ttest_ind(x, y, equal_var=False)
        if d < 0:
            greater = 1 - two_sided / 2  # "greater" is the alternative that x has a larger mean than y
            less = two_sided / 2  # "less" is the alternative that x has a larger mean than y
        elif d > 0:
            greater = two_sided / 2
            less = 1 - two_sided / 2
        else:
            greater = less = two_sided / 2
    else:
        less = greater = two_sided = 0
    return {'two_sided': two_sided, 'greater': greater, 'less': less}

def U_test(x, y):
    if x and y:
        d, two_sided = stats.ranksums(x, y)
        if d < 0:
            greater = 1 - two_sided / 2
            less = two_sided / 2
        elif d > 0:
            greater = two_sided / 2
            less = 1 - two_sided / 2
        else:
            greater = less = two_sided / 2
    else:
        less = greater = two_sided = 0

    return {'two_sided': two_sided, 'greater': greater, 'less': less}

folder_path = sys.argv[1]
file_name = sys.argv[2]
project_name = sys.argv[3]

# conn_diff = MySQLdb.connect("localhost", "edward", "mutter12er", "edward_Cuffdiff")
# cursor_diff = conn_diff.cursor()
# conn_links = MySQLdb.connect("localhost", "edward", "mutter12er", "edward_Cufflinks")
# cursor_links = conn_links.cursor()

level = file_name.split('.')[0]

sql_column = []
condition_pair_list = []
answer_dict = OrderedDict()
# Permutations = list(itertools.permutations(['1','2','3','4'], 2))  #for TCGA-ACC
# Permutations = list(itertools.permutations(['n','1','2','3','4'], 2))
Permutations = list(itertools.permutations(['n', '1', '2', '3', '4'], 2))  # for debug n =1
for combination in Permutations:
    c1 = 'normal' if combination[0] == 'n' else 'stage_%s' % combination[0]
    c2 = 'normal' if combination[1] == 'n' else 'stage_%s' % combination[1]
    condition_pair = "%s_%s" % (combination[0], combination[1])
    condition_pair_upper = "%s_%s" % (combination[0].upper(), combination[1].upper())
    file_path = "%s%s/%s" % (folder_path, condition_pair, file_name)

    if os.path.isfile(file_path):
        TEST_dict = OrderedDict()
        with open(file_path, 'r') as f:
            for row in f:
                row = row.strip().split('\t')
                GI = row[0]
                tracking_c = row[1]
                FPKM = row[6]
                if tracking_c == 'q1':
                    if GI not in answer_dict:
                        answer_dict[GI] = OrderedDict()
                        answer_dict[GI]["%s(%s)" % (c1, condition_pair_upper)] = [FPKM]

                    else:
                        if "%s(%s)" % (c1, condition_pair_upper) not in answer_dict[GI]:
                            answer_dict[GI]["%s(%s)" % (c1, condition_pair_upper)] = [FPKM]
                        else:
                            answer_dict[GI]["%s(%s)" % (c1, condition_pair_upper)].append(FPKM)

                    if GI not in TEST_dict:
                        TEST_dict[GI] = OrderedDict()
                        TEST_dict[GI]["%s(%s)" % (c1, condition_pair_upper)] = [FPKM]

                    else:
                        if "%s(%s)" % (c1, condition_pair_upper) not in TEST_dict[GI]:
                            TEST_dict[GI]["%s(%s)" % (c1, condition_pair_upper)] = [FPKM]
                        else:
                            TEST_dict[GI]["%s(%s)" % (c1, condition_pair_upper)].append(FPKM)

                if tracking_c == 'q2':
                    if GI not in answer_dict:
                        answer_dict[GI] = OrderedDict()
                        answer_dict[GI]["%s(%s)" % (c2, condition_pair_upper)] = [FPKM]

                    else:
                        if "%s(%s)" % (c2, condition_pair_upper) not in answer_dict[GI]:
                            answer_dict[GI]["%s(%s)" % (c2, condition_pair_upper)] = [FPKM]
                        else:
                            answer_dict[GI]["%s(%s)" % (c2, condition_pair_upper)].append(FPKM)

                    if GI not in TEST_dict:
                        TEST_dict[GI] = OrderedDict()
                        TEST_dict[GI]["%s(%s)" % (c2, condition_pair_upper)] = [FPKM]

                    else:
                        if "%s(%s)"%(c2,condition_pair_upper) not in TEST_dict[GI]:
                            TEST_dict[GI]["%s(%s)"%(c2,condition_pair_upper)] = [FPKM]
                        else:
                            TEST_dict[GI]["%s(%s)"%(c2,condition_pair_upper)] += [FPKM]

            diff_table_name = "%s_%s_%s" % (project_name, condition_pair_upper, level)
            sql_add_column = """ALTER TABLE `%s` ADD COLUMN `KS_test_twosided` varchar(100),
                                                ADD COLUMN `KS_test_greater` varchar(100),
                                                ADD COLUMN `KS_test_less` varchar(100),
                                                ADD COLUMN `T_test_twosided` varchar(100),
                                                ADD COLUMN `T_test_greater` varchar(100),
                                                ADD COLUMN `T_test_less` varchar(100),
                                                ADD COLUMN `U_test_twosided` varchar(100),
                                                ADD COLUMN `U_test_greater` varchar(100),
                                                ADD COLUMN `U_test_less` varchar(100)""" % diff_table_name

            for GI, stage_FPKM in TEST_dict.items():
                # cursor_diff.execute("SHOW COLUMNS FROM `%s`" % diff_table_name)
                # detect_column = [column[0] for column in cursor_diff.fetchall()]
                # if len(detect_column) == 14:
                    # cursor_diff.execute(sql_add_column)

                TEST_arg_FPKM = []
                for stage, FPKM_list in stage_FPKM.items():
                    TEST_arg_FPKM.append(list(map(float, FPKM_list)))

                print(condition_pair_upper, GI)
                # KS_test_result = KS_test(TEST_arg_FPKM[0],TEST_arg_FPKM[1])
                # T_test_result = T_test(TEST_arg_FPKM[0],TEST_arg_FPKM[1])
                # U_test_result = U_test(TEST_arg_FPKM[0],TEST_arg_FPKM[1])
                print('--------', TEST_arg_FPKM)
                KS_test_result = KS_test(TEST_arg_FPKM[1], TEST_arg_FPKM[0])
                T_test_result = T_test(TEST_arg_FPKM[1], TEST_arg_FPKM[0])
                U_test_result = U_test(TEST_arg_FPKM[1], TEST_arg_FPKM[0])
                # TEST_result = KS_test_result.values() + T_test_result.values() + U_test_result.values()
                TEST_result = []
                TEST_result.extend(KS_test_result.values())
                TEST_result.extend(T_test_result.values())
                TEST_result.extend(U_test_result.values())
                diff_update_sql = """UPDATE `%s` SET `KS_test_twosided`='%s',`KS_test_greater`='%s',`KS_test_less`='%s',
                                                    `T_test_twosided`='%s',`T_test_greater`='%s',`T_test_less`='%s',
                                                    `U_test_twosided`='%s',`U_test_greater`='%s',`U_test_less`='%s' 
                                                    WHERE `test_id`='%s'""" % (diff_table_name,
                                                                            TEST_result[0], TEST_result[1], TEST_result[2],
                                                                            TEST_result[3], TEST_result[4], TEST_result[5],
                                                                            TEST_result[6], TEST_result[7], TEST_result[8],
                                                                            GI)
                # cursor_diff.execute(diff_update_sql)
                # conn_diff.commit()

                print(TEST_result)

    if "%s(%s)" % (c1, condition_pair_upper) not in sql_column:
        sql_column.append("%s(%s)" % (c1, condition_pair_upper))

    if "%s(%s)" % (c2, condition_pair_upper) not in sql_column:
        sql_column.append("%s(%s)" % (c2, condition_pair_upper))
# conn_diff.close()

table_name = "%s_%s_FPKM_Cuffdiff" % (project_name, level)
sql_column = ['`%s_name` varchar(100)' % level[:-1]] + list(map(lambda x: '`' + x + '` LONGTEXT', sql_column))
sql_column = ','.join(sql_column)

# cursor_links.execute("DROP TABLE IF EXISTS `%s`" % (table_name))
# cursor_links.execute("CREATE TABLE `%s`(%s,primary key(`%s`))" % (table_name, sql_column, level[:-1] + "_name"))
# cursor_links.execute("ALTER TABLE `%s` ADD INDEX gene_isoform(`%s`)" % (table_name, level[:-1] + "_name"))
f_out = open(f"{current_path}/EXP_output/%s/%s"%(project_name,table_name),'w')
for GI, stage_FPKM in answer_dict.items():
    insert = [GI]
    for stage, FPKM_list in stage_FPKM.items():
        FPKM_sting = ','.join(FPKM_list)
        insert += [FPKM_sting]
    f_out.write('\t'.join(insert) + "\r\n")
#     cursor_links.execute("INSERT INTO `%s` VALUES%s" % (table_name, tuple(insert,)))
#     conn_links.commit()
# conn_links.close()
f_out.close()