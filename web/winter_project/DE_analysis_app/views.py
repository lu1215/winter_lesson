from django.shortcuts import render
import pandas as pd

# Create your views here.
def DE_page(request):
    return render(request, 'DE_analysis.html', locals())

def DE_cal(request):
    DE_level = request.POST["type"]
    cancer = request.POST["cancer"]
    primary_site, project = cancer.split("|")
    DE_filter = request.POST.getlist("DE_filter[]")
    print(DE_filter)
    ########## need to check condition1 and condition2 column data ##########
    table_data, table_col = DEscreener_getdata(DE_level, primary_site, project, DE_filter)
    df_DE = pd.DataFrame(table_data, columns=table_col)
    print(table_col)
    # df_DE = pd.read_csv(f"{current_path}/../static/data/DE_data/{selected_condition}.csv", skiprows=8)
    # print(df_DE)
    df_DE[list(df_DE.columns)[-2]] = df_DE[list(df_DE.columns)[-2]].apply(lambda x: round(1/x, 5))
    replace_col_name = list(df_DE.columns)[-1]
    print(list(df_DE.columns)[-1])
    df_DE.rename(columns={replace_col_name: "q-value"}, inplace=True)
    # print(df_DE)
    lefton = "name" if "name" in df.columns else "gene_name"
    df = df_DE if len(df) == 0 else pd.merge(df, df_DE, left_on=lefton, right_on="gene_name", how="inner")

def DEscreener_getdata(DE_level, primary_site, project, DE_filter):
    condition1 = DE_filter[0].split('|')[0]
    condition2 = DE_filter[1].split('|')[0]
    condition1_count = DE_filter[0].split('|')[2]
    condition2_count = DE_filter[1].split('|')[2]
    FC_select = DE_filter[2].strip()
    FC_input = DE_filter[3].strip()
    TEST_select = DE_filter[4]
    TESTstates_select = DE_filter[5].split(' (')[0].strip()
    TEST_input = DE_filter[6].strip()
    FC_input = float(FC_input) if FC_input else FC_input
    TEST_input = float(TEST_input) if TEST_input else TEST_input
    stage_to_num_dict = {
        'normal' :'N',
		'stage i' : '1',
		'stage ii' : '2',
		'stage iii' : '3',
		'stage iv' : '4',
	}
    table_name_cuffdiff = "%s_%s_%s_%s"%(project,stage_to_num_dict[condition1],stage_to_num_dict[condition2],DE_level)
    print(table_name_cuffdiff)
    if TEST_select == 'Cuffdiff DE test':
        TEST_column = "q_value"
    else:
        TEST_column = "%s_%s"%(TEST_select.replace(' ','_'),TESTstates_select.replace(' ','')) if TEST_select else ''
    diff_data,download_table_data, table_column = filter_logFC_Pvalue(TEST_column,table_name_cuffdiff,FC_select,FC_input,TEST_input,DE_level)
    if len(diff_data) != 0:
        return diff_data, table_column
    else:
        return [], table_column
    
def filter_logFC_Pvalue(TEST_column,table_name_cuffdiff,FC_select,FC_input,TEST_input,DE_level="genes"):
    # cursor = connections['edward_Cuffdiff'].cursor()
    # start = time.time()
    
    if TEST_input:
        if DE_level == 'genes':
            sql = "SELECT `test_id`,`value_1`,`value_2`,`log2(fold_change)`,`%s`,ROUND( (`value_1` / `value_2`), 3 ) AS `fold_change` FROM `%s`"%(TEST_column,table_name_cuffdiff)
        elif DE_level == 'isoforms':	
            sql = "SELECT `test_id`,`gene_id`,`value_1`,`value_2`,`log2(fold_change)`,`%s`,ROUND( (`value_1` / `value_2`), 3 ) AS `fold_change` FROM `%s`"%(TEST_column,table_name_cuffdiff)
    else:
        if DE_level == 'genes':
            sql = "SELECT `test_id`,`value_1`,`value_2`,`log2(fold_change)`,ROUND( (`value_1` / `value_2`), 3 ) AS `fold_change` FROM `%s`"%(table_name_cuffdiff)
        elif DE_level == 'isoforms':	
            sql = "SELECT `test_id`,`gene_id`,`value_1`,`value_2`,`log2(fold_change)`,ROUND( (`value_1` / `value_2`), 3 ) AS `fold_change` FROM `%s`"%(table_name_cuffdiff)

    if FC_input != '' and TEST_input == '':
        if FC_select == '≥':
            sql += " WHERE (ROUND(`value_1`,3) >= 0.001 AND ROUND(`value_2`,3) >= 0.001) AND (ROUND(`value_2`,3)/ROUND(`value_1`,3) >= %s )"%(FC_input)
        elif FC_select == '≤':
            sql += " WHERE (ROUND(`value_1`,3) >= 0.001 AND ROUND(`value_2`,3) >= 0.001) AND (ROUND(`value_2`,3)/ROUND(`value_1`,3) <= %s )"%(FC_input)
        elif FC_select == 'Either':
            sql += " WHERE (ROUND(`value_1`,3) >= 0.001 AND ROUND(`value_2`,3) >= 0.001) AND (pow(2,`log2(fold_change)`) >= %s OR pow(2,`log2(fold_change)`) <= %s OR `log2(fold_change)` = 'inf' OR `log2(fold_change)` = '-inf')"%(FC_input,FC_input)

    elif FC_input == '' and TEST_input != '':
        
        sql += " WHERE (ROUND(`value_1`,3) >= 0.001 AND ROUND(`value_2`,3) >= 0.001) AND `%s` <= %s"%(TEST_column,TEST_input)

    elif FC_input != '' and TEST_input != '':
        
        if FC_select == '≥':
            sql += " WHERE (ROUND(`value_1`,3) >= 0.001 AND ROUND(`value_2`,3) >= 0.001) AND (((ROUND(`value_2`,3)/ROUND(`value_1`,3)) >= %s ) AND `%s` <= %s)"%(FC_input,TEST_column,TEST_input)
        elif FC_select == '≤':
            sql += " WHERE (ROUND(`value_1`,3) >= 0.001 AND ROUND(`value_2`,3) >= 0.001) AND (((ROUND(`value_2`,3)/ROUND(`value_1`,3)) <= %s ) AND `%s` <= %s)"%(FC_input,TEST_column,TEST_input)
        elif FC_select == 'Either':
            sql += " WHERE (ROUND(`value_1`,3) >= 0.001 AND ROUND(`value_2`,3) >= 0.001) AND ((pow(2,`log2(fold_change)`) >= %s OR pow(2,`log2(fold_change)`) <= %s OR `log2(fold_change)` = 'inf' OR `log2(fold_change)` = '-inf') AND `%s` <= %s)"%(FC_input,FC_input,TEST_column,TEST_input)
    else:
        sql += " WHERE (ROUND(`value_1`,3) >= 0.001 AND ROUND(`value_2`,3) >= 0.001) "


    print(sql)
    cursor.execute(sql)
    table_data = list(map(list,cursor.fetchall()))
    download_table_data = []
    
    if table_data:
        if DE_level =='genes':
            if TEST_input:
                for row in table_data:
                    if float(row[-2])  >= 1 :
                        row[-2] = 1
                    else:
                        row[-2] = float(row[-2]) 
                    if (row[-3] == 'inf' or row[-3] == '-inf'):
                        row[-3] = float(row[-3])
                    else :
                        row[-3] = row[-1]
                    row[-4] = float(row[-4])
                    row[-5] = float(row[-5])
                    row.pop()
                    
                    download_table_data.append([row[0],round(row[2],3),round(row[1],3),round(round(row[2],3)/round(row[1],3),3),row[4]])
            else:
                for row in table_data:
                    if (row[-2] == 'inf' or row[-2] == '-inf'):
                        row[-2] = float(row[-2])
                    else :
                        row[-2] = row[-1]
                    row[-3] = float(row[-3])
                    row[-4] = float(row[-4])
                    row.pop()
                    
                    download_table_data.append([row[0],round(row[2],3),round(row[1],3),round(round(row[2],3)/round(row[1],3),3)])
                    
        else:
            path = '/home/edward/Django/edward_project/static/error_isoforms.txt'
            f = open(path, 'r')
            error_isoforms = eval(f.read())
            f.close()
            table_data = pd.DataFrame(table_data)
            table_data = table_data[~table_data[0].isin(error_isoforms)]
            table_data = table_data.values.tolist()
            if TEST_input:
                for row in table_data:
                    if float(row[-2])  >= 1 :
                        row[-2] = 1
                    else:
                        row[-2] = float(row[-2]) 
                    if (row[-3] == 'inf' or row[-3] == '-inf'):
                        row[-3] = float(row[-3])
                    else :
                        row[-3] = row[-1]
                    row[-4] = float(row[-4])
                    row[-5] = float(row[-5])
                    row.pop()
                    download_table_data.append([row[0],row[1],round(row[3],3),round(row[2],3),round(round(row[3],3)/round(row[2],3),3)])
            else:
                for row in table_data:
                    if (row[-2] == 'inf' or row[-2] == '-inf'):
                        row[-2] = float(row[-2])
                    else :
                        row[-2] = row[-1]
                    row[-3] = float(row[-3])
                    row[-4] = float(row[-4])
                    row.pop()
                    
                    download_table_data.append([row[0],row[1],round(row[3],3),round(row[2],3),round(round(row[3],3)/round(row[2],3),3)])
                    
        cursor.close()
        if TEST_input:
            table_data = sorted(table_data,key=operator.itemgetter(4, 3))
        else:
            table_data = sorted(table_data,key=operator.itemgetter(3, 2), reverse=True)
        end = time.time()
    
    return table_data,download_table_data