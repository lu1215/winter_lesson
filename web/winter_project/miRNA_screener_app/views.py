from django.shortcuts import render
import pandas as pd
from django.db import connections
import csv
import os
current_path = os.path.dirname(__file__)


# Create your views here.
def miRNA_page(request):
    with open(f'{current_path}/../static/data/miRNA/homo_mirna_list.csv', 'r') as file:
        reader = csv.reader(file)
        homo_miRNA_list = list(reader)
    return render(request, 'miRNA_screener.html', locals())

def miRNA_main(request):
    selected_miRNA = request.POST.getlist("selected_miRNA[]")
    set_select = request.POST["miRNA_set"]
    if set_select == "union":
        set_operation = "UNION"
    elif set_select == "intersect":
        set_operation = "INTERSECT"
    elif "difference" in set_select :
        # Set the value for the "difference" case
        if set_select == "difference2_1":
            selected_miRNA = [selected_miRNA[1], selected_miRNA[0]]
        set_operation = "DIFFERENCE"
    ## query mysql
    df_miRNA = pd.DataFrame(miRNAscreener_getdata(selected_miRNA, set_operation))
    df_miRNA.dropna(subset = ['gene_name'], inplace=True)
    if set_operation=="INTERSECT":
        for e in selected_miRNA:
            df_miRNA[e] = "o"
    elif set_operation=="DIFFERENCE":
        df_miRNA[selected_miRNA[0]] = "o"
        df_miRNA[selected_miRNA[1]] = "x"
    df = df_miRNA if len(df) == 0 else pd.merge(df, df_miRNA, left_on="name", right_on="gene_name", how="inner")

def miRNAscreener_getdata(miRNA_list: list, set_operation: str) -> list:
    cursor = connections['edward_miRNA'].cursor()
    sql_command = ""
    if len(miRNA_list) <= 1:
        sql_command = f'SELECT DISTINCT mirna_name,gene_name FROM Homo_sapiens_miRNA WHERE mirna_name = "{miRNA_list[0]}";'
    else:
        for i in range(len(miRNA_list)):
            if i != 0:
                sql_command += f" UNION "
            sql_command += f'SELECT mirna_name,gene_name FROM Homo_sapiens_miRNA WHERE mirna_name = "{miRNA_list[i]}"'
    print(sql_command)
    cursor.execute(sql_command)
    columns = [col[0] for col in cursor.description]
    rows = cursor.fetchall()
    result_list = [dict(zip(columns, row)) for row in rows]
    if set_operation == "INTERSECT":
        df = pd.DataFrame(result_list)
        counts = df['gene_name'].value_counts()
        df = df[df['gene_name'].isin(counts[counts == len(miRNA_list)].index)]
        df.to_csv("test.csv", index=False)
        df = df[['gene_name']]
        df.drop_duplicates(inplace=True)
        result_list = df.to_dict('records')
    elif set_operation == "UNION":
        df = pd.DataFrame(result_list)
        df['value'] = 1
        df = df.pivot_table(index='gene_name', columns='mirna_name', values='value', fill_value=0).reset_index()
        # 重新排序欄位
        df = df.sort_values(by='gene_name').reset_index(drop=True)
        # 將 0 與 1 轉換成 o 與 x
        df.iloc[:, 1:] = df.iloc[:, 1:].applymap(lambda x: 'o' if x == 1 else 'x')
        result_list = df.to_dict('records')
    elif set_operation == "DIFFERENCE":
        df = pd.DataFrame(result_list)
        df.drop_duplicates(subset=['gene_name'], keep=False, inplace=True)
        df = df[df["mirna_name"] == miRNA_list[0]]
        df = df[['gene_name']]
        result_list = df.to_dict('records')
    cursor.close()
    return result_list