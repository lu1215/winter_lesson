import sqlite3
import os
import pandas as pd

current_path = os.path.dirname(__file__)
def miRNAscreener_getdata(miRNA_list: list, set_operation: str) -> list:
    with sqlite3.connect(db_path, check_same_thread=False) as db_conn:
        db_path = f"{current_path}/../database/db"
        cursor = db_conn.cursor()
        sql_command = ""
        if len(miRNA_list) <= 1:
            sql_command = f'SELECT DISTINCT mirna_name,gene_name FROM Homo_sapiens_miRNA WHERE mirna_name = "{miRNA_list[0]}";'
        else:
            for i in range(len(miRNA_list)):
                if i != 0:
                    sql_command += f" UNION "
                    # sql_command += f"OR "
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
        # df.to_csv("test.csv", index=False)
        # 使用 pivot_table 將資料表進行轉換
        df = df.pivot_table(index='gene_name', columns='mirna_name', values='value', fill_value=0).reset_index()
        # 重新排序欄位
        df = df.sort_values(by='gene_name').reset_index(drop=True)
        # 將 0 與 1 轉換成 o 與 x
        df.iloc[:, 1:] = df.iloc[:, 1:].applymap(lambda x: 'o' if x == 1 else 'x')
        result_list = df.to_dict('records')
        # result.to_csv("test.csv", index=False)
    elif set_operation == "DIFFERENCE":
        df = pd.DataFrame(result_list)
        df.drop_duplicates(subset=['gene_name'], keep=False, inplace=True)
        df = df[df["mirna_name"] == miRNA_list[0]]
        # df.to_csv("test.csv", index=False)
        df = df[['gene_name']]
        result_list = df.to_dict('records')
    return result_list