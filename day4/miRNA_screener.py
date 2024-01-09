import sqlite3
import os
import pandas as pd

current_path = os.path.dirname(__file__)
def miRNAscreener_getdata(miRNA_list: list, set_operation: str) -> list:
    db_path = f"{current_path}/../database/db"
    with sqlite3.connect(db_path, check_same_thread=False) as db_conn:
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
        # df.drop_duplicates(inplace=True)
        counts = df['gene_name'].value_counts()
        df = df[df['gene_name'].isin(counts[counts == len(miRNA_list)].index)]
        df = df[['gene_name']]
        df.drop_duplicates(inplace=True)
        for e in miRNA_list:
            df[e] = "o"
        result_list = df.to_dict('records')
    elif set_operation == "UNION":
        df = pd.DataFrame(result_list)
        df['value'] = 1
        df = df.pivot_table(index='gene_name', columns='mirna_name', values='value', fill_value=0).reset_index()
        df = df.sort_values(by='gene_name').reset_index(drop=True)
        df.iloc[:, 1:] = df.iloc[:, 1:].applymap(lambda x: 'o' if x == 1 else 'x')
        result_list = df.to_dict('records')
    elif set_operation == "DIFFERENCE":
        df = pd.DataFrame(result_list)
        df.drop_duplicates(subset=['gene_name'], keep=False, inplace=True)
        df = df[df["mirna_name"] == miRNA_list[0]]
        df = df[['gene_name']]
        df[miRNA_list[0]] = "o"
        df[miRNA_list[1]] = "x"
        result_list = df.to_dict('records')
    return result_list

if __name__ == "__main__":
    miRNA_list = ["hsa-miR-34a-3p", "hsa-miR-21-5p"]
    set_operation = "INTERSECT"
    data = miRNAscreener_getdata(miRNA_list, set_operation)
    print(pd.DataFrame(data))