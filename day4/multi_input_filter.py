import sqlite3
import os
import pandas as pd


current_path = os.path.dirname(__file__)
def multi_input_getdata(*input_dict):
    # example input {"miRNA_list": ["hsa-miR-21-5p", "hsa-miR-34a-5p"]}, {"chromosome_number": [6]} ...
    db_path = f"{current_path}/../database/db"
    with sqlite3.connect(db_path, check_same_thread=False) as db_conn:
        cursor = db_conn.cursor()
        sql_command = f'SELECT mirna_name,gene_name FROM Homo_sapiens_miRNA WHERE '
        for count, element in enumerate(input_dict):
            print(element)
            if count != 0: sql_command += " AND "
            key = list(element.keys())[0]
            values = list(element.values())[0]
            print(len(values))
            for idx in range(len(values)):
                value = values[idx]
                print(key, value)
                if len(values) == 1:
                    sql_command += f'( {key} = "{value}") '
                elif idx == 0 :
                    sql_command += f'( {key} = "{value}" OR '
                elif idx == len(values) - 1:
                    sql_command += f'{key} = "{value}")'
                else:
                    sql_command += f" {key} = '{value}' OR "
        print(sql_command)
        cursor.execute(sql_command)
        columns = [col[0] for col in cursor.description]
        result = [dict(zip(columns, row)) for row in cursor.fetchall()]
        return result
        
if __name__ == "__main__":
    data = multi_input_getdata({"miRNA_name":["hsa-miR-196a-5p"]}, {"gene_name":["CYP26B1", "USP28"]})
    # data = multi_input_getdata({"miRNA_name":["hsa-miR-34a-5p", "hsa-miR-21-5p"]}, {"chromosome": ["chr2"]})
    print(pd.DataFrame(data))
    # {"field1": "mirna_name"}