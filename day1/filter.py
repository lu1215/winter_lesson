import sqlite3
import os

def advanced_filter(chromosome_number:int, type_of_gene:str):
    current_path = os.path.dirname(__file__)
    # this need to adjust when upload to server
    db_path = f"{current_path}/../database/db"
    # db_path = f"{current_path}/db.sqlite3"
    with sqlite3.connect(db_path, check_same_thread=False) as db_conn:
        cursor = db_conn.cursor()
        query = f"SELECT `Symbol` FROM `NCBI_gene_info_20180201` WHERE `chromosome` = '{chromosome_number}' AND `type_of_gene` = '{type_of_gene}'"
        print(query)
        cursor.execute(query)
        result = cursor.fetchall()
    return result

if __name__ == "__main__":
    print(advanced_filter(1, 'protein-coding'))