import sqlite3
import os
import itertools 


def advanced_filter(chromosome_number:list, type_of_gene:list):
    current_path = os.path.dirname(__file__)
    target = list(itertools.product(chromosome_number, type_of_gene))
    
    # this need to adjust when upload to server
    db_path = f"{current_path}/../database/db"
    # db_path = f"{current_path}/db.sqlite3"
    with sqlite3.connect(db_path, check_same_thread=False) as db_conn:
        cursor = db_conn.cursor()
        query = f"SELECT `Symbol` FROM `NCBI_gene_info_20180201` WHERE "
        for idx, e in enumerate(target):
                if idx == len(target) - 1:
                    query += f"(`chromosome` = '{e[0]}' AND `type_of_gene` = '{e[1]}')"
                else:
                    query += f"(`chromosome` = '{e[0]}' AND `type_of_gene` = '{e[1]}') OR"
        print(query)
        cursor.execute(query)
        result = cursor.fetchall()
    return result

if __name__ == "__main__":
    print(advanced_filter([22, 9], ['ncRNA']))