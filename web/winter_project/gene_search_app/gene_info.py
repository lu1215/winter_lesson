import sqlite3
import os

class gene_info:
    def __init__(self, gene_name:str):
        self.gene_name = gene_name
        self.gene_trans_db_name = "hg38_gene_transcripts_20180130"
        self.Synonyms_db_name = "NCBI_gene_info_20180201"
        self.trans_db_name = "NCBI_transcript_info_20180209"
        self.mRNA_transcript = self.query_db(self.gene_trans_db_name, ["`transcripts`"], ["`gene`"], self.gene_name)
        self.Synonyms, self.Chromosome, self.Description, self.Other_Designations = self.query_db(self.Synonyms_db_name, ["`Synonyms`","`chromosome`","`description`","`Other_designations`"], ["`Symbol`"], self.gene_name)
        self.transcript_info = self.query_db_for_transcript(self.trans_db_name, "*", self.mRNA_transcript[0].split(','))

    def query_db(self, db_table_name:str, db_column_name:list, targets:list, db_column_value:str, to_dict_list:bool=False):
        current_path = os.path.dirname(__file__)
        # this need to adjust when upload to server
        db_path = f"{current_path}/../../../database/db"
        print(db_path)
        # db_path = f"{current_path}/db.sqlite3"
        db_column_name = ','.join(db_column_name)
        targets = ','.join(targets)
        # https://nkust.gitbook.io/python/sqlite-liao-cao-zuo-jie
        with sqlite3.connect(db_path, check_same_thread=False) as db_conn:
            cursor = db_conn.cursor()
            query = f"SELECT {db_column_name} FROM `{db_table_name}` WHERE {targets} = '{db_column_value}'"
            print(query)
            cursor.execute(query)
            # if result:
                # self.mRNA_transcripts = result
            if to_dict_list:
                columns = [col[0] for col in cursor.description]
                result = [dict(zip(columns, row)) for row in cursor.fetchall()]
            else:
                result = cursor.fetchall()[0]
        return result
    
    def query_db_for_transcript(self, db_table_name:str, db_column_name:list = "*", transcript_list:list = []):
        current_path = os.path.dirname(__file__)
        # this need to adjust when upload to server
        db_path = f"{current_path}/../../../database/db"
        # db_path = f"{current_path}/db.sqlite3"
        db_column_name = ','.join(db_column_name)
        with sqlite3.connect(db_path, check_same_thread=False) as db_conn:
            cursor = db_conn.cursor()
            query = f"SELECT {db_column_name} FROM `{db_table_name}` WHERE "
            for idx, transcript in enumerate(transcript_list):
                if idx == len(transcript_list)-1:
                    query += f"`transcript_name` = '{transcript}'"
                else:
                    query += f"`transcript_name` = '{transcript}' OR "
            print(query)
            cursor.execute(query)
            columns = [col[0] for col in cursor.description]
            result = [dict(zip(columns, row)) for row in cursor.fetchall()]
        return result
    
    def get_gene_info(self):
        return {
            "Symbol":self.gene_name,
            "mRNA transcript(s)": ','.join(self.mRNA_transcript),
            "Synonyms": self.Synonyms,
            "Chromosome": self.Chromosome,
            "Description": self.Description,
            "Other Designations": self.Other_Designations
        }
    
    def get_transcript_info(self):
        return self.transcript_info

if __name__ == "__main__":
    data = gene_info("BRCA1")
    print(data.get_gene_info())
    print(data.get_transcript_info())