import sqlite3
import os
# from tabulate import tabulate 
def search(gene_name):
    #connect to database
    conn = sqlite3.connect('Hw/database.db')
    cursor = conn.cursor()

    #find transripts
    query = "SELECT field2 FROM hg38_gene_transcripts_20180130 WHERE field1 = ?"
    cursor.execute(query,(gene_name,))
    result_1 = cursor.fetchall()
    result_1 = result_1[0][0]
    transcripts = result_1.split(',')

    #find gene info
    query = "SELECT field6,field8,field10,field15 FROM NCBI_gene_info_20180201 WHERE field4 = ?"
    cursor.execute(query,(gene_name,))
    result_2 = cursor.fetchall()

    #find transcript info
    trans_info_list = []
    for transcript in transcripts:
        query = "SELECT field2,field3 FROM NCBI_transcript_info_20180209 WHERE field1 = ?"
        cursor.execute(query,(transcript,))
        result_3 = cursor.fetchall()
        trans_info = {
        "transcript_name":transcript,
        "definition":result_3[0][0],
        "transcript_variant":result_3[0][1]
        }
        trans_info_list.append(trans_info)
        
    data_class = [{
        "Symbol":gene_name,
        "Synonyms":result_2[0][0],
        "Chromosome":result_2[0][1],
        "Description":result_2[0][2],
        "Other Designations":result_2[0][3]
    }]
    
    return data_class,trans_info_list

def search_bonus(chromosome_number , type_of_gene):
    #connect to database
    conn = sqlite3.connect('../database/db')
    cursor = conn.cursor()
    
    #find gene info
    data_class_list = []
    query = """
            SELECT Symbol FROM NCBI_gene_info_20180201 
            WHERE chromosome IN ({}) 
            AND type_of_gene IN ({})
            """.format(','.join(['?'] * len(chromosome_number)), ','.join(['?'] * len(type_of_gene)))

            
    cursor.execute(query,chromosome_number+type_of_gene)
    result_bonus = cursor.fetchall()
    for row in result_bonus:
        data_class = {
            "Symbol": row[0]
        }
        data_class_list.append(data_class)
    return data_class_list

if __name__ == "__main__":
    input_1 = input("Enter : ")
    input_2 = input("Enter : ")
    
    number_list = input_1.split(",")
    type_list = input_2.split(",")
    data_class_list = search_bonus(number_list,type_list)
    
    # print(tabulate(data_class_list, headers="keys", tablefmt="plain"))
    print(data_class_list)
    