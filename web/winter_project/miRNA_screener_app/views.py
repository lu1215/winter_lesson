from django.shortcuts import render
import pandas as pd
from django.db import connections
import csv
import os
from miRNA_screener_app.miRNA_screener import miRNAscreener_getdata
from django.http import JsonResponse
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
    elif set_select == "intersection":
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
    df_miRNA.rename(columns={"gene_name": "name"},inplace=True)
    print(df_miRNA)
    return JsonResponse({"result": df_miRNA.to_dict('records')})