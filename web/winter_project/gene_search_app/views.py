from django.shortcuts import render
from gene_search_app.gene_info import gene_info
from django.http import HttpResponse, JsonResponse
import pandas as pd
# from django import 


# Create your views here.
def search_page(request):
    genes = pd.read_csv('static/data/gene_names/all_genes.csv')["gene"].tolist()
    return render(request, "gene_search.html", locals())

def searching(request):
    # print(request.POST)
    search_val = request.POST["search_val"]
    result = gene_info(search_val)
    print(result.get_gene_info())
    print(result.get_transcript_info())
    return JsonResponse({"gene_info":[result.get_gene_info()], "transcript_info":result.get_transcript_info()})
