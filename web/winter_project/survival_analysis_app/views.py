from django.shortcuts import render
from django.http import HttpResponse, JsonResponse
from survival_analysis_app.survival_analysis import *


# Create your views here.
def survival_analysis_page(request):
    return render(request, "survival_analysis.html", locals())

def sur_filtering(request):
    # search_type = request.POST["type"]
    # cancer = request.POST["cancer"]
    # stage = request.POST["stage"]
    # high_percent = request.POST["high_percent"]
    # low_percent = request.POST["low_percent"]
    # pvalue = request.POST["pvalue"]
    # primary_site, project = cancer.split("|")
    
    # p_value = survival_plot_realtime(search_type, cancer, stage, high_percent, low_percent, pvalue)
    
    input_type = request.POST["type"]
    cancer = request.POST["cancer"]
    stage = request.POST["stage"]
    high_percent = request.POST["high_percent"]
    low_percent = request.POST["low_percent"]
    input_pvalue = request.POST["pvalue"]
    primary_site, project = cancer.split("|")

    table_name = '%s_%s_FPKM_Cufflinks'%(project,input_type)
    column_table = "%s|%s"%(stage,table_name)
    cursor = connections['edward_Cufflinks'].cursor()
    primary_key = 'gene_name' if input_type == 'genes' else 'isoform_name'
    stage_list = ['stage_1','stage_2','stage_3','stage_4']
    all_cancer_data = get_allcancer_data(column_table, input_type)
    result_list = []