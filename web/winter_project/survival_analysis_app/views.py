from django.shortcuts import render
from django.http import HttpResponse, JsonResponse
from survival_analysis_app.survival_analysis_v3 import *
from django.views.decorators.csrf import csrf_exempt, csrf_protect
from survival_analysis_app.class_list import *
from django.template.defaulttags import register


# this function can be called in html
@register.filter
def get_item(dictionary, key):
    return dictionary.get(key)

@register.filter(name='divide')
def divide(value,arg):
    try:
        return float(value) / float(arg)
    except (ValueError, ZeroDivisionError):
        return None
# Create your views here.
def survival_analysis_page(request):
    return render(request, "survival_analysis.html", locals())

def sur_filtering(request):
    input_type = request.POST["type"]
    cancer = request.POST["cancer"]
    stage = request.POST["stage"]
    high_percent = request.POST["high_percent"]
    low_percent = request.POST["low_percent"]
    input_pvalue = request.POST["pvalue"]
    primary_site, project = cancer.split("|")
    table_name = '%s_%s_FPKM_Cufflinks'%(project,input_type)
    column_table = "%s|%s"%(stage,table_name)
    # df_pvalue, _ = cal_pvalue_main(request)
    print(input_type, cancer, stage, high_percent, low_percent, input_pvalue)
    _, result_list = cal_pvalue_main(input_type, cancer, stage, high_percent, low_percent, input_pvalue)
    return JsonResponse({"result":result_list})

@csrf_exempt
def detail_page(request):
    stage_dict = {
        'normal' :'normal',
        'stage i' : 'stage_1',
        'stage ii' : 'stage_2',
        'stage iii' : 'stage_3',
        'stage iv' : 'stage_4',
    }
    primary_site, project = request.GET["cancer"].split("|")
    input_type = request.GET["type"]
    name = request.GET["name"]
    Low_Percentile = request.GET["lp"]
    High_Percentile = request.GET["hp"]
    stage = request.GET["stage"]
    max_days = survival_max_days(project, name, input_type, stage)["max_survival_days"]
    img_str, p_value= survival_plot_realtime(project, primary_site, input_type, name ,'0', Low_Percentile, High_Percentile, stage)
    gene = Detail.belong_gene(input_type,name) if input_type == 'isoforms' else name
    genelink_dict = Filter.gene_link([gene],input_type)
    all_transcript1 = Detail.all_transcript(gene)
    all_transcript = []
    path = f'{current_path}/../static/data/error_isoforms.txt'
    f = open(path, 'r')
    error_isoforms = eval(f.read())
    f.close()
    for i in all_transcript1:
        if(i not in error_isoforms):
            all_transcript.append(i)
    NCBI_gene_summary = Summary.NCBI_gene_summary(gene)
    NCBI_transcript_summary = Summary.NCBI_transcript_summary(all_transcript)
    primary_key = input_type[:-1]+'_name'
    stage_dict_gather = {
        'normal' :'normal',
        'stage i' : 'stage_1',
        'stage ii' : 'stage_2',
        'stage iii' : 'stage_3',
        'stage iv' : 'stage_4',
    }
    stage_Arabic_dict = {
        'normal' :'N',
        'stage i' : '1',
        'stage ii' : '2',
        'stage iii' : '3',
        'stage iv' : '4',
    }
    primary_stage = OrderedDict()
    primary_condition = OrderedDict()
    project_list = []
    primary_condition = OrderedDict(sorted(primary_condition.items()))
    primary_stage = OrderedDict(sorted(primary_stage.items()))
    condition_list = list(primary_condition.values())
    group_name = list(primary_stage.values())
    return render(request, 'survival_analysis_detail.html', locals())