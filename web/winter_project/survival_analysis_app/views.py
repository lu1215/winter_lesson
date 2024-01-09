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
    img_str= survival_plot_realtime(project, primary_site, input_type, name ,'0', Low_Percentile, High_Percentile, stage)
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
    # project_primary_dict = Filter.project_primary_dict()[0]
    # primary_list_all = list(project_primary_dict.values())
    # number_of_sample = Filter.project_primary_dict()[1]
    primary_key = input_type[:-1]+'_name'
    # TEST_column = "%s_%s"%(TEST_select.replace(' ','_'),TESTstates_select.replace(' ','').lower()) if TEST_select else ''

    # condition_FC_qv_list = Filter.filter(input_type,name,primary_key,list(project_primary_dict.keys()),FC_select,FC_input,TEST_column,TEST_input)
    # primary_list = sorted(list(set([project_primary_dict[x[0]] for x in condition_FC_qv_list])))
        
    # for cancer in primary_list:
    #     condition_FC_qv_list_cancer = []
    #     for row in condition_FC_qv_list:
    #         if project_primary_dict[row[0]] == cancer:
    #             condition_FC_qv_list_cancer += [row]
        # download_table = Filter.download_table(cancer,condition_FC_qv_list_cancer,name,FC_select,FC_input,TEST_select,TESTstates_select,TEST_input)
        # request.session["DE_Conditions_%s_%s"%(cancer,"0")] = download_table
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
    # primary_list = []
    # for row in condition_FC_qv_list:
    #     primary_project = "%s|%s"%(project_primary_dict[row[0]],row[0])
    #     c1_pair = "%s(%s_%s)"%(stage_dict_gather[row[1]],stage_Arabic_dict[row[1]],stage_Arabic_dict[row[2]])
    #     c2_pair = "%s(%s_%s)"%(stage_dict_gather[row[2]],stage_Arabic_dict[row[1]],stage_Arabic_dict[row[2]])
    #     if primary_project not in primary_stage:
    #         primary_stage[primary_project] = [[c1_pair,c2_pair]]
    #         primary_condition[primary_project] = [row[1],row[2]]
    #     else:
    #         primary_stage[primary_project] += [[c1_pair,c2_pair]]
    #         primary_condition[primary_project] += [row[1],row[2]]
    #     if row[0] not in project_list:
    #         project_list += [row[0]]
    primary_condition = OrderedDict(sorted(primary_condition.items()))
    primary_stage = OrderedDict(sorted(primary_stage.items()))

    condition_list = list(primary_condition.values())
    group_name = list(primary_stage.values())
    # condition_list = list(map(lambda x:sorted(x),list(primary_stage.values())))
    # boxplot_data,primary_nodata = Boxplot.boxplot_data(primary_stage,input_type,name)
    # print(boxplot_data)
    return render(request, 'survival_analysis_detail.html', locals())