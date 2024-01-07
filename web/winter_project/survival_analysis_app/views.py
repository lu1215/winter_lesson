from django.shortcuts import render
from django.http import HttpResponse, JsonResponse
from survival_analysis_app.survival_analysis_v2 import *


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
    _, result_list = cal_pvalue_main(input_type, cancer, stage, high_percent, low_percent, input_pvalue)
    return JsonResponse({"result":result_list})

def get_allcancer_data(column_table,search_by):
	stage_dict = {
		'stage i' : 'stage_1',
		'stage ii' : 'stage_2',
		'stage iii' : 'stage_3',
		'stage iv' : 'stage_4',
	}
	current_path = os.path.dirname(__file__)
	db_path = f"{current_path}/../../../database/db"
	# db_path = f"{current_path}/db.sqlite3"
	cancer = "liver_cancer|TCGA-LIHC"
	table_name = 'TCGA-LIHC_genes_FPKM_Cufflinks'
	# column_table = "%s|%s"%(stage,table_name)
	# db_column_name = ','.join(db_column_name)
	with sqlite3.connect(db_path, check_same_thread=False) as db_conn:
		cursor = db_conn.cursor()
		# cursor = connections['edward_Cufflinks'].cursor()
		primary_key = 'gene_name' if search_by == 'genes' else 'isoform_name'
		# column = primary_key
		# column = stage_dict[column_table.split('|')[0]] if column_table.split('|')[0] != 'all stage' else ','.join(stage_dict.values())
		column = f"{primary_key}, {stage_dict[column_table.split('|')[0]] if column_table.split('|')[0] != 'all stage' else ','.join(stage_dict.values())}"
		# table_name = column_table.split('|')[1]
		cursor.execute("SELECT %s FROM `%s`"%(column,table_name))
		result = list(cursor.fetchall())
	return result