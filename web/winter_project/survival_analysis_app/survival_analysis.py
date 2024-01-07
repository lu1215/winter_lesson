import numpy as np
from argparse import ArgumentParser
import os
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import csv
import sqlite3
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import base64
import io
matplotlib.use('Agg')


file_path = os.path.dirname(os.path.abspath(__file__))
class Survival_plot():
	def survival_data_realtime(project, column_table, search_by, GT_input):
		stage_dict = {
			'stage i' : 'stage_1',
			'stage ii' : 'stage_2',
			'stage iii' : 'stage_3',
			'stage iv' : 'stage_4',
		}
		primary_key = 'gene_name' if search_by == 'genes' else 'isoform_name'
		column = f"{primary_key}, {stage_dict[column_table.split('|')[0]] if column_table.split('|')[0] != 'all stage' else ','.join(stage_dict.values())}"
		table_name = column_table.split('|')[1]
		db_path = f"{file_path}/../database/db"
		with sqlite3.connect(db_path, check_same_thread=False) as db_conn:
			cursor = db_conn.cursor()
			query = f"SELECT {column} FROM `{table_name}`"
			cursor.execute(query)
			columns = [col[0] for col in cursor.description]
			result = cursor.fetchall()
		# df = pd.read_csv(f"{file_path}/data/{project}_{search_by}_FPKM_Cufflinks.csv")
		df = pd.DataFrame(result, columns=columns)
		table_name = column_table.split('|')[1]
		# primary_key = 'gene_name' if search_by == 'genes' else 'isoform_name'
		selected_columns = [primary_key]
		if column_table.split('|')[0] != 'all stage':
			print(column_table.split('|')[0])
			# selected_columns += stage_dict[column_table.split('|')[0]] if column_table.split('|')[0] != 'all stage' else ','.join(stage_dict.values())
			selected_columns += [stage_dict[column_table.split('|')[0]]]
		else:
			data_columns = df.columns.tolist()
			data_columns = data_columns[1:]
			selected_columns = [primary_key]
			for e in data_columns:
				if e != 'normal' and e !='all_stage' and e != 'pvalue(normal_allstage)' and e != 'pvalue(allstage_normal)':
					selected_columns.append(e)
			# selected_columns = ','.join(selected_columns)
		df = df[selected_columns]
		df_result = df[df[f"{primary_key}"] == GT_input]
		result = df_result.drop(columns= [f"{primary_key}"]).values.tolist()[0]
		return result

	def survival_plot(T1,E1,T2,E2,GT_input,primary_site,random_id,Low_Percentile,High_Percentile,survival_days,survival_select):
		# from lifelines.estimation import KaplanMeierFitter
		plot_path = f"{file_path}/plot_result/survival_plot_{primary_site.replace(' ','_').replace('(','').replace(')','')}_{GT_input}_{int(survival_days)}_{survival_select}.png"
		try:
			os.remove(plot_path)
		except OSError:
			pass
		kmf = KaplanMeierFitter()
		dpi = 100
		logrank_result = logrank_test(T1, T2, E1, E2)
		logrank_p_value = logrank_result.p_value
		logrank_test_statistic = logrank_result.test_statistic
		if T2 != [] and E2 != []:
			kmf.fit(T2, event_observed=E2, label='Low Expression (n={}, {}%)'.format(len(T2),Low_Percentile))
			ax = kmf.plot(ci_show=False,color='green',show_censors=True, figsize=(1200/dpi, 800/dpi))
		if T1 != [] and E1 != []:
			kmf.fit(T1, event_observed=E1, label='High Expression (n={}, {}%)'.format(len(T1),High_Percentile))
			kmf.plot(ax=ax,ci_show=False,color='red',show_censors=True) 
		font = {'family' : 'verdana'}
		matplotlib.rc('font', **font)
		plt.subplots_adjust(left=0.06, right=0.94, top=0.94, bottom=0.06)
		plt.title("%s"%GT_input)
		ax.text(0.1, 0.2, 'Low Percentile: {}%\nHigh Percentile: {}%\nDays: {}\nCondition: {}'.format(Low_Percentile,High_Percentile,survival_days,survival_select), horizontalalignment='left',verticalalignment='center', transform=ax.transAxes, bbox={'facecolor':'#F5F5F5', 'pad':5},fontsize=10)
		ax.text(0.1, 0.1, 'logrank pvalue: {}'.format(round(logrank_p_value,3)), horizontalalignment='left',verticalalignment='center', transform=ax.transAxes, bbox={'facecolor':'#F5F5F5', 'pad':5},fontsize=10)
		ax.set_ylim(ymin=0)
		ax.set_xlabel('Days',fontsize=12)
		ax.set_ylabel('Survival probability',fontsize=12)
		ax.grid(True)
		gridlines = ax.get_xgridlines() + ax.get_ygridlines()
		for line in gridlines:
			line.set_linestyle('-.')
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
		plt.show()
		plt.savefig(plot_path,dpi=dpi)
		img_str = io.BytesIO()
		plt.savefig(img_str,dpi=300,format = 'png')
		plt.clf()
		plt.close()
		img_str = base64.b64encode(img_str.getvalue()).decode("utf-8").replace("\n", "")
		plt.close()
		return logrank_p_value, img_str

	def survival_download(T1,E1,T2,E2,high_case,low_case,high_FPKM,low_FPKM,GT_input,primary_site,random_id,Low_Percentile,High_Percentile,survival_select):
		output_data = [
						['Query: %s'%(GT_input)],
						['Primary site: %s'%(primary_site)],
						['Low Percentile: %s'%(Low_Percentile)],
						['High Percentile: %s'%(High_Percentile)],
						['Condition: %s'%(survival_select)],
						[],
						['Patient','Days','Status','Expression','Group']
					]
		for idx,row in enumerate(low_case):
			Status = 'Alive' if E2[idx] == False else 'Dead'
			output_data += [[row,T2[idx],Status,low_FPKM[idx],'Low']]
		for idx,row in enumerate(high_case):
			Status = 'Alive' if E1[idx] == False else 'Dead'
			output_data += [[row,T1[idx],Status,high_FPKM[idx],'High']]
		print('survival_download')

		# download file
		filename=f"Survival_Profile_{primary_site}_{Low_Percentile}_{High_Percentile}.csv"
		with open(f"{file_path}/csv_result/{filename}", "w") as f:
			writer = csv.writer(f)
			for e in output_data:
				writer.writerow(e)
		
		return output_data

def survival_plot_realtime(request:dict):
	project = request['project']
	# ex: TCGA-COAD
	primary_site = request['primary_site']
	# ex: Adrenal Gland Adrenocortical Carcinoma
	search_by = request['search_by']
	GT_input = request['GT_input']
	random_id = request['random_id']
	Low_Percentile = request['Low_Percentile']
	High_Percentile = request['High_Percentile']
	survival_days = request['survival_days']
	survival_select = request['survival_select']
	# stage : all stage, stage i, stage ii, stage iii, stage iv
	table_name = '%s_%s_FPKM_Cufflinks'%(project,search_by)
	column_table = "%s|%s"%(survival_select,table_name)
 #### patched by t50504
	survival_data = Survival_plot.survival_data_realtime(project, column_table,search_by,GT_input)
	survival_str = ""
	case_id_list = []
	FPKM_list = [float(y.split("|")[0]) for x in survival_data for y in x.split(',')]
	low_quartile = np.percentile(FPKM_list, float(Low_Percentile))
	high_quartile = np.percentile(FPKM_list, 100-float(High_Percentile))
	T1 = [] #high 存活天數
	E1 = [] #high 是否死亡
	T2 = []
	E2 = []
	high_case = []
	low_case = []
	high_FPKM = []
	low_FPKM = []
	for stage in survival_data:
		for info in stage.split(','):
			FPKM = float(info.split('|')[0])
			case_id = info.split('|')[1]
			survival_times = float(info.split('|')[2]) if info.split('|')[2] != 'None' else info.split('|')[2] #存活天數
			# print(case_id,survival_times)
			survival_events = False if info.split('|')[3] == 'alive' else True #是否死亡
			if FPKM > high_quartile and (survival_times != 0 and survival_times != 'None') and survival_times <= float(survival_days):
				T1 += [survival_times]
				E1 += [survival_events]
				case_id_list += [case_id]
				high_case += [case_id]
				high_FPKM += [FPKM]
			elif FPKM < low_quartile and (survival_times != 0 and survival_times != 'None') and survival_times <= float(survival_days):
				T2 += [survival_times]
				E2 += [survival_events]
				case_id_list += [case_id]
				low_case += [case_id]
				low_FPKM += [FPKM]
	if (T2 != [] and E2 != []) and (T1 != [] and E1 != []):
		pvalue, img_str = Survival_plot.survival_plot(T1,E1,T2,E2,GT_input,primary_site,random_id,Low_Percentile,High_Percentile,max(T1+T2)+1,survival_select)
		survival_download = Survival_plot.survival_download(T1,E1,T2,E2,high_case,low_case,high_FPKM,low_FPKM,GT_input,primary_site,random_id,Low_Percentile,High_Percentile,survival_select)
	else:
		pvalue = 1
		survival_str = ' Survival analysis is not available for '+GT_input+' since more than half of the samples have zero expression.'
	return pvalue, img_str

def survival_max_days(project, GT_input, search_by, survival_select) -> float:
	table_name = '%s_%s_FPKM_Cufflinks'%(project,search_by)
	column_table = "%s|%s"%(survival_select,table_name)
	survival_data = "".join(Survival_plot.survival_data_realtime(project, column_table,search_by,GT_input)).split(",")
	survival_days = [float(x.split("|")[2]) for x in survival_data if x.split("|")[2] != 'None']
	max_survival_days = max(survival_days)
	return max_survival_days

def sur_main(input_name, Low_Percentile, High_Percentile, survival_select):
	## cancer ➔ project table
	# ex : python survival_analysis.py -p TCGA-ACC --primary_site Adrenal_Gland_Adrenocortical_Carcinoma -t genes -n KIF23 --Low_Percentile 50 --High_Percentile 50 --survival_days 4628 --survival_select all_stage
	# ex : python survival_analysis.py -p TCGA-ACC --primary_site Adrenal_Gland_Adrenocortical_Carcinoma -t isoforms -n NM_000014 --Low_Percentile 50 --High_Percentile 50 --survival_days 4673 --survival_select all_stage
	# ex : python survival_analysis.py -n KIF23 --Low_Percentile 50 --High_Percentile 50 --survival_select all_stage
	# stage argument :stage_i, stage_ii, stage_iii, stage_iv
	# parser = ArgumentParser()
	# # parser.add_argument("-p", "--project")
	# # parser.add_argument("--primary_site")
	# # parser.add_argument("-t", "--type", help="input type isoforms or genes")
	# parser.add_argument("-n", "--name", help="input isoform or gene name")
	# parser.add_argument("--Low_Percentile")
	# parser.add_argument("--High_Percentile")
	# # parser.add_argument("--survival_days")
	# parser.add_argument("--survival_select")
	# args = parser.parse_args()
	# input_project = parser.project
	input_project = "TCGA-LIHC"
	input_primary_site = "Liver_cancer"
	# input_primary_site = input_primary_site.replace("_", " ")
	input_type = "genes"
	# input_name = args.name
	# Low_Percentile = args.Low_Percentile
	# High_Percentile = args.High_Percentile
	# survival_select = args.survival_select.replace("_", " ")
	survival_select = survival_select.replace("_", " ")
	survival_days = survival_max_days(input_project, input_name, input_type, survival_select)
	# if survival_max_days(input_project, input_name, input_type, survival_select)+5 < float(survival_days) or 0 > float(survival_days):
	# 	print(f"maxmium days is {survival_max_days(input_project, input_name, input_type, survival_select)}")
	# 	print("input days error")
	# else: 
	plot_arg = {
		# 'project':args.project,
		'project':input_project,
		'primary_site': input_primary_site,
		'search_by': input_type,
		'GT_input': input_name,
		"random_id": "",
		'Low_Percentile': Low_Percentile,
		'High_Percentile': High_Percentile,
		'survival_days': survival_days,
		'survival_select': survival_select,
	}
	return survival_plot_realtime(plot_arg)

## 前端再進行說明
def get_allcancer_data(column_table,search_by):
	stage_dict = {
		'stage i' : 'stage_1',
		'stage ii' : 'stage_2',
		'stage iii' : 'stage_3',
		'stage iv' : 'stage_4',
	}
	current_path = os.path.dirname(__file__)
	db_path = f"{current_path}/../database/db"
	# db_path = f"{current_path}/db.sqlite3"
	cancer = "liver_cancer|TCGA-LIHC"
	primary_site, project = cancer.split("|")
	table_name = '%s_genes_FPKM_Cufflinks'%(project)
	column_table = "%s|%s"%(stage,table_name)
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
		print(column)
		table_name = column_table.split('|')[1]
		cursor.execute("SELECT %s FROM `%s`"%(column,table_name))
		result = list(cursor.fetchall())
	return result

if __name__ == "__main__":
	# stage_list = ['stage_1','stage_2','stage_3','stage_4', 'all stage']
	# data = get_allcancer_data('TCGA-LIHC_genes_FPKM_Cufflinks|all_stage','genes')
	sur_main()