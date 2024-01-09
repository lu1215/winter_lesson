from django.db import connections
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import os
import matplotlib.pyplot as plt
import matplotlib
from collections import OrderedDict
import time
import operator
import pandas as pd
import math
import base64
import io
import sqlite3
current_path = os.path.dirname(__file__)


stage_list = ['stage_1','stage_2','stage_3','stage_4']
stage_dict = {
	'stage i' : 'stage_1',
	'stage ii' : 'stage_2',
	'stage iii' : 'stage_3',
	'stage iv' : 'stage_4',
}

class Survival_plot():
	def survival_data_default(column_table_list,search_by,GT_input):
		cursor = connections['edward_Cufflinks'].cursor()
		primary_key = 'gene_name' if search_by == 'genes' else 'isoform_name'
		survival_result = []
		
		stage_list = ['stage_1','stage_2','stage_3','stage_4']
		table_name = column_table_list.split('|')[1]
		stage = column_table_list.split('|')[0]
		cursor.execute("SHOW columns FROM `{}`".format(table_name))
		data_columns = [column[0] for column in cursor.fetchall()]
		data_columns = data_columns[1:]
		colmns = []
		for i in data_columns:
			if i != 'normal' and i !='all_stage' and i != 'pvalue(normal_allstage)' and i != 'pvalue(allstage_normal)':
				colmns.append('`'+i+'`')

		str_colmns = ','.join(stage_list)
		print("SELECT %s FROM `%s` WHERE `%s` = '%s'"%(','.join(colmns),table_name,primary_key,GT_input))
		cursor.execute("SELECT %s FROM `%s` WHERE `%s` = '%s'"%(','.join(colmns),table_name,primary_key,GT_input))
		
		# cursor.execute("SELECT `stage_1`,`stage_2`,`stage_3`,`stage_4` FROM `%s` WHERE `%s` = '%s'"%(table_name,primary_key,GT_input))
		result = cursor.fetchall()[0]
		survival_result += [list(result)]

		return survival_result

	def survival_data_realtime(column_table,search_by,GT_input):

		table_name = column_table.split('|')[1]
		cursor = connections['edward_Cufflinks'].cursor()
		primary_key = 'gene_name' if search_by == 'genes' else 'isoform_name'
		if column_table.split('|')[0] != 'all stage':
			column = stage_dict[column_table.split('|')[0]] if column_table.split('|')[0] != 'all stage' else ','.join(stage_dict.values())
		else:
			cursor.execute("SHOW columns FROM `{}`".format(table_name))
			data_columns = [column[0] for column in cursor.fetchall()]
			data_columns = data_columns[1:]
			column=[]
			for i in data_columns:
				if i != 'normal' and i !='all_stage' and i != 'pvalue(normal_allstage)' and i != 'pvalue(allstage_normal)':
					column.append('`'+i+'`')
			column = ','.join(column)
		print("SELECT %s FROM `%s` WHERE `%s` = '%s'"%(column,table_name,primary_key,GT_input))
		cursor.execute("SELECT %s FROM `%s` WHERE `%s` = '%s'"%(column,table_name,primary_key,GT_input))
		print("SELECT %s FROM `%s` WHERE `%s` = '%s'"%(column,table_name,primary_key,GT_input))
		result = list(cursor.fetchall()[0])
		return result
		
	def survival_plot(T1,E1,T2,E2,GT_input,primary_site,random_id,Low_Percentile,High_Percentile,survival_days,survival_select):
		plot_path = f"{current_path}/../static/Survival_plot/survival_plot_%s_%s.png"%(primary_site.replace(' ','_').replace("(",'').replace(")",''),random_id)
		try:
			os.remove(plot_path)
		except OSError:
			pass
		kmf = KaplanMeierFitter()
		dpi = 100
		# T1 = [232,565,205,4967,370,276,59,1649,819,590,82,522,3364,640,211,223,783,1556,3981,1005,638,1949,615,823,484,1884,899,588,539,562,344,88,413,789,273,311,1639,510,278,37,20,122,851,118,2423,1048,540,547,685,610,237,89,406,603,572,272,415,303,377,599,328,248,272,508,539,62,579,1108,2109,251,974,384,154,633,173,2027,467,1110,382,5041,1582,1830,491,128,90,372,474,522,200,4343,253,272,480,646,758,408,400,251,2828,1947]
		# E1 = [True,True,True,False,True,False,False,False,True,True,False,False,False,False,True,True,False,True,False,True,False,False,True,True,False,False,False,False,False,False,True,True,True,False,True,True,False,True,True,False,True,True,False,True,False,False,False,True,True,False,True,False,True,False,False,True,True,True,False,True,True,True,True,True,True,True,True,False,False,True,True,False,True,False,True,False,True,False,False,False,False,False,False,True,True,False,True,True,True,False,True,True,False,False,False,True,True,False,True,False] #是否死亡
		# T2 = [13,1326,906,474,1090,76,469,2177,95,55,680,544,149,544,144,2008,542,129,2954,734,415,460,182,385,467,1845,324,578,949,240,3183,399,20,262,612,2380,361,453,428,370,168,1029,1912,366,369,105,1561,642,2703,859,333,410,1348,3011,1350,455,700,921,1792,690,641,832,577,798,324,512,117,565,691,57,1869,1064,389,415,384,1454,28,1072,1621,368,385,2139,864,56,2020,345,1604,1804,365,388,508,364,224,2312,142,997,1522,67,294,1529]
		# E2 = [False,False,False,False,False,True,False,False,False,False,True,True,True,True,True,False,False,False,True,True,False,True,True,True,True,False,True,False,True,False,True,False,False,True,True,False,False,True,False,False,True,False,False,False,False,False,False,False,False,True,False,False,True,False,False,False,False,False,False,True,False,False,True,False,True,False,False,True,False,True,True,True,False,False,False,False,False,False,False,False,True,False,False,True,False,False,False,True,False,True,False,False,False,False,True,False,False,False,True,False]
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
		return output_data
		
class Others():
	def primary_project_stage():
		cursor = connections['default'].cursor()

		primary_project_stage = {}
		cursor.execute("SELECT DISTINCT `primary_site`,`project`,`condition1`,`condition2` FROM `Mutual_Relationship`")
		result2 = cursor.fetchall()
		for row in result2:
			primary_project = "%s|%s"%(row[0],row[1])
			if primary_project not in primary_project_stage:
				primary_project_stage[primary_project] = [row[2],row[3]]
			else:
				if row[2] not in primary_project_stage[primary_project]:
					primary_project_stage[primary_project] += [row[2]]
				if row[3] not in primary_project_stage[primary_project]:
					primary_project_stage[primary_project] += [row[3]]
		return primary_project_stage

class Detail():
	def belong_gene(search_by,GT_input):
		cursor = connections['default'].cursor()
		sql = "SELECT `gene` FROM `hg38_gene_transcripts_20180130` WHERE `transcripts` = '%s' OR `transcripts` LIKE '%s' OR `transcripts` LIKE '%s' OR `transcripts` LIKE '%s'"%(GT_input,"%,"+GT_input,GT_input+",%","%,"+GT_input+",%")
		cursor.execute(sql)
		result = cursor.fetchall()[0][0]
		return result
	def all_transcript(gene):
		# cursor_default = connections['default'].cursor()
		db_path = f"{current_path}/../../../database/db"
		with sqlite3.connect(db_path, check_same_thread=False) as db_conn:
			cursor_default = db_conn.cursor()
			cursor_default.execute("SELECT `transcripts` FROM `hg38_gene_transcripts_20180130` WHERE `gene` = '%s'"%(gene))
		transcript_list = sorted(cursor_default.fetchall()[0][0].split(','))
		return transcript_list

	def diff_data_single(TEST_column,table_name_cufflinks,gene_isoform,DE_level):
		cursor = connections['edward_Cuffdiff'].cursor()
		if TEST_column:
			if DE_level == 'genes':
				sql = "SELECT `test_id`,`value_1`,`value_2`,`log2(fold_change)`,`%s` FROM `%s` WHERE `test_id` = '%s'"%(TEST_column,table_name_cufflinks,gene_isoform)
			elif DE_level == 'isoforms':
				sql = "SELECT `test_id`,`gene_id`,`value_1`,`value_2`,`log2(fold_change)`,`%s` FROM `%s` WHERE `test_id` = '%s'"%(TEST_column,table_name_cufflinks,gene_isoform)
		else:
			if DE_level == 'genes':
				sql = "SELECT `test_id`,`value_1`,`value_2`,`log2(fold_change)` FROM `%s` WHERE `test_id` = '%s'"%(table_name_cufflinks,gene_isoform)
			elif DE_level == 'isoforms':
				sql = "SELECT `test_id`,`gene_id`,`value_1`,`value_2`,`log2(fold_change)` FROM `%s` WHERE `test_id` = '%s'"%(table_name_cufflinks,gene_isoform)

		cursor.execute(sql)
		result = list(cursor.fetchall()[0])
		if TEST_column:
			result[-1] = float(result[-1])
			result[-2] = pow(2,float(result[-2]))
		else:
			result[-1] = pow(2,float(result[-1]))
		
		return result

	def diff_data_all_transcript(TEST_column,transcript_list,table_name_cuffdiff_isoforms):

		cursor_Cuffdiff = connections['edward_Cuffdiff'].cursor()
		if TEST_column:
			if len(transcript_list) == 1:
				sql = "SELECT `test_id`,`value_1`,`value_2`,`log2(fold_change)`,`%s` FROM `%s` WHERE `test_id` IN ('%s') ORDER BY `%s`"%(TEST_column,table_name_cuffdiff_isoforms,transcript_list[0],'test_id')
			else:
				sql = "SELECT `test_id`,`value_1`,`value_2`,`log2(fold_change)`,`%s` FROM `%s` WHERE `test_id` IN %s ORDER BY `%s`"%(TEST_column,table_name_cuffdiff_isoforms,tuple(transcript_list),'test_id')
		else:
			if len(transcript_list) == 1:
				sql = "SELECT `test_id`,`value_1`,`value_2`,`log2(fold_change)` FROM `%s` WHERE `test_id` IN ('%s') ORDER BY `%s`"%(table_name_cuffdiff_isoforms,transcript_list[0],'test_id')
			else:
				sql = "SELECT `test_id`,`value_1`,`value_2`,`log2(fold_change)` FROM `%s` WHERE `test_id` IN %s ORDER BY `%s`"%(table_name_cuffdiff_isoforms,tuple(transcript_list),'test_id')
		cursor_Cuffdiff.execute(sql)
		result = list(map(list,cursor_Cuffdiff.fetchall()))
		if TEST_column:
			for row in result:
				row[-1] = float(row[-1])
				row[-2] = pow(2,float(row[-2]))
				if row[0] == transcript_list[0]:
					result.insert(0, result.pop(result.index(row))) # move isoform to first
		else:
			for row in result:
				row[-1] = pow(2,float(row[-1]))
				if row[0] == transcript_list[0]:
					result.insert(0, result.pop(result.index(row))) # move isoform to first
		return result
	def struture_info_all_transcript(transcript_list):
		cursor = connections['default'].cursor()
		if len(transcript_list) == 1:
			sql = "SELECT * FROM `Isoform_struture_info` WHERE `isoform_name` IN ('%s') ORDER BY `isoform_name`"%transcript_list[0]
		else:
			sql = "SELECT * FROM `Isoform_struture_info` WHERE `isoform_name` IN {0} ORDER BY `isoform_name`".format(tuple(transcript_list))
		cursor.execute(sql)
		result = list(map(list,cursor.fetchall()))
		for row in result:
			if row[0] == transcript_list[0]:
				result.insert(0, result.pop(result.index(row))) # move isoform to first
		return result
			
class Summary():
	def NCBI_gene_summary(gene):
		db_path = f"{current_path}/../../../database/db"
		with sqlite3.connect(db_path, check_same_thread=False) as db_conn:
			cursor_summary = db_conn.cursor()
			# cursor_summary = connections['default'].cursor()
			cursor_summary.execute("SELECT `Symbol`,`Synonyms`,`chromosome`,`description`,`Other_designations` FROM `NCBI_gene_info_20180201` WHERE `Symbol` = '%s'"%(gene))
			result_rows = cursor_summary.fetchall()
			summary_result = result_rows[0] if result_rows else ''
		return summary_result

	def NCBI_transcript_summary(transcript_list):
		db_path = f"{current_path}/../../../database/db"
		with sqlite3.connect(db_path, check_same_thread=False) as db_conn:
			cursor_default = db_conn.cursor()
			sql = "SELECT * FROM `NCBI_transcript_info_20180209` WHERE "
			for idx, transcript in enumerate(transcript_list):
				if idx == len(transcript_list)-1:
					sql += "`transcript_name` = '%s'"%transcript
				else:
					sql += "`transcript_name` = '%s' OR "%transcript
			cursor_default.execute(sql)
			result = cursor_default.fetchall()
		NCBI_transcript_summary = {}
		for row in result:
			NCBI_transcript_summary[row[0]] = row
		return NCBI_transcript_summary
		
class Filter():
	def project_primary_dict():
		db_path = f"{current_path}/../../../database/db"
		with sqlite3.connect(db_path, check_same_thread=False) as db_conn:
			cursor = db_conn.cursor()
			cursor = connections['default'].cursor()
			cursor.execute("SELECT DISTINCT `primary_site`,`project` FROM `Mutual_Relationship` WHERE `project` != 'TCGA-BRCA'") ## Breast 先關掉 還沒有資料
			result = cursor.fetchall()
		project_primary_dict = OrderedDict()
		for row in result:
			if row[0] not in project_primary_dict:
				project_primary_dict[row[1]] = row[0]

		number_of_sample = {}
		cursor.execute("SELECT DISTINCT `project`,`condition1`,`#_of_normal`,`#_of_stage_1`,`#_of_stage_2`,`#_of_stage_3`,`#_of_stage_4` FROM `Mutual_Relationship` WHERE `project` != 'TCGA-BRCA'")## Breast 先關掉 還沒有資料
		result2 = cursor.fetchall()
		for row in result2:
			if row[0] not in number_of_sample:
				number_of_sample[row[0]] = {
					'normal':row[2],
					'stage i':row[3],
					'stage ii':row[4],
					'stage iii':row[5],
					'stage iv':row[6],
				}
		return [project_primary_dict,number_of_sample]
	
	def filter(search_by,GT_input,primary_key,table_list,FC_select,FC_input,TEST_column,TEST_input):

		stage_dict = {
			'N':'normal',
			'1':'stage i',
			'2':'stage ii',
			'3':'stage iii',
			'4':'stage iv',
		}

		cursor = connections['edward_Cuffdiff'].cursor()
		result_dict = OrderedDict()
		table_list = list(map(lambda x:"%s_gather_%s"%(x,search_by),table_list))

		for table_name in table_list:
			
			table_detect = cursor.execute("show tables like '%s'"%table_name)
			if table_detect == 1:
				sql = "SELECT * FROM `%s` WHERE `%s` = '%s'"%(table_name,primary_key,GT_input)
				print(sql)
				row_count = cursor.execute(sql)
				if row_count == 0:
					table_condition_FC_qv_list = []
					break
		
				result = cursor.description

				result_dict_pre =  [OrderedDict(zip([col[0] for col in result], row)) for row in cursor.fetchall()][0]
				result_dict[table_name] = OrderedDict()
				for condition_pre,value in result_dict_pre.items():
					if condition_pre != primary_key:
						condition = '_'.join(condition_pre.split('_')[:-1]) if 'test' not in condition_pre else '_'.join(condition_pre.split('_')[:-3])
						value_type = condition_pre.split('_')[-1] if 'test' not in condition_pre else '_'.join(condition_pre.split('_')[-3:])
						if condition not in result_dict[table_name]:
							result_dict[table_name][condition] = OrderedDict()
							value = float("inf") if value == 'inf' else float(value)
							value = float("-inf") if value == '-inf' else float(value)
							value = int(value) if value.is_integer() else value
							result_dict[table_name][condition][value_type] = value
						else:
							if value_type not in result_dict[table_name][condition]:
								value = float("inf") if value == 'inf' else float(value)
								value = float("-inf") if value == '-inf' else float(value)
								value = int(value) if value.is_integer() else value
								result_dict[table_name][condition][value_type] = value

		table_condition_FC_qv_list = []
		for table_name,info in result_dict.items():
			cursor.execute('SELECT COUNT(*) FROM `{}`'.format(table_name))
			data_count = list(cursor.fetchall()[0])[0]
			print(TEST_input)
			project = table_name.split('_')[0]
			for condition,value in info.items():

				if FC_input != '' and TEST_input != '':
					
					if FC_select == '≥':
						if pow(2,value['FC']) >= FC_input and value[TEST_column] <= float(TEST_input) :
							condition_fullname1 = stage_dict[condition.split('_')[0]]
							condition_fullname2 = stage_dict[condition.split('_')[1]]
							if value[TEST_column] >= 1:
								value[TEST_column] = 1
							else :
								value[TEST_column] = value[TEST_column]
							table_condition_FC_qv_list += [[project,condition_fullname1,condition_fullname2,value['value1'],value['value2'],pow(2,value['FC']),value[TEST_column]]]
					elif FC_select == '≤':
						if pow(2,value['FC']) <= FC_input and value[TEST_column] <= float(TEST_input) :
							condition_fullname1 = stage_dict[condition.split('_')[0]]
							condition_fullname2 = stage_dict[condition.split('_')[1]]
							if value[TEST_column] >= 1:
								value[TEST_column] = 1
							else :
								value[TEST_column] = value[TEST_column]
							table_condition_FC_qv_list += [[project,condition_fullname1,condition_fullname2,value['value1'],value['value2'],pow(2,value['FC']),value[TEST_column]]]	
					elif FC_select == 'Either':
						if (pow(2,value['FC']) >= FC_input or pow(2,value['FC']) <= FC_input) and value[TEST_column] <= float(TEST_input) :
							condition_fullname1 = stage_dict[condition.split('_')[0]]
							condition_fullname2 = stage_dict[condition.split('_')[1]]
							if value[TEST_column] >= 1:
								value[TEST_column] = 1
							else :
								value[TEST_column] = value[TEST_column]
							table_condition_FC_qv_list += [[project,condition_fullname1,condition_fullname2,value['value1'],value['value2'],pow(2,value['FC']),value[TEST_column]]]

				elif FC_input != '' and TEST_input == '':
					if FC_select == '≥':
						if pow(2,value['FC']) >= FC_input:
							condition_fullname1 = stage_dict[condition.split('_')[0]]
							condition_fullname2 = stage_dict[condition.split('_')[1]]
							table_condition_FC_qv_list += [[project,condition_fullname1,condition_fullname2,value['value1'],value['value2'],pow(2,value['FC'])]]

					elif FC_select == '≤':
						if pow(2,value['FC']) <= FC_input:
							condition_fullname1 = stage_dict[condition.split('_')[0]]
							condition_fullname2 = stage_dict[condition.split('_')[1]]
							table_condition_FC_qv_list += [[project,condition_fullname1,condition_fullname2,value['value1'],value['value2'],pow(2,value['FC'])]]

					elif FC_select == 'Either':
						if (pow(2,value['FC']) >= FC_input or pow(2,value['FC']) <= FC_input):
							condition_fullname1 = stage_dict[condition.split('_')[0]]
							condition_fullname2 = stage_dict[condition.split('_')[1]]
							table_condition_FC_qv_list += [[project,condition_fullname1,condition_fullname2,value['value1'],value['value2'],pow(2,value['FC'])]]

				elif FC_input == '' and TEST_input != '':
					
					if value[TEST_column] <= float(TEST_input) :
						condition_fullname1 = stage_dict[condition.split('_')[0]]
						condition_fullname2 = stage_dict[condition.split('_')[1]]
						if value[TEST_column] >= 1:
							value[TEST_column] = 1
						else :
							value[TEST_column] = value[TEST_column]
						table_condition_FC_qv_list += [[project,condition_fullname1,condition_fullname2,value['value1'],value['value2'],pow(2,value['FC']),value[TEST_column]]]

				else:
					condition_fullname1 = stage_dict[condition.split('_')[0]]
					condition_fullname2 = stage_dict[condition.split('_')[1]]
					table_condition_FC_qv_list += [[project,condition_fullname1,condition_fullname2,value['value1'],value['value2'],pow(2,value['FC'])]]
		print(table_condition_FC_qv_list)
		table_condition_FC_qv_list = sorted(table_condition_FC_qv_list, key=operator.itemgetter(1, 2))#sort by condition1 and condition2
		# print(table_condition_FC_qv_list)
		return table_condition_FC_qv_list

	def gene_link(gene_list,DE_level):
		# cursor = connections['default'].cursor()
		db_path = f"{current_path}/../../../database/db"
		with sqlite3.connect(db_path, check_same_thread=False) as db_conn:
			cursor = db_conn.cursor() 
			cursor.execute("SELECT `Gene_link`,`Symbol`,`Synonyms` FROM `NCBI_gene_info_20180201`")
			result = cursor.fetchall()
		genelink_dict_Symbol = {}
		for row in result:
			genelink_dict_Symbol[row[1]] = row[0]
		genelink_dict_Synonyms = {}
		for row in result:
			for Synonyms in row[1].split('|'):
				if Synonyms not in genelink_dict_Synonyms and Synonyms != '-':
					genelink_dict_Synonyms[Synonyms] = row[0]
		genelink_dict = {}
		if DE_level == 'genes':
			for gene in gene_list:
				if gene in genelink_dict_Symbol:
					genelink = genelink_dict_Symbol[gene]
				elif gene in genelink_dict_Synonyms:
					genelink = genelink_dict_Synonyms[gene]
				else:
					genelink = "https://www.ncbi.nlm.nih.gov/gene/?term=%s"%gene

				genelink_dict[gene] = genelink

		elif DE_level == 'isoforms':
			for gene in gene_list:
				if gene in genelink_dict_Symbol:
					genelink = genelink_dict_Symbol[gene]
				elif gene in genelink_dict_Synonyms:
					genelink = genelink_dict_Synonyms[gene]
				else:
					genelink = "https://www.ncbi.nlm.nih.gov/gene/?term=%s"%gene

				genelink_dict[gene] = genelink

		return genelink_dict

	def filter_logFC_Pvalue(TEST_column,table_name_cuffdiff,FC_select,FC_input,TEST_input,DE_level):
		# if DE_level == 'genes':
		#	 sql = "SELECT `test_id`,`value_1`,`value_2`,pow(2,`log2(fold_change)`),`q_value` FROM `%s`"%(table_name_cuffdiff)
		# elif DE_level == 'isoforms':	
		#	 sql = "SELECT `test_id`,`gene_id`,`value_1`,`value_2`,pow(2,`log2(fold_change)`),`q_value` FROM `%s`"%(table_name_cuffdiff)

		# if FC_input != '' and qvalue_input == '':
		#	 if FC_select == '≥':
		#		 sql += " WHERE pow(2,`log2(fold_change)`) >= %s OR `log2(fold_change)` = 'inf'"%(FC_input)
		#	 elif FC_select == '≤':
		#		 sql += " WHERE pow(2,`log2(fold_change)`) <= %s OR `log2(fold_change)` = '-inf'"%(FC_input)
		#	 elif FC_select == 'Either':
		#		 sql += " WHERE pow(2,`log2(fold_change)`) >= %s OR pow(2,`log2(fold_change)`) <= %s OR `log2(fold_change)` = 'inf' OR `log2(fold_change)` = '-inf'"%(FC_input,FC_input)

		# elif FC_input == '' and qvalue_input != '':
		#	 sql += " WHERE `q_value` <= %s"%(qvalue_input)

		# elif FC_input != '' and qvalue_input != '':
		#	 if FC_select == '≥':
		#		 sql += " WHERE (pow(2,`log2(fold_change)`) >= %s OR `log2(fold_change)` = 'inf') AND `q_value` <= %s"%(FC_input,qvalue_input)
		#	 elif FC_select == '≤':
		#		 sql += " WHERE (pow(2,`log2(fold_change)`) <= %s OR `log2(fold_change)` = '-inf') AND `q_value` <= %s"%(FC_input,qvalue_input)
		#	 elif FC_select == 'Either':
		#		 sql += " WHERE (pow(2,`log2(fold_change)`) >= %s OR pow(2,`log2(fold_change)`) <= %s OR `log2(fold_change)` = 'inf' OR `log2(fold_change)` = '-inf') AND `q_value` <= %s"%(FC_input,FC_input,qvalue_input)
		# else:
		#	 sql = sql
		cursor = connections['edward_Cuffdiff'].cursor()
		start = time.time()
		
		if TEST_input:
			
			if DE_level == 'genes':
				sql = "SELECT `test_id`,`value_1`,`value_2`,`log2(fold_change)`,`%s`,ROUND( (`value_1` / `value_2`), 3 ) AS `fold_change` FROM `%s`"%(TEST_column,table_name_cuffdiff)
			elif DE_level == 'isoforms':	
				sql = "SELECT `test_id`,`gene_id`,`value_1`,`value_2`,`log2(fold_change)`,`%s`,ROUND( (`value_1` / `value_2`), 3 ) AS `fold_change` FROM `%s`"%(TEST_column,table_name_cuffdiff)
		else:
			if DE_level == 'genes':
				sql = "SELECT `test_id`,`value_1`,`value_2`,`log2(fold_change)`,ROUND( (`value_1` / `value_2`), 3 ) AS `fold_change` FROM `%s`"%(table_name_cuffdiff)
			elif DE_level == 'isoforms':	
				sql = "SELECT `test_id`,`gene_id`,`value_1`,`value_2`,`log2(fold_change)`,ROUND( (`value_1` / `value_2`), 3 ) AS `fold_change` FROM `%s`"%(table_name_cuffdiff)

		if FC_input != '' and TEST_input == '':
			if FC_select == '≥':
				sql += " WHERE (ROUND(`value_1`,3) >= 0.001 AND ROUND(`value_2`,3) >= 0.001) AND (ROUND(`value_2`,3)/ROUND(`value_1`,3) >= %s )"%(FC_input)
			elif FC_select == '≤':
				sql += " WHERE (ROUND(`value_1`,3) >= 0.001 AND ROUND(`value_2`,3) >= 0.001) AND (ROUND(`value_2`,3)/ROUND(`value_1`,3) <= %s )"%(FC_input)
			elif FC_select == 'Either':
				sql += " WHERE (ROUND(`value_1`,3) >= 0.001 AND ROUND(`value_2`,3) >= 0.001) AND (pow(2,`log2(fold_change)`) >= %s OR pow(2,`log2(fold_change)`) <= %s OR `log2(fold_change)` = 'inf' OR `log2(fold_change)` = '-inf')"%(FC_input,FC_input)

		elif FC_input == '' and TEST_input != '':
			
			sql += " WHERE (ROUND(`value_1`,3) >= 0.001 AND ROUND(`value_2`,3) >= 0.001) AND `%s` <= %s"%(TEST_column,TEST_input)

		elif FC_input != '' and TEST_input != '':
			
			if FC_select == '≥':
				sql += " WHERE (ROUND(`value_1`,3) >= 0.001 AND ROUND(`value_2`,3) >= 0.001) AND (((ROUND(`value_2`,3)/ROUND(`value_1`,3)) >= %s ) AND `%s` <= %s)"%(FC_input,TEST_column,TEST_input)
			elif FC_select == '≤':
				sql += " WHERE (ROUND(`value_1`,3) >= 0.001 AND ROUND(`value_2`,3) >= 0.001) AND (((ROUND(`value_2`,3)/ROUND(`value_1`,3)) <= %s ) AND `%s` <= %s)"%(FC_input,TEST_column,TEST_input)
			elif FC_select == 'Either':
				sql += " WHERE (ROUND(`value_1`,3) >= 0.001 AND ROUND(`value_2`,3) >= 0.001) AND ((pow(2,`log2(fold_change)`) >= %s OR pow(2,`log2(fold_change)`) <= %s OR `log2(fold_change)` = 'inf' OR `log2(fold_change)` = '-inf') AND `%s` <= %s)"%(FC_input,FC_input,TEST_column,TEST_input)
		else:
			sql += " WHERE (ROUND(`value_1`,3) >= 0.001 AND ROUND(`value_2`,3) >= 0.001) "


		print(sql)
		cursor.execute(sql)
		table_data = list(map(list,cursor.fetchall()))
		end = time.time()
		download_table_data = []
		start = time.time()
		
		if table_data:
			if DE_level =='genes':
				if TEST_input:
					for row in table_data:
						if float(row[-2])  >= 1 :
							row[-2] = 1
						else:
							row[-2] = float(row[-2]) 
						if (row[-3] == 'inf' or row[-3] == '-inf'):
							row[-3] = float(row[-3])
						else :
							row[-3] = row[-1]
						row[-4] = float(row[-4])
						row[-5] = float(row[-5])
						row.pop()
						
						download_table_data.append([row[0],round(row[2],3),round(row[1],3),round(round(row[2],3)/round(row[1],3),3),row[4]])
				else:
					for row in table_data:
						if (row[-2] == 'inf' or row[-2] == '-inf'):
							row[-2] = float(row[-2])
						else :
							row[-2] = row[-1]
						row[-3] = float(row[-3])
						row[-4] = float(row[-4])
						row.pop()
						
						download_table_data.append([row[0],round(row[2],3),round(row[1],3),round(round(row[2],3)/round(row[1],3),3)])
						
			else:
				path = '/home/edward/Django/edward_project/static/error_isoforms.txt'
				f = open(path, 'r')
				error_isoforms = eval(f.read())
				f.close()
				table_data = pd.DataFrame(table_data)
				table_data = table_data[~table_data[0].isin(error_isoforms)]
				table_data = table_data.values.tolist()
				if TEST_input:
					for row in table_data:
						if float(row[-2])  >= 1 :
							row[-2] = 1
						else:
							row[-2] = float(row[-2]) 
						if (row[-3] == 'inf' or row[-3] == '-inf'):
							row[-3] = float(row[-3])
						else :
							row[-3] = row[-1]
						row[-4] = float(row[-4])
						row[-5] = float(row[-5])
						row.pop()
						download_table_data.append([row[0],row[1],round(row[3],3),round(row[2],3),round(round(row[3],3)/round(row[2],3),3)])
				else:
					for row in table_data:
						if (row[-2] == 'inf' or row[-2] == '-inf'):
							row[-2] = float(row[-2])
						else :
							row[-2] = row[-1]
						row[-3] = float(row[-3])
						row[-4] = float(row[-4])
						row.pop()
						
						download_table_data.append([row[0],row[1],round(row[3],3),round(row[2],3),round(round(row[3],3)/round(row[2],3),3)])
						
			cursor.close()
			if TEST_input:
				table_data = sorted(table_data,key=operator.itemgetter(4, 3))
			else:
				table_data = sorted(table_data,key=operator.itemgetter(3, 2), reverse=True)
			end = time.time()
		
		return table_data,download_table_data

	def download_table(diff_data,select_filter_value):

		primary_site = select_filter_value[1].split('|')[0]
		project = select_filter_value[1].split('|')[1]
		condition1 = select_filter_value[2].split('|')[0]
		condition2 = select_filter_value[3].split('|')[0]
		DE_level = select_filter_value[0]
		FC_select = select_filter_value[4].split(' (')[0].strip()
		FC_input = select_filter_value[5].strip()
		TEST_select = select_filter_value[6]
		TESTstates_select = select_filter_value[7].split(' (')[0].strip()
		TEST_input = select_filter_value[8].strip()

		output_data = [
						['%s Differential Expression %s'%(len(diff_data),DE_level.title())],
						['Primary site',primary_site],
						['Condition1',condition2],
						['Condition2',condition1],
						['Differential Expression level',"DE %s"%DE_level],
						]
		if FC_input != '' and TEST_input == '':
			output_data += [
								['Fold Change (%s/%s)'%(condition2,condition1),"%s %s Fold"%(FC_select,FC_input)],
								['Hypothesis Test',"-"],
								[],
								[DE_level[:-1].title(),'Condition1 Avg FPKM','Condition2 Avg FPKM','Fold Change']
							]
		elif FC_input == '' and TEST_input != '':
			output_data += [
								['Fold Change (%s/%s)'%(condition2,condition1),"-"],
								['Hypothesis Test',"%s (%s) significance level: %s"%(TEST_select,TESTstates_select,TEST_input)],
								[],
								[DE_level[:-1].title(),'Condition2 Avg FPKM','Condition1Avg FPKM','Fold Change',"%s (%s)"%(TEST_select.replace('_',' '),TESTstates_select)]
							]
		elif FC_input != '' and TEST_input != '':
			output_data += [
								['Fold Change (%s/%s)'%(condition2,condition1),"%s %s Fold"%(FC_select,FC_input)],
								['Hypothesis Test',"%s (%s) significance level: %s"%(TEST_select,TESTstates_select,TEST_input)],
								[],
								[DE_level[:-1].title(),'Condition2 Avg FPKM','Condition1 Avg FPKM','Fold Change',"%s (%s)"%(TEST_select.replace('_',' '),TESTstates_select)]
							]
		else:
			output_data += [
								['Fold Change (%s/%s)'%(condition2,condition1),"-"],
								['Hypothesis Test',"-"],
								[],
								[DE_level[:-1].title(),'Condition2 Avg FPKM','Condition1 Avg FPKM','Fold Change']
							]
		print(diff_data)
		output_data = output_data + diff_data
		# for row in diff_data:
		#	 output_data += [row]

		return output_data

	def heatmap_data(diff_data,project,DE_level,condition1,condition2):

		stage_Arabic_dict = {
			'normal' :'normal',
			'stage i' : 'stage_1',
			'stage ii' : 'stage_2',
			'stage iii' : 'stage_3',
			'stage iv' : 'stage_4',
		}
		stage_dict = {
			'normal' :'N',
			'stage i' : '1',
			'stage ii' : '2',
			'stage iii' : '3',
			'stage iv' : '4',
		}
		gene_list = [x[0] for x in diff_data]
		column1 = "%s(%s_%s)"%(stage_Arabic_dict[condition1],stage_dict[condition1],stage_dict[condition2])
		column2 = "%s(%s_%s)"%(stage_Arabic_dict[condition2],stage_dict[condition1],stage_dict[condition2])
		table_name = "%s_%s_FPKM_Cuffdiff"%(project,DE_level)
		if len(gene_list) == 1:
			sql = "SELECT `%s_name`,`%s`,`%s` FROM `%s` WHERE `%s_name` IN ('%s') ORDER BY `%s_name`"%(DE_level[:-1],column1,column2,table_name,DE_level[:-1],gene_list[0],DE_level[:-1])
		else:
			sql = "SELECT `%s_name`,`%s`,`%s` FROM `%s` WHERE `%s_name` IN %s ORDER BY `%s_name`"%(DE_level[:-1],column1,column2,table_name,DE_level[:-1],tuple(gene_list),DE_level[:-1])
		cursor = connections['edward_Cufflinks'].cursor()
		cursor.execute(sql)
		result = list(map(list,cursor.fetchall()))
		result_log2_add1 = []
		for row in result:
			GI = row[0]
			FPKM_list1 = row[1].split(',')
			FPKM_list2 = row[2].split(',')
			FPKM_list1_new = map(lambda x:str(math.log2(float(x)+1)),FPKM_list1)
			FPKM_list2_new = map(lambda x:str(math.log2(float(x)+1)),FPKM_list2)
			result_log2_add1 += [[GI,','.join(FPKM_list1_new),','.join(FPKM_list2_new)]]

		return result_log2_add1
	
class Boxplot():
	def boxplot_data(primary_stage,search_by,GT_input):

		cursor = connections['edward_Cufflinks'].cursor()
		primary_key = 'gene_name' if search_by == 'genes' else 'isoform_name'
		primary_nodata = []
		primary_FPKM = OrderedDict()
		for primary_project,stage_list in primary_stage.items():

			project = primary_project.split('|')[1]
			table_name = '%s_%s_FPKM_Cuffdiff'%(project,search_by)
			cursor.execute("SELECT * FROM `%s` WHERE `%s` = '%s'"%(table_name,primary_key,GT_input))

			# result = cursor.fetchall()
			result = cursor.description
			T = cursor.fetchall()
			
			if len(list(T)) !=0:
				print("SELECT * FROM `%s` WHERE `%s` = '%s'"%(table_name,primary_key,GT_input))
				primary_FPKM[primary_project] = [dict(zip([col[0] for col in result], row)) for row in T][0]
			else:
				primary_nodata.append(primary_project)


		EXP_result = []
		for primary_project,stage_list in primary_stage.items():
			if primary_project not in primary_nodata:
				EXP_pre = []
				for stage in stage_list:
					stage1 = stage[0]
					stage2 = stage[1]
					EXP_pre += [list(map(float,primary_FPKM[primary_project][stage1].split(',')))]
					EXP_pre += [list(map(float,primary_FPKM[primary_project][stage2].split(',')))]
				EXP_result += [EXP_pre]
		
		return EXP_result,primary_nodata