import pandas as pd 
import numpy as np
import pymysql
from statsmodels.stats.multitest import multipletests
def Correction(data,state):
	data=data.fillna('||')
	# nan_index=data[data=='||'].index.tolist()
	data=data[~data.isin(['||'])].reset_index(drop=True)
	a=list(multipletests(data,alpha=1, method= state)[1])
	

	return a
def creat_table(df,table_name):
	conn = pymysql.connect("localhost","edward","mutter12er","edward_Cuffdiff")
	cursor = conn.cursor()
	sql_str = 'CREATE TABLE `'+str(table_name)+'`'
	column_list = list(df.columns)
	column_str = str(column_list)
	column_str = column_str.replace("'",'`').replace(',',' LONGTEXT,').replace('[','(').replace(']',' LONGTEXT)')
	
	cursor.execute(sql_str+column_str+'ENGINE=InnoDB  DEFAULT CHARSET=utf8 AUTO_INCREMENT=1;')	
	conn.commit()
def csv2mysql(df,table_name) :
	#連接帳號密碼在Putty上
	conn = pymysql.connect("localhost","edward","mutter12er","edward_Cuffdiff")

	#指標
	cursor = conn.cursor()

	values=df.values.tolist()

	s = ",".join(['%s' for _ in range(len(df.columns))])

	cursor.executemany('INSERT INTO `{}` VALUES ({})'.format(table_name,s),values)

	# 提交，不然無法儲存新建或者修改的資料
	conn.commit()

conn = pymysql.connect("localhost","edward","mutter12er","edward_Cuffdiff")
cursor = conn.cursor()
cursor.execute("SELECT TABLE_NAME  FROM INFORMATION_SCHEMA.TABLES WHERE TABLE_TYPE = 'BASE TABLE' AND TABLE_SCHEMA='edward_Cuffdiff' ")
data = cursor.fetchall()
sql_data = list(set(pd.DataFrame(list(data))[0]))
g = []
e = []
cancer = []
print(len(sql_data))
for i in sql_data:
	# if '(old)' in i and '_gather' not in i:
	g.append(i)

# for i in g:
# 	x = i.split('_')
# 	if x[1] + '_' + x[2] not in e:
# 		e.append(x[1] + '_' + x[2])
# 	if x[0] not in cancer:
# 		cancer.append(x[0])
# cancer_list = []
# e = sorted(e)
# data_list = []
# for i in cancer:
# 	genes = [i,'genes']
# 	isoforms = [i,'isoforms']
# 	for k in e:
# 		if i+'_'+k+'_genes(old)' not in g:
# 			print(i+'_'+k+'_genes')
# 			genes.append('-')
# 		if i+'_'+k+'_genes(old)' in g:
# 			cursor.execute('SELECT COUNT(*) FROM `{}_{}_genes(old)` '.format(i,k))
# 			data_count = list(cursor.fetchall()[0])[0]
# 			genes.append(data_count)
# 		if i+'_'+k+'_isoforms(old)' not in g:
# 			isoforms.append('-')
# 		if i+'_'+k+'_isoforms(old)' in g:
# 			cursor.execute('SELECT COUNT(*) FROM `{}_{}_isoforms(old)` '.format(i,k))
# 			data_count = list(cursor.fetchall()[0])[0]
# 			isoforms.append(data_count)
# 	data_list.append(genes)
# 	data_list.append(isoforms)
# columns = ['cancer','select'] + e
# data1 = pd.DataFrame(data_list,columns=columns)
# data1.to_excel('個數.xlsx',index=False)
# for i in g:
# 	cursor.execute('SELECT `isoform_name` FROM `{}` '.format(i))
# 	data = pd.DataFrame(list(cursor.fetchall()))
# 	data = set(list(data[0]))
# 	if len(gene_list) != 0:
# 		g = list(set(data)-set(gene_list))
# 		print(g)
# 		gene_list = gene_list + g
# 	else:
# 		gene_list = gene_list + list(data)
# data1 = pd.DataFrame({'isoforms':gene_list})
# data1.to_excel('isoforms.xlsx',index=False)
# for i in g:
# 	cancer = i.split('_')
	
# 	cursor.execute('SELECT COUNT(*) FROM `{}_genes_FPKM_Cuffdiff` '.format(cancer[0]))
# 	data_count = list(cursor.fetchall()[0])[0]
# 	cursor.execute('SELECT COUNT(*) FROM `{}_isoforms_FPKM_Cuffdiff` '.format(cancer[0]))
# 	data_count1 = list(cursor.fetchall()[0])[0]
# 	t = [cancer[0],data_count,data_count1]
# 	cancer_list.append(t)
# df = pd.DataFrame(cancer_list,columns=['cancer','genes','isoforms'])
# df.to_excel('count.xlsx',index=False)
# for i in e:
# 	cursor.execute("DROP TABLE `{}`".format(i))

	
for i in g:
	print(i)
	cursor.execute("SELECT * FROM `{}` WHERE 1".format(i))
	data = pd.DataFrame(list(cursor.fetchall()))
	cursor.execute("SHOW columns FROM `"+i+"`")
	data_columns = [column[0] for column in cursor.fetchall()]
	data.columns = data_columns
	if 'gather' in i:
		data_columns_list = data.columns.values.tolist()[4:]
		print(data_columns_list)
		for k in data_columns_list:
			if ('value1' not in k ) and ('value2' not in k) and ('FC' not in k):
				
				dd = list(data[k][data[k]=='nan'])
				data[k] = data[k].str.replace('nan','2')	
				data[k] = Correction(data[k].astype(float),'fdr_bh')
				data[k] = np.where(np.in1d(data[k],dd)==True,'nan',data[k])
		cursor.execute("RENAME TABLE `{}` TO `{}(old)`".format(i,i))
		creat_table(data,table_name=i)
		csv2mysql(data,table_name=i)
	else:
		dd = list(data['test_id'][data['KS_test_twosided']=='nan'])
		data['KS_test_twosided'] = data['KS_test_twosided'].str.replace('nan','2')		
		data['KS_test_twosided'] = Correction(data['KS_test_twosided'].astype(float),'fdr_bh')
		data['KS_test_twosided'] = np.where(np.in1d(data['test_id'],dd)==True,'nan',data['KS_test_twosided'])

		dd = list(data['test_id'][data['KS_test_greater']=='nan'])
		data['KS_test_greater'] = data['KS_test_greater'].str.replace('nan','2')
		data['KS_test_greater'] = Correction(data['KS_test_greater'].astype(float),'fdr_bh')
		data['KS_test_greater'] = np.where(np.in1d(data['test_id'],dd)==True,'nan',data['KS_test_greater'])

		dd = list(data['test_id'][data['KS_test_less']=='nan'])
		data['KS_test_less'] = data['KS_test_less'].str.replace('nan','2')
		data['KS_test_less'] = Correction(data['KS_test_less'].astype(float),'fdr_bh')
		data['KS_test_less'] = np.where(np.in1d(data['test_id'],dd)==True,'nan',data['KS_test_less'])

		dd = list(data['test_id'][data['T_test_twosided']=='nan'])
		data['T_test_twosided'] = data['T_test_twosided'].str.replace('nan','2')
		data['T_test_twosided'] = Correction(data['T_test_twosided'].astype(float),'fdr_bh')
		data['T_test_twosided'] = np.where(np.in1d(data['test_id'],dd)==True,'nan',data['T_test_twosided'])

		dd = list(data['test_id'][data['T_test_greater']=='nan'])
		data['T_test_greater'] = data['T_test_greater'].str.replace('nan','2')
		data['T_test_greater'] = Correction(data['T_test_greater'].astype(float),'fdr_bh')
		data['T_test_greater'] = np.where(np.in1d(data['test_id'],dd)==True,'nan',data['T_test_greater'])

		dd = list(data['test_id'][data['T_test_less']=='nan'])
		data['T_test_less'] = data['T_test_less'].str.replace('nan','2')
		data['T_test_less'] = Correction(data['T_test_less'].astype(float),'fdr_bh')
		data['T_test_less'] = np.where(np.in1d(data['test_id'],dd)==True,'nan',data['T_test_less'])

		dd = list(data['test_id'][data['U_test_twosided']=='nan'])
		data['U_test_twosided'] = data['U_test_twosided'].str.replace('nan','2')
		data['U_test_twosided'] = Correction(data['U_test_twosided'].astype(float),'fdr_bh')
		data['U_test_twosided'] = np.where(np.in1d(data['test_id'],dd)==True,'nan',data['U_test_twosided'])

		dd = list(data['test_id'][data['U_test_greater']=='nan'])
		data['U_test_greater'] = data['U_test_greater'].str.replace('nan','2')
		data['U_test_greater'] = Correction(data['U_test_greater'].astype(float),'fdr_bh')
		data['U_test_greater'] = np.where(np.in1d(data['test_id'],dd)==True,'nan',data['U_test_greater'])

		dd = list(data['test_id'][data['U_test_less']=='nan'])
		data['U_test_less'] = data['U_test_less'].str.replace('nan','2')
		data['U_test_less'] = Correction(data['U_test_less'].astype(float),'fdr_bh')
		data['U_test_less'] = np.where(np.in1d(data['test_id'],dd)==True,'nan',data['U_test_less'])
		cursor.execute("RENAME TABLE `{}` TO `{}(old)`".format(i,i))
		creat_table(data,table_name=i)
		csv2mysql(data,table_name=i)