import scipy
from statsmodels.stats.multitest import multipletests
import pandas as pd
import os


current_path = os.path.dirname(__file__)
def enrichment(seq_data:str, correction:str = "None", p_limit:float = 1):
    # seq_data = request.POST["seq"]
    # correction = request.POST["Correction"]
    # p_limit = float(request.POST["p_limit"])
    # print(type(p_limit), p_limit)

    def fisher(A,B,C,D) :
        T = int(A)      # 交集清單 1 剩下的交集處 
        S = int(B)      # 輸入genes數 18 輸入一總數
        G = int(C)      # miRNA target數 1230 篩選的樣本 
        F = int(D)      # homo 總 genes數 26380 輸入二總數

        S_T = S-T
        G_T = G-T
        F_G_S_T = F-G-S+T

        oddsratio, pvalue_greater = scipy.stats.fisher_exact( [ [T,G_T] , [S_T,F_G_S_T]] ,'greater')
        oddsratio, pvalue_less = scipy.stats.fisher_exact( [ [T,G_T] , [S_T,F_G_S_T]] ,'less')

        return pvalue_greater,G,F

    # ================================================================================================
    # ================================================================================================
    # query_domain = "go_f_map_id"
    # domain_data = pd.read_csv("data/{}.csv".format(query_domain))
    domain_data = pd.read_csv(f"{current_path}/../static/data/enrichment_data/miRNA_domain_map_id.csv")
    # domain_data = pd.read_csv("miRNA_domain_map_id.csv")

    length = len(domain_data)

    input_list = seq_data.replace(" ","").replace("'","").replace("[","").replace("]","").split(",")

    D = [26380 for n in range(length)]
    C = list(domain_data["count"])
    B = [len(input_list) for n in range(length)]
    list_A = list(range(length))
    A = list(range(length))
    test = list(range(length))
    # for n in range(len(domain_data["go_f"])):
    for n in range(len(domain_data["mirna_name"])):

        list_A[n] = list(set(input_list)&set(domain_data["gene_name"][n].replace("'","").replace(" ","").replace("[","").replace("]","").split(",")))
        # print(domain_data["gene_name"][n])
        A[n] = len(list_A[n])

        test[n] = fisher(A[n], B[n], C[n], D[n])[0]

    cut_off = 0.01
    P_value_corr_FDR = multipletests(test,alpha=cut_off, method= "fdr_bh")
    P_value_corr_Bon = multipletests(test,alpha=cut_off, method= "bonferroni")

    # result = pd.DataFrame({
    #     "Domain_id":domain_data["go_f"], "P-value":test, "FDR":P_value_corr_FDR[1], "Bonferroni":P_value_corr_Bon[1], "A": A, "B":B, "C":C, "D":D,
    # })
    result = pd.DataFrame({
        "Domain_id":domain_data["mirna_name"],"gene_name":domain_data["gene_name"],"P-value":test,"FDR":P_value_corr_FDR[1],"Bonferroni":P_value_corr_Bon[1], "A": A, "B":B, "C":C, "D":D,
    })

    # expected_ratio = "{} / {} ( {} %)".format(C,D,(C/D))

    # print(len(P_value_corr_FDR))
    # print("[ [小於cut off 回傳True], [校正後的P-value], [corrected alpha for Sidak method], [corrected alpha for Bonferroni method] ]")

    # print(result)
    # print(result[result["FDR"]<=0.01])
    if correction == "None":
        result = result[result["P-value"] <= p_limit].reset_index(drop=True)
    else:
        result = result[result[correction] <= p_limit].reset_index(drop=True)

    observed_ratio = result["A"]/result["B"]
    expected_ratio = result["C"]/result["D"]
    result.insert(7,"observed_ratio", observed_ratio)
    result.insert(10,"expected_ratio", expected_ratio)
    column_names = ["Domain_id","gene_name","P-value","FDR","Bonferroni","A","B","observed_ratio","C","D","expected_ratio"]
    print(result)
    result_return = result.to_dict('records')
    return {"result": result_return}

if __name__ == "__main__":
    seq_data = "TMEM158,CTAGE6,DCP1A,MARCKS"
    correction = "FDR"
    p_limit = 0.01
    data = enrichment(seq_data, correction, p_limit)
    print(data)