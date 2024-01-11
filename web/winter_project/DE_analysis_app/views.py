from django.shortcuts import render
import pandas as pd

# Create your views here.
def DE_page(request):
    return render(request, 'DE_analysis.html', locals())

def DE_cal(request):
    DE_level = request.POST["type"]
    cancer = request.POST["cancer"]
    primary_site, project = cancer.split("|")
    DE_filter = request.POST.getlist("DE_filter[]")
    print(DE_filter)
    ########## need to check condition1 and condition2 column data ##########
    table_data, table_col = DEscreener_getdata(DE_level, primary_site, project, DE_filter)
    df_DE = pd.DataFrame(table_data, columns=table_col)
    print(table_col)
    # df_DE = pd.read_csv(f"{current_path}/../static/data/DE_data/{selected_condition}.csv", skiprows=8)
    # print(df_DE)
    df_DE[list(df_DE.columns)[-2]] = df_DE[list(df_DE.columns)[-2]].apply(lambda x: round(1/x, 5))
    replace_col_name = list(df_DE.columns)[-1]
    print(list(df_DE.columns)[-1])
    df_DE.rename(columns={replace_col_name: "q-value"}, inplace=True)
    # print(df_DE)
    lefton = "name" if "name" in df.columns else "gene_name"
    df = df_DE if len(df) == 0 else pd.merge(df, df_DE, left_on=lefton, right_on="gene_name", how="inner")