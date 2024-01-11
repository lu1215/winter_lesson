def survival_plot(T1,E1,T2,E2,GT_input,primary_site,random_id,Low_Percentile,High_Percentile,survival_days,survival_select):
    plot_path = f"{current_path}/survival_plot_%s_%s.png"%(primary_site.replace(' ','_').replace("(",'').replace(")",''), GT_input)
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