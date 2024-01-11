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