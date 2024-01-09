function download_csv(){
    window.location.href = `{% static 'cancer_project/static/data/csv_result/Survival_Profile_${primary_site.tostring()}_${High_Percentile.tostring()}_${Low_Percentile.tostring()}.csv' %}`;
}
$(document).ready(function(){
    var project = $("#project").text().trim();
    var primary_site = $("#primary_site").text().trim();
    var input_type = $("#input_type").text().trim();
    var Low_Percentile = $("#Low_Percentile").text().trim();
    var High_Percentile = $("#High_Percentile").text().trim();
    var stage = $("#stage").text().trim();
    var name = $("#name").text().trim();
    var days = $("#days").text().trim();
    var dynamicURL = `../../static/data/csv_result/Survival_Profile_${primary_site}_${name}_${High_Percentile}_${Low_Percentile}.csv`;
    document.getElementById('download_link').href = dynamicURL;
    // {% static 'cancer_project/static/data/csv_result/Survival_Profile_${primary_site.tostring()}_${High_Percentile.tostring()}_${Low_Percentile.tostring()}.csv' %}
    // if('{{survival_str|safe}}' == ''){
        // var max_survival_days_realtime = {{ max_survival_times|safe }}[0];
        // var random_id = "{{ random_id }}";
        // var now = new Date().getTime();
        // var primary_project_stage = {{ primary_project_stage|safe }}
        // // console.log({{ boxplot_data|safe }})
        // // console.log({{ primary_site_list|safe }})
        // // console.log({{ boxplot_data|safe }})
        // // console.log(primary_project_stage)
        // // console.log({{ project_list|safe }})
        // // console.log(max_survival_days_realtime)
        // var project = {{ project_list|safe }}[0]
        // var primary_site = {{ primary_site_list|safe }}[0]
        // var max_survival_times = {{ max_survival_times|safe }}[0]
        // var stages_list = primary_project_stage[primary_site+"|"+project]

        // condition_list = condition_list.map(x => x.charAt(0).toUpperCase() + x.slice(1));
    // $('#boxplot').append('<div id="'+project+'"></div><br><br>');    
    // var survival_toolbar = '<br><br>'+
    //                         '<h3 style="width:100%;height: 41.2px;line-height: 41.2px;">Survival Plot of '+"<a href='{{ genelink_dict|get_item:GT_input }}' target='_blank'>{{ GT_input }}</a>"+' in <font style="color:red">'+primary_site+'</font> cancer</h3>'+
    //                         '<hr style="width: 100%; color: black; height: 1px; background-color:#4d4d4d;"/>';
    
    // // $(`#${project}`).append(survival_toolbar); 

    // var survival_select_option = "";
    // for (var i = 0; i < stages_list.length; i++) {
    //     if(stages_list[i] != 'normal'){
    //         survival_select_option += '<option value="'+stages_list[i]+'">Only '+stages_list[i]+'</option>'
    //     }
    // }
    // var primary_site = primary_site.replace("(",'').replace(")",'');
    // var primary_site_img = primary_site.replace(/\s+/g, "_");
    // var survival_plot_path = "/edward_static/Survival_plot/survival_plot_"+primary_site_img+"_"+random_id+".png?"+now;
    // $('#'+project).append('<div style="border:1px solid; border-radius:4px" align="center"><div class="input-group" style="width: 85%;z-index:0;height: 42px;margin-top: 15px;margin-left: 5%;">'+
    //                                 '<span class="input-group-addon" style="border-right:0;width:25%;padding: 2px 8px;">'+
    //                                     '<span style="font-weight: bold">Lower Percentile</span><hr style="width: 100%; color: black; height: 2px; background-color:#8B0000;" />'+
    //                                     '<input type="text" id="survival_input_low_'+primary_site+'" name="input_low" class="form-control" placeholder="Low Percentile" value="50" style="width:100%;height: 22px;">'+
    //                                 '</span>'+
    //                                 '<span class="input-group-addon" style="border-right:0;width:25%;padding: 2px 8px;">'+
    //                                     '<span style="font-weight: bold">High Percentile</span><hr style="width: 100%; color: black; height: 2px; background-color:#8B0000;" />'+
    //                                     '<input type="text" id="survival_input_high_'+primary_site+'" name="input_high" class="form-control" placeholder="High Percentile" value="50" style="width:100%;height: 22px;">'+
    //                                 '</span>'+
    //                                 '<span class="input-group-addon" style="width:25%;padding: 2px 8px;">'+
    //                                     '<span style="font-weight: bold">Days</span><hr style="width: 100%; color: black; height: 2px; background-color:#8B0000;" />'+
    //                                     '<input type="text" id="survival_days_'+primary_site+'" name="input_days" class="form-control" placeholder="Days" value='+max_survival_times+' style="width:100%;height: 22px;">'+
    //                                 '</span>'+
    //                                 '<span class="input-group-addon" style="border-left:0;padding: 2px 8px;width: 25%">'+
    //                                     '<span style="font-weight: bold">Samples Included</span><hr style="width: 100%; color: black; height: 2px; background-color:#8B0000;" />'+
    //                                     '<select id="survival_select_'+primary_site+'" style="height: 22px;width: 100%;" name="samples_select">'+
    //                                         '<option value="all stage">All Cancer Stages</option>'+
    //                                         survival_select_option+
    //                                     '</select>'+
    //                                 '</span>'+
    //                                 '<span class="input-group-btn">'+
    //                                     '<button id="'+primary_site+'" name="survival_submit" value="'+project+'|'+"{{ search_by }}"+'|'+"{{ GT_input }}"+'|'+random_id+'" class="btn btn-success" style="height: 100%;width: 100%;border-bottom-left-radius: 0;border-top-left-radius: 0;">Submit</button>'+
    //                                 '</span>'+
    //                                 '<div class="tooltip2" style="float: right; margin-right: 5%; margin-top: 30%;margin-left: 10px;"><span class="glyphicon glyphicon-question-sign" aria-hidden="true" style="font-size: 30px;"></span><span class="tooltiptext2">Patients are sorted by expression of your gene/transcript. Lower percentile refers to the lower slice you want. Upper percentile is the upper slice. I would recommend comparing bottom third versus top third, or bottom quartile versus top quartile, so 33:33, or 25:25. If there is a small number of patients 50:50 is likely the best option.</span></div>'+
    //                             '</div>'+
    //                             '<img id="img_'+primary_site_img+'" src="'+survival_plot_path+'" class="img-responsive" style="width: 90%;"><br><a href="' + "{% url 'download_realtime' %}?filename=Survival_Profile_"+primary_site+"_50_50_"+"{{ random_id }}"+'.csv" id="'+primary_site+'_download">'+
    //                                 '<button type="button" class="btn btn-danger" style="margin-bottom: 15px;"><span class="glyphicon glyphicon-download-alt" aria-hidden="true"></span> Download Data</button></a></div><br><br>');
    
    // $(".main-svg").css({"border": "1px solid","border-radius":"4px"})
    

    // $("[name='samples_select']").change(function(){

        
    //     var project = document.getElementById(primary_site).value.split('|')[0];
    //     var search_by = document.getElementById(primary_site).value.split('|')[1];
    //     var GT_input = document.getElementById(primary_site).value.split('|')[2];
    //     var survival_select = document.getElementById("survival_select_"+primary_site).value;
    //     $.ajax({
    //             url:"{% url 'survival_max_days' %}",
    //             data :  {
    //                 'project': project,
    //                 'search_by': search_by,
    //                 'GT_input': GT_input, 
    //                 'survival_select':survival_select,

    //                 csrfmiddlewaretoken: '{{ csrf_token }}'
    //             },
    //             type:'POST',

    //             datatype : 'text',

                                    
    //             success:function(data){
    //                 max_survival_days_realtime = data.max_survival_days
    //                 $("[name='input_days']").val(max_survival_days_realtime)
    //             },
                                            
    //             beforeSend:function(){

    //             },

    //             complete:function(){
    //                 $.unblockUI();
    //             },        
    //             error:function(xhr, ajaxOptions, thrownError){ 
    //                 alert(xhr.status); 
    //                 alert(thrownError);
    //             }    
    //         });
    // })

    $('#download_btn').click(function(){

    });



    $('#submit_plot').click(function(){
        // alert("submit_plot clicked!");
        project = $("#project").text().trim();
        primary_site = $("#primary_site").text().trim();
        input_type = $("#input_type").text().trim();
        Low_Percentile = $("#low_percent").val();
        High_Percentile = $("#high_percent").val();
        // var stage_ele = $("#stage_select options:selected");
        // console.log(stage_ele);

        stage = document.getElementById("stage_select").value;
        // var stage = stage_ele.value;
        name = $("#name").text().trim();
        // days = $("#days").val();
        
        var dynamicURL = `../../static/data/csv_result/Survival_Profile_${primary_site}_${name}_${High_Percentile}_${Low_Percentile}.csv`;
        document.getElementById('download_link').href = dynamicURL;
        // var select_elements = document.querySelectorAll('[id="select_element"]');
        // for(var i=0; i<select_elements.length; i++) {
        //     select_value.push(select_elements[i].value);
        //     console.log(i,select_elements[i].value);
        // }
        // console.log(select_value[0]);

        console.log(stage);
        console.log(days);
        console.log(Low_Percentile);
        console.log(High_Percentile);
        $.ajax({
            headers: { 'X-CSRFToken': csrf_token },
            type: 'POST',
            url:'/survival_analysis/sur_plot/',
            dataType : 'json',
            data : {
                'project':project, 
                'primary_site': primary_site,
                'type':input_type,
                'name':name,
                'stage': stage,
                'High_Percentile': High_Percentile,
                'Low_Percentile': Low_Percentile,
                // 'days': days,
            },
            beforeSend:function(){
                var count=0
                tID= setInterval(timedCount , 50);
                    function timedCount() {
                    count=count+0.05;
                    swal({
                        title: "Running...",
                        text: "It may take several minutes.\nPlease be patient.\n \nRunning time: "+parseInt(count)+" seconds\nClick anywhere of the page \nif the running time does not change",                       
                        button: false,
                    });
                };
            },           
            success:function(response){
                console.log(response.img_str);
                $("#sur_img").attr("src", '');
                $("#sur_img").attr("src", `data:image/png;base64,${response.img_str}`);
                clearInterval(tID);
                delete tID
                swal.close();
                // $('#gap').html('<br><br><br>')
                // $('#output').html(data);
                console.log(response.img_str);
                // $("#sur_img").attr("src", `data:image/png;base64,${response.img_str}`);
            },
            error:function(xhr, ajaxOptions, thrownError){ 
                clearInterval(tID);
                delete tID
                swal.close();
                swal('error')
            }       
        });                                                                                                      
    });

    // $('[name="survival_submit"]').click(function(){

    //     var primary_site = event.target.id
    //     var primary_site = primary_site.replace("(",'').replace(")",'');
    //     // var primary_site  = primary_site0.replace("(",'[').replace(")",']');
    //     console.log(primary_site);
    //     var project = document.getElementById(primary_site).value.split('|')[0];
    //     var search_by = document.getElementById(primary_site).value.split('|')[1];
    //     var GT_input = document.getElementById(primary_site).value.split('|')[2];
    //     var random_id = document.getElementById(primary_site).value.split('|')[3];
    //     var Low_Percentile = parseFloat(document.getElementById("survival_input_low_"+primary_site).value);
    //     var High_Percentile = parseFloat(document.getElementById("survival_input_high_"+primary_site).value);
    //     var survival_days = parseFloat(document.getElementById("survival_days_"+primary_site).value);
    //     var survival_select = document.getElementById("survival_select_"+primary_site).value;
        
    //     if(Low_Percentile >= 0 && Low_Percentile <= 100 && High_Percentile >= 0 && High_Percentile <= 100 && Low_Percentile+High_Percentile <= 100 ){
    //         if(survival_days > 0 && survival_days <= max_survival_days_realtime){
    //             $('#'+primary_site).css('cursor','wait');
    //             $('html').css('cursor','wait');

    //             $.ajax({
    //                 url:"{% url 'survival_plot_realtime' %}",
    //                 data :  {

    //                     'primary_site':primary_site,
    //                     'project':project,
    //                     'search_by':search_by,
    //                     'GT_input':GT_input,
    //                     'random_id':random_id,
    //                     'Low_Percentile':Low_Percentile,
    //                     'High_Percentile':High_Percentile,
    //                     'survival_days':survival_days,
    //                     'survival_select':survival_select,

    //                     csrfmiddlewaretoken: '{{ csrf_token }}'
    //                 },
    //                 type:'POST',

    //                 datatype : 'text',

                                        
    //                 success:function(data){
    //                     now = new Date().getTime();
    //                     var primary_site_img = primary_site.replace(/\s+/g, "_");
    //                     var survival_plot_path = "/edward_static/Survival_plot/survival_plot_"+primary_site_img+"_"+random_id+".png?"+now;
    //                     $('#img_'+primary_site_img).fadeOut('200', function() {
    //                         document.getElementById('img_'+primary_site_img).src = survival_plot_path;
    //                         $('#img_'+primary_site_img).fadeIn('200');
    //                         $('#'+primary_site).css('cursor','pointer');
    //                         $('html').css('cursor','default');
    //                     })
    //                     $('#'+primary_site+"_download").attr("href", "{% url 'download_realtime' %}?filename=Survival_Profile_"+Low_Percentile+"_"+High_Percentile+"_"+"{{ random_id }}"+".csv")
    //                 },
                                                
    //                 beforeSend:function(){

    //                 },

    //                 complete:function(){
    //                     $.unblockUI();
    //                 },        
    //                 error:function(xhr, ajaxOptions, thrownError){ 
    //                     alert(xhr.status); 
    //                     alert(thrownError);
    //                 }    
    //             });   
    //         }
    //         else{  
    //             if(survival_days < 0){
    //                 swal('Days must be greater than 0!')
    //             }
    //             else if(survival_days >= max_survival_days_realtime){
    //                 swal('Days must be less than '+max_survival_days_realtime+'!')
    //             }
    //         }
    //     }
    //     else{
    //         swal('Invalid percentile!')
    //     } 
    // });
// }
});