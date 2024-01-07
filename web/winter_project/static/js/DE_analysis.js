$(document).ready(function(){
    $('#search_btn').click(function(){
        // get input type
        var search_type = document.querySelector('input[name="type"]:checked').value;
        // get input cancer
        var select_elements = document.querySelectorAll('[name="select_element"]');
        var select_cancer = [];
        for(var i=0; i<select_elements.length; i++) {
            select_cancer.push(select_elements[i].value);
            console.log(i,select_elements[i].value);
        }
        //// DE screener
        var DE_switch = $('#switch_DE').prop('checked');
        // var DE_condition = document.querySelector('input[name="DE_tmp"]:checked').value;
        var DE_filter_elements = document.querySelectorAll('[name="DE_filter_elements"]');
        var DE_filter = [];
        for(var i=0; i<DE_filter_elements.length; i++) {
            DE_filter.push(DE_filter_elements[i].value);
            console.log(i,DE_filter_elements[i].value);
        }
        var DE_info_dict = new Object();
        DE_info_dict.condition1 = document.getElementById("condition1_pre").value.split('|')[0];
        DE_info_dict.condition2 = document.getElementById("condition2_pre").value.split('|')[0];
        console.log(DE_info_dict.condition1, DE_info_dict.condition2);
        var switch_dict = new Object();
        console.log(survival_switch, miRNA_switch, DE_switch);
        switch_dict.survival = false;
        switch_dict.miRNA = false;
        switch_dict.DE = true;
        var switch_string = JSON.stringify(switch_dict)
        if( survival_switch==false && miRNA_switch==false && DE_switch==false){
            swal('Please open at least 1 screener');
        }
        else if(select_cancer.includes('') == false ){
            $.ajax({
                headers: { 'X-CSRFToken': csrf_token },
                type: 'POST',
                // url:'/survival_analysis/cal_pvalue/',
                url:'/screener/screener_cal_result_gene/',
                dataType : 'json',
                data : {
                    'switch_dict':switch_string,
                    'type':search_type,
                    'cancer':select_cancer[0],
                    // p value
                    'stage': stage,
                    'high_percent': high_percent,
                    'low_percent': low_percent,
                    "pvalue":p_value,
                    // miRNA
                    'selected_miRNA[]':selected_miRNA,
                    'miRNA_set':miRNA_set,
                    // DE
                    'DE_filter[]':DE_filter,
                    // 'DE_condition':DE_condition,
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
                    console.log(response)
                    clearInterval(tID);
                    delete tID
                    swal.close();
                    $(`#output_table`).html();
                    $('#gap').html('<br><br><br>')
                    // $('#output').html(data);
                    // $('html,body').animate({scrollTop:$('#output').offset().top},800);
                    console.log(response.result);
                    construct_prediction_datatable(
                        "output_table", search_type, response.result, 
                        high_percent, low_percent, stage, 
                        select_cancer[0], switch_dict, selected_miRNA ,DE_info_dict);
                    $('#input_type_td').text(search_type);
                    $('#primary_site_td').text(select_cancer[0]);
                    $('#high_percent_td').text(high_percent + "%");
                    $('#low_percent_td').text(low_percent + "%");
                    $('#input_pvalue_td').text(p_value);
                    $('#output_message').text(`${response.result.length} ${search_type} met the input criteria`);
                    $('#screener_type_td').text(response.screener_type);
                    var resultElement = document.getElementById('result_block');
                    resultElement.style.display = 'block';
                },
                error:function(xhr, ajaxOptions, thrownError){ 
                    alert(thrownError);
                    clearInterval(tID);
                    delete tID
                    swal.close();
                    swal('error')
                }       
            });    
        }
        else{
            swal('Input Error!')
        }                                                                                                       
    });
});