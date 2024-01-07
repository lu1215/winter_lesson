$(document).ready(function(){
    $('#search_btn').click(function(){
        var search_type = document.querySelector('input[name="type"]:checked').value;
        var select_elements = document.querySelectorAll('[name="select_element"]');
        var select_cancer = [];
        for(var i=0; i<select_elements.length; i++) {
            select_cancer.push(select_elements[i].value);
            console.log(i,select_elements[i].value);
        }
        //// p-value screener
        // get stage
        var stage = document.getElementById("stage_select").value;
        // get percentile
        var high_percent = document.getElementById("high_percentile").value;
        var low_percent = document.getElementById("low_percentile").value;
        // get input pvalue
        var p_value = document.getElementById("p_value").value; 

        $.ajax({
            headers: { 'X-CSRFToken': csrf_token },
            type: 'POST',
            // url:'/survival_analysis/cal_pvalue/',
            url:'/survival_analysis/search_ajax/',
            dataType : 'json',
            data : {
                'type':"genes",
                'cancer':select_cancer[0],
                // p value
                'stage': stage,
                'high_percent': high_percent,
                'low_percent': low_percent,
                "pvalue":p_value,
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
                // $(`#output_table`).html();
                $('#gap').html('<br><br><br>')
                // $('#output').html(data);
                // $('html,body').animate({scrollTop:$('#output').offset().top},800);
                console.log(response.result);
                construct_prediction_datatable(
                    "output_table", search_type, response.result, 
                    high_percent, low_percent, stage, 
                    select_cancer[0], {'survival': true, 'miRNA': false, 'DE': false}, selected_miRNA ,{});
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
    });
});