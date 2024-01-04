function Tag_len_to_adjust_set() {
    var previousLength = 0;
    // Set up an interval to check the length at regular intervals (every 1000 milliseconds)
    var intervalId = setInterval(function() {
        var currentLength = $('#miRNA_input_area').jsonTagEditor('getTags')[0].tags.length;
        // Check if the length has changed
        if (currentLength !== previousLength) {
            if (currentLength < 2) {
                $("#set_operation_area").css("display", "none"); 
            } else {
                $("#set_operation_area").css("display", "block");
                if(currentLength != 2){
                    // alert('Please input 2 miRNA');
                    // c2_select.find("option[value='" + selectedValueC1 + "']").prop("disabled", true);
                    $(".diff_set").hide();
                    // $("#miRNA_difference2_1").hide();
                    $("#miRNA_union").prop("checked", true);
                }
                else{
                    $(".diff_set").show();
                    // $("#miRNA_difference1_2").show();
                    // $("#miRNA_difference2_1").show();
                }
            }
            // Update the previous length for the next check
            previousLength = currentLength;
        }
    }, 200); // Adjust the interval time as needed (e.g., every second)
}

$(document).ready(function(){
    Tag_len_to_adjust_set();
    // for miRNA screener autocomplete data
    var homo_miRNA_list = document.getElementById('homo_miRNA_list').getAttribute('data-json').replace(/[\[\]\'()]/g, '').split(",");
    var miRNA_type = "homo";
    var miRNA_list;
    // miRNA type switch
    switch (miRNA_type) {
        case 'homo':
            miRNA_list = homo_miRNA_list;
            break;
        default:
            alert('沒有符合的條件');
    }

    // miRNA screener setting
    $('#miRNA_input_area').jsonTagEditor({
        autocomplete: {
            minLength: 7,
            delay: 0, // show suggestions immediately
            position: { collision: 'flip' }, // automatic menu position up/down
            source: miRNA_list
        },
        forceLowercase: false,
        placeholder: 'miRNA name(s)'
    });

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
        //// miRNA screener
        var miRNA_switch = $("#switch_miRNA").prop('checked');
        $('#miRNA_input_area').jsonTagEditor('getTags')[0].tags.length;
        var selected_miRNA = [];
        for(var i = 0; i < $('#miRNA_input_area').jsonTagEditor('getTags')[0].tags.length; i++){
            selected_miRNA.push($('#miRNA_input_area').jsonTagEditor('getTags')[0].tags[i].value);
        }
        // selected_miRNA = getSelectedCheckboxes('#miRNA_checkbox');
        var miRNA_set = document.querySelector('input[name="miRNA_set"]:checked').value;

        $.ajax({
            headers: { 'X-CSRFToken': csrf_token },
            type: 'POST',
            // url:'/survival_analysis/cal_pvalue/',
            url:'/screener/screener_cal_result_gene/',
            dataType : 'json',
            data : {
                // miRNA
                'selected_miRNA[]':selected_miRNA,
                'miRNA_set':miRNA_set,
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
                    select_cancer[0], {'survival': true, 'miRNA': false, 'DE': false}, selected_miRNA ,DE_info_dict);
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