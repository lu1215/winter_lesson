//// add condition1 and condition2 select option
function addoption(response, select_item){
    select_item.empty();
    var default_option = $("<option>").val("default").text("Select an option");
    select_item.append(default_option);
    for (var i = 0; i < response.condition_list.length; i++) {
        var condition = response.condition_list[i];
        // Accessing the number_of_sample using the condition.1 key
        var numberOfSample = response.number_of_sample[condition[1]];
        // Check if the condition is met
        if (numberOfSample !== '1') {
            // Create a new option element
            var option = document.createElement("option");
            // Set the value attribute
            option.value = condition[1] + "|" + condition[0] + "|" + numberOfSample;
            // Set the inner text (displayed text)
            option.innerText = condition[1].charAt(0).toUpperCase() + condition[1].slice(1) + " (n=" + numberOfSample + ")";
            // Append the option to the select element
            select_item.append(option);
        }
    }
    
}

//// Statistical Significance options setting
function select_function(){
	var select_test = document.getElementById("TEST_select").value;
	var test_state = document.getElementById("TESTstates_select");
	if(select_test=='Cuffdiff DE test'){
		test_state[0].style.display="none";
		test_state[1].style.display="none";
		test_state[2].style.removeProperty('display');
		test_state[2].selected = 'selected';
	}
	else{
		test_state[0].selected = 'selected';
		test_state[0].style.removeProperty('display');
		test_state[1].style.removeProperty('display');
		test_state[2].style.display="none";
	}
}

$(document).ready(function(){
    //// construct the condition1 and condition2 select option
    $('#primary_site_div').change(function(){
        // $('#c1_div').html("<select id='condition1_pre' style='width: 100%' disabled><option></option></select>"); 
        var project_select = document.getElementById("primary_select").value.split('|')[1];
        $.ajax({
            headers: { 'X-CSRFToken': csrf_token },
            type: 'POST',
            url:'/screener/primary_site_realtime/',
            dataType : 'json',
            data :{
                'project_select':project_select,
            },
            success:function(response){
                $('#condition1_pre').removeAttr("disabled");
                $('#condition2_pre').removeAttr("disabled");
                addoption(response, $('#condition1_pre'));
                addoption(response, $('#condition2_pre'));
            },
            beforeSend:function(){},
            complete:function(){},        
            error:function(xhr, ajaxOptions, thrownError){ 
                alert(xhr.status); 
                alert(thrownError);
            }       
        });  
    });
    
    //// avoid user select same condition
    // use $(document).on("change", "#c1_div", function () {}) to avoid iwe cant listen element changed by ajax
    $(document).on("change", "#c1_div", function () {
        var selectedValueC1 = $("#condition1_pre").val();
        var c2_select = $("#condition2_pre");
        c2_select.find("option[value='" + selectedValueC1 + "']").prop("disabled", true);
        c2_select.find("option").not("option[value='" + selectedValueC1 + "']").prop("disabled", false);
        // $(this).find("option[value='" + "default" + "']").prop("disabled", true);
    });

    $(document).on("change", "#c2_div", function () {
        var selectedValueC2 = $("#condition2_pre").val();
        var c1_select = $("#condition1_pre");
        c1_select.find("option[value='" + selectedValueC2 + "']").prop("disabled", true);
        c1_select.find("option").not("option[value='" + selectedValueC2 + "']").prop("disabled", false);
        // $(this).find("option[value='" + "default" + "']").prop("disabled", true);
    });

    //// swap two select value
    document.getElementById('reverseButton').addEventListener('click', function() {
        var select1 = document.getElementById('condition1_pre');
        var select2 = document.getElementById('condition2_pre');
        // Swap the values of the select elements
        var temp = select1.value;
        select1.value = select2.value;
        select2.value = temp;
    });
});
