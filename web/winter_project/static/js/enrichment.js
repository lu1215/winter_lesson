$.ajaxSetup({
	headers: { 'X-CSRFToken': csrf_token },
	type: 'POST',
});

$(document).ready(function(){
	$("#comparison_table").hide();
    $("#btn_cal").click(function(){
        var seq = $("#input_seq_content").val();
        var Correction = $("#c_type").val();
        var p_limit = $("#p-value").val();
		var c_type = $("#c_type").val();

        $.ajax({
			headers: { 'X-CSRFToken': csrf_token },
			url: '/enrichment_app/enrichment_ajax/',
			type: 'POST',
			dataType: "json",
			data: {
				seq:seq,
                Correction:Correction,
                p_limit:p_limit,
			},
			// async: false, 
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
			success: function(response){
                clearInterval(tID);
                delete tID
                swal.close();
                data = response.enrichment_result
				$("#relate_table").DataTable({
					destroy : true,
					data: data,
					columns:[
						{ title: 'Domain_id', data: 'Domain_id' },
						// { title: 'P-value', data: 'P-value' },
						{ title: '-log10(corrected P-value)' },
						{ title: 'FDR' , data: 'FDR'},
						{ title: 'Bonferroni', data: 'Bonferroni'},
						{ title: 'Expected Ratio'},
						{ title: 'Observed Ratio'},
						{ title: 'Fold Enrichment'},
						
						// { title: 'Comparsion'},
					],
					columnDefs: [
						// {
						// 	// 指定第一列，從0開始，0表示第一列，1表示第二列……
						// 	targets: 6,
						// 	render: function(data, type, row, meta) {
						// 		return `<button class="btn_comparison">${row[0]}</button>`
						// 	},
						// },
						{
							// 指定第一列，從0開始，0表示第一列，1表示第二列……
							targets: 1,
							render: function(data, type, row, meta) {
								return -Math.log10(row['P-value']);
							},
						},
                        {
							// 指定第一列，從0開始，0表示第一列，1表示第二列……
							targets: 4,
							render: function(data, type, row, meta) {
								return `${row['B']} / ${row['D']} = ${row['expected_ratio'].toFixed(4)}`
							},
						},
                        {
							// 指定第一列，從0開始，0表示第一列，1表示第二列……
							targets: 5,
							render: function(data, type, row, meta) {
								return `${row['A']} / ${row['B']} = ${row['observed_ratio'].toFixed(4)}`
							},
						},
						{
							// 指定第一列，從0開始，0表示第一列，1表示第二列……
							targets: 6,
							render: function(data, type, row, meta) {
								var value = row['observed_ratio'] / row['expected_ratio']
								return Number(value.toFixed(4));
								// return `${}`
							},
						},
                        // `${dataset[i][6]}/${dataset[i][9]} = ${Math.round((dataset[i][10])*100 *100)/100}%`,
					],
				});

				//只出現自己選擇的過濾評分項
				if(c_type == "FDR"){
					$("#relate_table").DataTable().column(1).visible(false);
					$("#relate_table").DataTable().column(2).visible(true);
					$("#relate_table").DataTable().column(3).visible(false);
				}
				else if(c_type == "Bonferroni"){
					$("#relate_table").DataTable().column(1).visible(false);
					$("#relate_table").DataTable().column(2).visible(false);
					$("#relate_table").DataTable().column(3).visible(true);
				}
				else{
					$("#relate_table").DataTable().column(1).visible(true);
					$("#relate_table").DataTable().column(2).visible(false);
					$("#relate_table").DataTable().column(3).visible(false);
				}
            },
            error: function(){
				alert('enrichment data error');
			},
        });
    });
});