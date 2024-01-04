$(document).ready(function(){
    $('#gene_search_btn').click(function(){
        var search_value = $("#search-bar").val(); // Get the input value from the search-bar element
        console.log(search_value);
        $.ajax({
            headers: { 'X-CSRFToken': csrf_token },
            type: 'POST',
            url:'/gene_search/search_ajax/',
            dataType : 'json',
            data : {
                'search_val':search_value,
            },
            success: function(response){
                $(`#gene_table`).DataTable({
                    "initComplete": function (settings, json) {  
                        $(`#gene_table`).wrap("<div style='overflow:auto; width:100%;position:relative;'></div>");            
                    },
                    destroy : true,
                    data: response.gene_info,
                    columns:[
                        { title:`Symbol`, data: "Symbol"},
                        { title:"mRNA transcript(s)", data: "mRNA transcript(s)"},
                        { title:"Synonyms", data: "Synonyms"},
                        { title:"Chromosome", data: "Chromosome"},
                        { title:"Description", data: "Description"},
                        { title:"Other Designations", data: "Other Designations"},
                    ],
                });
                $(`#transcript_table`).DataTable({
                    "initComplete": function (settings, json) {  
                        $(`#transcript_table`).wrap("<div style='overflow:auto; width:100%;position:relative;'></div>");            
                    },
                    destroy : true,
                    data: response.transcript_info,
                    columns:[
                        { title:`transcript_name`, data: "transcript_name"},
                        { title:"definition", data: "definition"},
                        { title:"transcript_variant", data: "transcript_variant"},
                    ],
                });
                $("#result_area").show();
            },
            error: function(xhr, ajaxOptions, thrownError){
                console.log(data);
                console.log(xhr.status);
                console.log(thrownError);
            }
        });
    });
});