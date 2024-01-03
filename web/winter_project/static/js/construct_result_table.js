function construct_prediction_datatable(name, search_type, data, high_percentile, low_percentile, stage, cancer, switch_dict, selected_miRNA, DE_info_dict)
{
    //// Dynamically add/remove datatable column&data using ajax.
    // without if condition will raise error about reinitialise(Cannot read property 'aDataSort' of undefined)
    if ( $.fn.DataTable.isDataTable( `#${name}` ) ) {
        $(`#${name}`).DataTable().destroy();
        $(`#${name}`).empty();
    };
    // colspan column
    var headers = '<thead><tr><th colspan="1">Name</th>'
    // <th colspan="2">Title</th><th colspan="2">1</th><th colspan="2">2</th><th colspan="2">3</th><th colspan="2">4</th><th colspan="2">5</th><th colspan="2">6</th><th colspan="2">7</th><th colspan="2">8</th></tr>';
    // default columns
    var columns = [{ title:`${search_type}`, data: "name"}];
    var columnDefs = []
    var target_index = 1;
    // create survival analysis columns and columnDefs
    if(switch_dict.survival == true){
        headers += '<th colspan="2">survival analysis</th>';
        columns.push(
            { title:"survival analysis logrank p-value", data: "logrank_p_value"},
            { title:"survival analysis detail"},
        );
        columnDefs.push(
            {
                // 指定第一列，從0開始，0表示第一列，1表示第二列……
                targets: target_index ,
                autoWidth: true,
                render: function(data, type, row, meta) {
                    // console.log(typeof data);
                    if (type === 'display' && typeof data === 'number') {
                        return Number(data.toFixed(6)).toExponential();
                    }
                    return data;
                },
            },
            {
                // 指定第一列，從0開始，0表示第一列，1表示第二列……
                targets: target_index + 1,
                autoWidth: true,
                render: function(data, type, row, meta) {
                    // if id has space or : or . it can not be use
                    // var title = meta.settings.aoColumns[meta.col-2].sTitle;
                    // console.log(title)
                    // var name = row[title];
                    var name = row["name"]; // 後面的key是用data:後的名稱
                    return `<a href="/survival_analysis/detail/?cancer=${cancer}&type=${search_type}&name=${name}&hp=${high_percentile}&lp=${low_percentile}&stage=${stage}" target="_blank">view</a>`
                },
            },
        );
        target_index += 2;
    }
    // create miRNA columns and columnDefs
    if(switch_dict.miRNA == true){
        headers += `<th colspan="${selected_miRNA.length}">miRNA</th>`;
        selected_miRNA.forEach((miRNA) => {
            columns.push(
                { title:`${miRNA}`, data: `${miRNA}`}
            );
        });
        target_index += selected_miRNA.length;
    }
    // create DE columns and columnDefs
    if(switch_dict.DE == true){
        headers += `<th colspan="4">DE</th>`;
        columns.push(
            { title:`${DE_info_dict["condition1"]} Avg FPKM`, data: `value_1`},
            { title:`${DE_info_dict["condition2"]} Avg FPKM`, data: `value_2`},
            { title:`Fold Change`, data: `fold_change`},
            { title:`q-value`, data: `q-value`},
        );
        columnDefs.push(
            {
                // 指定第一列，從0開始，0表示第一列，1表示第二列……
                targets: target_index +3,
                autoWidth: true,
                render: function(data, type, row, meta) {
                    // console.log(typeof data);
                    if (type === 'display' && typeof data === 'number') {
                        return Number(data.toFixed(6)).toExponential();
                    }
                    return data;
                },
            },
        );
        target_index += 4;
    }
    // $(`#${name}`).append(headers);
    $(`#${name}`).DataTable({
        // 'scrollX':true, will split header and data, it will cause header and data misalign
        // datatable size will auto Resizes with the window
        "initComplete": function (settings, json) {  
            $(`#${name}`).wrap("<div style='overflow:auto; width:100%;position:relative;'></div>");            
        },
        // bSort: false,
        order: [[0, 'asc']],
        destroy : true,
        data: data,
        columns: columns,
        // [
        //     // "visible": false because these data are for download csv column
        //     { title:`${search_type}`, data: "name"},
        //     { title:"survival analysis logrank p-value", data: "logrank_p_value"},
        //     { title:"survival analysis detail"},
        //     // { title:"max_time", data:"max_time"},
        // ],
        columnDefs: columnDefs,
        // [
        //     {
        //         // 指定第一列，從0開始，0表示第一列，1表示第二列……
        //         targets: 2,
        //         autoWidth: true,
        //         render: function(data, type, row, meta) {
        //             // if id has space or : or . it can not be use
        //             // var title = meta.settings.aoColumns[meta.col-2].sTitle;
        //             // console.log(title)
        //             // var name = row[title];
        //             var name = row["name"]; // 後面的key是用data:後的名稱
        //             return `<a href="/survival_analysis/detail/?cancer=${cancer}&type=${search_type}&name=${name}&hp=${high_percentile}&lp=${low_percentile}&stage=${stage}" target="_blank">view</a>`
        //         },
        //     },
        //     {
        //         // 指定第一列，從0開始，0表示第一列，1表示第二列……
        //         targets: 1,
        //         autoWidth: true,
        //         render: function(data, type, row, meta) {
        //             // console.log(typeof data);
        //             if (type === 'display' && typeof data === 'number') {
        //                 return Number(data.toFixed(6)).toExponential();
        //             }
        //             return data;
        //         },
        //     },
        // ],
        dom: 'Bfrtip',
        buttons: [{
            extend: 'csv',
            text: 'Export CSV',
            title: `${cancer}_${search_type}_${stage}_${high_percentile}_${low_percentile}`,
            exportOptions: {
                columns: [0, 1]
                // columns: [0]
            }
        }],
        
    });
}
