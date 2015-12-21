function resetDisplay() {
    $('#matrix').empty();
    $('#reads_canvas').empty();
    $('#thresholdDiv').css('display', 'block');
}

function updateCell(geneID, tumors_error_list, threshold) {
    console.log("[LOG] updateCell received geneID: " + geneID + ", size(list): " + tumors_error_list.length);
    //console.log("[REMOVE] List: " + tumors_binary_error_list);
    var geneTagID = "gene_" + geneID;
    var geneRow = $("#" + geneTagID);

    for (var i = 0; i < tumors_error_list.length; i++) {
        if (tumors_error_list[i] > threshold)
            $('.cell', geneRow).eq(i).css('background-color', 'rgba(126, 0, 0, ' + (tumors_error_list[i] - threshold) + ')');
        else
            $('.cell', geneRow).eq(i).css('background-color', 'rgba(173, 229, 171, ' + (threshold - tumors_error_list[i]) + ')');
    }   
}

function removeGene(geneID, className){
    console.log("[LOG] removeGene received geneID: " + geneID);
    var geneTagID = "gene_" + geneID;
    $("#" + geneTagID).addClass(className).css('display', 'none');
}

function makeGeneMutated(geneID){
    console.log("[LOG] makeGeneMutated received geneID: " + geneID);
    var geneTagID = "gene_" + geneID;
    $("#" + geneTagID).css('display', 'none');
}

$(function(){
    $('#result').hide();
    $('#kill_ping').click(function() {
        $.ajax({
            url: "/kill_proc",
            cache: false,
            success: function(html){
                window.frames[0].stop();
                $("#result").html(html);
                $("#result").show('slow');
            }
        });
        return false;
    });
});

$(function(){
    $('#ping').click(function() {
        $('#result').hide();
    });
});

function createMatrix(genes_cnt, tumors_cnt) {
    $( document ).ready(function() {
        var nrow = genes_cnt;
        var ncol = tumors_cnt;

        // Header Row:
        var tumorHeaders = document.createElement("div");
        $(tumorHeaders).attr({
            'id': 'tumors_header_container'
        });
        var gh = document.createElement("div");
        $(gh).attr({
            'id': 'geneheader_placeholder' ,
            'class': 'gene_header'
        });
        $(tumorHeaders).append(gh);
        var header;

        for (var c = 0; c < ncol; c++){
            header = document.createElement("div");
            $(header).attr({
                'id': 'tumorheader_' + c,
                'class': 'tumor_header vertical_text'
            });
            header.innerHTML = "T" + c;
            $(tumorHeaders).append( header );
        }
        $( "#matrix" ).append( tumorHeaders );

        for (var r = 0; r < nrow; r++) {
            var row = document.createElement("div");
            $(row).attr({
                'id': 'gene_' + r,
                'class': 'geneRow'
            });




            // First, add the gene header
            header = document.createElement("div");
            $(header).attr({
                'id': 'geneheader_' + r,
                'class': 'gene_header'
            }).
            click(geneSelector);
            header.innerHTML = 'Gene ' + r;
            $(row).append(header);

            for (c = 0; c < ncol; c++){
                var col = document.createElement("div");
                $(col).attr({
                    'id': 'g_' + r + '_t_' + c,
                    'class': 'cell'
                }).
                click(cellClick);
                $(row).append( col );
            }
            //$(row).append('<div style="clear:both"></div>');
            $( "#matrix" ).append( row );
        }
        //console.log( "document loaded" );
    });
}


$(function() {
    // When the testform is submitted…
    $("#testform").submit(function() {
        // post the form values via AJAX…
        var postdata = {
            gene: $("#gene").val(),
            tumor: $("#tumor").val()
        } ;
        console.log('Post data ' + JSON.stringify(postdata))
        $.post('/showReads', postdata, function(data) {
            // the following function is called once the result is ready
            ;
        });
        return false ;
    });
});


function showReadsForGene(jsonObj) {
    var MutationStatusTranslation = ['undecided', 'mutated', 'normal'];

    // 0 --> undecided
    // 1 --> mutated
    // 2 --> normal
    var mutationStatus = 0;  // will be changed later
    var oppositeMutationStatus = 0;


    /* ********************************
    // First, clean up the previous reads 
    **********************************/
    $('#reads_canvas').empty();

    var baseWidth = 20;
    borderWidth = 1;    
    reference_gene = jsonObj.reference_gene;


    /* ********************************
        // Compute mutation rate per base 
    **********************************/
    var gene_length = reference_gene.gene.length
    var readCnt = new Array(gene_length);
    var mutationCnt = new Array(gene_length);
    var mutationRatio = new Array(gene_length);
    var maxMutationRatio = -1;

    // first initialize the arrays
    for (var i = 0; i < gene_length; i++) {
        readCnt[i] = 0;
        mutationCnt[i] = 0;
    }

    for (var i = 0; i < jsonObj.reads_list.length; i++){
        read = jsonObj.reads_list[i];
        read_relative_start_pos = read.start_position - reference_gene.gene_begin_position;

        for (var p = 0; p < read['read'].length; p++){
            if (p + read_relative_start_pos >= gene_length)
                break;

            readCnt[read_relative_start_pos + p]++;
            if (reference_gene.gene[read_relative_start_pos + p] != read['read'][p])
                mutationCnt[read_relative_start_pos + p]++;
        }
    }
    for (var i = 0; i < gene_length; i++){
        if (readCnt[i] == 0) mutationRatio[i] = 0;
        else mutationRatio[i] = mutationCnt[i] / readCnt[i];
        mutationRatio[i] = parseFloat(mutationRatio[i].toFixed(2));    // keep two decimal
    }

    for (var i = 0; i < gene_length; i++)
        if (mutationRatio[i] > maxMutationRatio)
            maxMutationRatio = mutationRatio[i];

    console.log('Threshold ' + jsonObj.threshold)

    if (maxMutationRatio > jsonObj.threshold) {
        mutationStatus = 1;
        oppositeMutationStatus = 2
    }
    else {
        mutationStatus = 2;
        oppositeMutationStatus = 1;
    }



    /* ********************************
    // Display stats
    **********************************/
    var buttons_stats = document.createElement("div");
    $(buttons_stats).attr({'id': 'buttons_stats' });


    var stats = document.createElement("div");
    $(stats).attr({'id': 'stats'});

    var stat_ul = document.createElement("ul");

    var gene_id_li = document.createElement("li");
    $(gene_id_li).html('<b>Gene ID:</b> ' + reference_gene.gene_id);

    var tumor_id_li = document.createElement("li");
    $(tumor_id_li).html('<b>Tumor ID:</b> ' + jsonObj.tumor_id);

    var max_mut_rate_li = document.createElement("li");
    $(max_mut_rate_li).html('<b>Max Mut. Ratio:</b> ' + maxMutationRatio);

    var mutation_status_li = document.createElement("li");
    $(mutation_status_li).html('<b>Mut. Status:</b> ' + MutationStatusTranslation[mutationStatus]);

    $(stat_ul).
    append(gene_id_li).
    append(tumor_id_li).
    append(max_mut_rate_li).
    append(mutation_status_li)

    $(stats).append(stat_ul);
    $(buttons_stats).append(stats)


    /* ********************************
    // Display fix buttons
    **********************************/
    var fixGeneButton = document.createElement("div");

    $(fixGeneButton).attr({
        'class': 'fix_gene myButton',
        'id': 'fix_gene_' + reference_gene.gene_id + '_tumor_' + jsonObj.tumor_id,
        'value': mutationStatus
    }).
    html('Treat this gene <b>' + MutationStatusTranslation[oppositeMutationStatus] + '</b>').
    click(function(){
        console.log($(this).text());
        $(this).addClass('not_active').unbind('click');
        var postdata = {
            gene: reference_gene.gene_id,
            tumor: jsonObj.tumor_id,
            newstate: oppositeMutationStatus
        } ;
        $.post('/fixGeneMutation', postdata, function(data) {
            // now change the color of the cell
            if (oppositeMutationStatus == 2)
                removeGene(reference_gene.gene_id, 'fixed');
            else makeGeneMutated(reference_gene.gene_id);
        });
    });
    $(buttons_stats).append(fixGeneButton);

    $('#reads_canvas').append(buttons_stats);



    /* ********************************
        // Display Reference Gene
    **********************************/
    var referenceGeneDiv = document.createElement("div");
    $(referenceGeneDiv).
    attr({'id': 'reference_gene'}).
    css('width', gene_length * baseWidth + 4);

    for (var i = 0; i < reference_gene.gene.length; i++){
        var referenceBaseDiv = document.createElement("div");
        var title =
            'Mut: ' + (mutationRatio[i] * 100) + '%' +
            ' (Pos: ' + (i + reference_gene.gene_begin_position) + 
            ', #R: ' + readCnt[i] +
            ', #M: ' + mutationCnt[i] + ')'


        $(referenceBaseDiv).
        attr({
            'class': 'reference_base_div',
            'id': 'reference_base_' + i,
            'title': title}).
        css('width', baseWidth).
        text(reference_gene.gene[i]).
        tooltip().
        hover(function() {$( this ).css('border', 'solid 1px');},
              function(){$( this ).css('border-style', 'none');});

        if (jsonObj.threshold > mutationRatio[i])
            $(referenceBaseDiv).addClass('below_threshold_base').
            css('background-color', 'rgba(6, 255, 6, ' + (0.5 - mutationRatio[i]) + ')');
        else $(referenceBaseDiv).addClass('above_threshold_base').
        css('background-color', 'rgba(255, 119, 69, ' + (mutationRatio[i]) + ')');


        $(referenceGeneDiv).append(referenceBaseDiv);
    }
    $('#reads_canvas').append(referenceGeneDiv);



    /* ********************************
    // Add an empty div
    **********************************/
    var emptyDiv = document.createElement("div");
    $(emptyDiv).css('height', 10);
    $('#reads_canvas').append(emptyDiv);


    /* ********************************
    // Display Reads
    **********************************/
    var allReadsDiv = document.createElement("div");
    $(allReadsDiv).attr({'id': 'reads'});    

    var level = 0;
    while (true) {
        // Find reads that belond to this level
        readsInThisLevel = jsonObj.reads_list.filter(function(readObj){
            return readObj.level  == level;
        });
        if (readsInThisLevel.length == 0)
            break;

        console.log('In level ' + level + ': found ' + readsInThisLevel.length + ' objects');

        var levelDiv = document.createElement("div");
        $(levelDiv).attr({
            'class': 'level_div',
            'id': 'level_' + level
        });

        // add the reads to their corresponding level div
        for (var i = 0; i < readsInThisLevel.length; i++){
            read = readsInThisLevel[i];
            read_relative_start_position = read.start_position - reference_gene.gene_begin_position

            if (i == 0){
                // first read in this level.
                // we add a div to fill in the gap up to the read

                var gapFiller = document.createElement("div");
                $(gapFiller)
                    .attr({'class': 'gap_filler'})
                    .css('width', read_relative_start_position * baseWidth - borderWidth);
                $(levelDiv).append(gapFiller);
            }
            else {
                // add a gap filler from the last read in that level up to this read
                var posDiff = read.start_position - (readsInThisLevel[i-1].start_position + readsInThisLevel[i-1].read.length);
                var gapFiller = document.createElement("div");
                $(gapFiller).attr({
                    'class': 'gap_filler'
                }).css('width', posDiff * baseWidth  - 2 * borderWidth);
                $(levelDiv).append(gapFiller);
            }


            // time to add the actual read
            var readDiv = document.createElement("div");
            for (var c = 0; c < read['read'].length; c++){
                var baseSpan = document.createElement("div");
                $(baseSpan)
                    .attr({
                    'class': 'read_base',
                    'id': 'level_' + level + '_read_' + i + '_base_' + c 
                }).css('width', baseWidth)
                    .text(read['read'][c]);
                if (c + read.start_position >= reference_gene.gene_begin_position + reference_gene.gene.length)
                    $(baseSpan).addClass('unimportant_base');
                else if (read['read'][c] == reference_gene.gene[read_relative_start_position + c])
                    $(baseSpan).addClass('normal_base')
                        .css('background-color', 'rgba(182, 255, 176, ' + read.quality + ')'); 
                else 
                    $(baseSpan).addClass('mutated_base')
                        .css('background-color', 'rgba(227, 86, 60, ' + read.quality + ')'); 

                $(readDiv).append(baseSpan);
            }

            $(readDiv)
                .attr({
                'class': 'read_string',
                'id': 'level_' + level + '_read_' + i,
                'title': 'Read quality: ' + read.quality
            }).
            css('width', (read['read'].length * baseWidth) ).
            css('overflow', 'hidden').
            tooltip().
            hover(function() {$( this ).css('border', 'solid 2px');},
                  function(){$( this ).css('border', 'solid 1px');});
            $(levelDiv).append(readDiv);
        }
        $(allReadsDiv).append(levelDiv);


        var emptyRow = document.createElement("div");
        $(emptyRow).css('height', 5);
        $(allReadsDiv).append(emptyRow);
        level ++;
    }
    $('#reads_canvas').append(allReadsDiv);
};


$(function(){
    $("#threshold").on("change", function(){
        $('#threshold_number').text((this.value / 100).toFixed(2))});

});

function hideGenes(){
    var postdata = {};
    $.post('/hideGenes', postdata, function(genesIDs) {
        // the following function is called once the result is ready
        console.log('received ' + genesIDs.length + ' to hide');
        for (var i = 0; i < genesIDs.length; i++)
            removeGene(genesIDs[i], 'hide');
    });
    return false ;

}

function showGenes(){
    $('.geneRow.hide').removeClass('hide').css('display', 'block' );
}

function sleep(milliseconds) {
    var start = new Date().getTime();
    for (var i = 0; i < 1e7; i++) {
        if ((new Date().getTime() - start) > milliseconds){
            break;
        }
    }
    alert("woke up!");
}

function cellClick(){
    // find gene and tumor id by splitting the element ID (g_{geneID}_t_{tumorID}), e.g. g_0_t_9
    var s = $(this).attr('id');
    myArr = s.split('_')

    console.log('Clicked on gene ' + myArr[1] + ' and tumor ' + myArr[3]);
    var postdata = {
        gene: myArr[1],
        tumor: myArr[3]
    };
    $.post('/showReads', postdata, function(json_obj) {
        showReadsForGene(json_obj);
    });
}

function geneSelector(){
    geneID = parseInt($(this).attr('id').split('_')[1]);
    console.log('geneSelector() called for gene ' + geneID);

    $('#correlationBox').empty().
    css('display', 'inline-block').
    html('<p>Gene <b>' + geneID + '</b> selected.</p> </br> Top (anti-)correlated genes:');

    var postdata = {gene: geneID};
    $.post('/getCorrelatedGenes', postdata, function(relatedGenesIDs) {
        $('#correlationBox').append('<ul>');
        text = '<ul>';

        // the following function is called once the result is ready
        console.log('received ' + relatedGenesIDs.length + ' as correlated genes');

        for (var i = 0; i < relatedGenesIDs.length; i++) {
            text += '<li><div id = "corr_1_1" onclick="goToCorrelation(' + relatedGenesIDs[i] + ')">Gene ' + relatedGenesIDs[i] +'</div></li>'
        }
        $('#correlationBox').append(text)
    });
}

function goToCorrelation(geneID){
    var contactTopPosition = $("#gene_" + geneID).position().top;
    $("#matrix").animate({scrollTop: contactTopPosition});

//    $('html, body').animate({
//        'scrollTop' : $("#gene_" + geneID).position().top
//    });
}