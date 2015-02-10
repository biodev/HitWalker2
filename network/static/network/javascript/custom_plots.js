
jsplot_functions = {
    'waterfall_plot':draw_waterfall_plot
}

function make_waterfall_plot(obj)
{
    var clicked_gene = $(obj).attr("data-gene");
    var post_data = {type:$(obj).attr("data-type"), sample:$(obj).html()};
    
    //first look to see if the plot exists
    
    var exist_plot = d3.select('[name='+post_data.type+'_'+post_data.sample+']');
    
    ///if so then simply highlight the appropriate gene
    if (exist_plot[0][0] != null)
    {
       
        exist_plot.selectAll("g .x.axis .tick.major text")
                            .classed("highlight", function(d,i)
                            {
                               
                                var temp_gene = d3.set(d3.select(this).attr("data-genel").split(','));
                                //console.log(JSON.stringify(temp_gene));
                                if (temp_gene.has(clicked_gene) || d3.select(this).attr("class") == "plot text highlight")
                                {
                                    return(true);
                                }
                                else
                                {
                                    return(false);
                                }
                            });
                            
        exist_plot.selectAll("rect.bar")
                            .classed("highlight", function(d,i)
                                  {
                                    //console.log(JSON.stringify(d));
                                    var temp_gene = d3.set(d.db_name);
                                     if (temp_gene.has(clicked_gene) || d3.select(this).attr("class") == "bar highlight")
                                     {
                                        return(true);
                                     }
                                     else
                                     {
                                        return(false);
                                     }
                                  });
                            
        //clicked_gene needs to be added to the datasource so that it can be remembered on refresh
        
        //not going to add this to history for now as it doesn't really impact node content, simply highlights
        //var use_panel = panel_num_from_id (exist_plot.attr("id"));
        
        all_vis[exist_plot.attr("id")].graph.gene.push(clicked_gene);
    }
    else
    {
        ///if not then do the query to create it
        //first need to do an ajax call to get the appropriate data
        
        $.post("/HitWalker2/get_data/", post_data, function(data, status, xhr)
                           {
                                if (status == "success")
                                {
                                    parsed_data = JSON.parse(data);
                                    parsed_data.graph['gene'] = [clicked_gene];
                                    add_to_image(parsed_data, 'waterfall_plot', {type:post_data.type ,sample:post_data.sample, gene:clicked_gene});
                                
                                }
                                else
                                {
                                    console.log(status)
                                }
                            }, "text");
    }
    
    $('[data-toggle="popover"]').popover('destroy');
        $(popover_ref).popover('destroy');
        popover_ref = null;
        d3.selectAll('.popover').remove();
    
}

function draw_waterfall_plot(vis, dta,w,h, shiftKey)
{
    //from http://bl.ocks.org/mbostock/3885705
    
    var margin = {left:30, top:5};
    w = w - 30;
    h = h-60;
    
    var x = d3.scale.ordinal()
        .rangeBands([0, w], .25, 1);
    
    var y = d3.scale.linear()
        .range([h, 0]);
        
    var xAxis = d3.svg.axis()
        .scale(x)
        .orient("bottom");
    
    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left");
    
    var sub_vis = vis.append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");
    
    var sel_x_name = d3.set();
    var x_dict = d3.map();
    
    dta.data.forEach(function(d,i,ar)
                {
                   var gene_set = d3.set(d.db_name);
                   dta.gene.forEach(function(d2, i2, ar2)
                   {
                        if (gene_set.has(d2))
                        {
                            sel_x_name.add(d.x_name);
                        }
                   })
                    
                    x_dict.set(d.x_name, d.db_name);
                })
    
    x.domain(dta.data.map(function(d) { return d.x_name; }));
    y.domain([0, d3.max(dta.data, function(d) { return d.value; })+20]);
    
    //might also make text invisble or shrink it relative to the number of genes being plotted
    
    sub_vis.append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + h + ")")
        .call(xAxis)
        .selectAll("text")
            .attr("y", -2)
            .attr("x", -10)
            //.attr("dy", ".100em")
            .attr("transform", "rotate(-90)")
            .style({'text-anchor':'end','font-family':'sans-serif', 'font-size':'8px'})
            //.style({'text-anchor':'end', 'font':'8px sans-serif'})
            .attr("class", function(d,i)
                  {
                        if (sel_x_name.has(d))
                        {
                            return("plot text highlight")
                        }
                        else
                        {
                            return("plot text");
                        }
                  })
            //used for updating graph if user selects sample and graph type already displayed.
            .attr("data-genel", function(d,i)
                  {
                        return(x_dict.get(d));
                  });
    sub_vis.append("g")
        .attr("class", "y axis")
        .call(yAxis)
        .selectAll("text").style({'font-family':'sans-serif', 'font-size':'10px'});
        //for adding a Y axis eventually...
      //.append("text")
      //  .attr("transform", "rotate(-90)")
      //  .attr("y", 6)
      //  .attr("dy", ".71em")
      //  .style("text-anchor", "end")
      //  .text("Value");
    
    sub_vis.selectAll(".bar")
        .data(dta.data)
        .enter().append("rect")
        .attr("class", function(d,i)
              {
                var gene_set = d3.set(d.db_name);
                var should_high = false;
                
                dta.gene.forEach(function(d2, i2, ar2)
                                 {
                                    if(gene_set.has(d2))
                                    {
                                        should_high = true;
                                    }
                                 });
                if (should_high)
                {
                    return("bar highlight")
                }
                else
                {
                    return("bar")
                }
              })
        .attr("x", function(d) { return x(d.x_name); })
        .attr("width", x.rangeBand())
        .attr("y", function(d) { return y(d.value); })
        .attr("height", function(d) { return h - y(d.value); });
        
        // //again used for updating graph if user selects sample and graph type already displayed.
        //.attr("data-genel", function(d,i)
        //          {
        //                return(x_dict.get(d));
        //          });;
    
    }