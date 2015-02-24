
 //probably need to change the organization of graph_stack such that it points to the current graph and not the previous...
//var move_in_stack = function(direction)
//{
//    if ((direction == -1 && stack_pos > 0) || (direction == 1 && stack_pos < (graph_stack.length - 1)))
//    {
//        //console.log(graph_stack[0])
//        stack_pos += direction;
//        var temp_graph = JSON.parse(graph_stack[stack_pos]);
//        var new_links = [];
//        
//        var node_map = d3.map();
//        
//        for (i=0; i < temp_graph.nodes.length;i++)
//        {
//            node_map.set(temp_graph.nodes[i].id, i);
//        }
//        
//        for (i=0; i < temp_graph.links.length;i++)
//        {
//            new_links.push({'source':node_map.get(temp_graph.links[i].source.id), 'target':node_map.get(temp_graph.links[i].target.id), 'attributes':temp_graph.links[i].attributes});
//        }
//        
//        graph_obj.nodes = temp_graph.nodes;
//        graph_obj.links = new_links;
//        
//        //remove_menu_items("node_att_div");
//        //remove_menu_items("filter_div");
//        //remove_menu_items("text_div");
//        
//        selected_node = d3.map();
//        update_graph();
//    }
//}

 var keyflip = function () {
    
    shiftKey = d3.event.shiftKey || d3.event.metaKey;
    //console.log(shiftKey)
    if (shiftKey == null)
    {
        
        return false;
    }
    else
    {
        return shiftKey
    }
}



//var updateLink = function()
//{
//    this.attr("x1", function(d) {
//                            return d.source.x;
//                    }).attr("y1", function(d) {
//                            return d.source.y;
//                    }).attr("x2", function(d) {
//                            return d.target.x;
//                    }).attr("y2", function(d) {
//                            return d.target.y;
//                    });
//}

var updateLink = function()
{
    
    this.attr("d", function(d,i)
              {
                    if (d.edge_total == 1 || d.cur_edge == 1)
                    {
                        return ("M"+d.source.x+","+d.source.y+"L"+d.target.x+","+d.target.y);
                    }
                    else
                    {
                        
                        if (d.cur_edge % 2 == 0)
                        {
                            cor_dir = 1 * 5 * Math.ceil(d.cur_edge/3)
                        }
                        else
                        {
                            cor_dir = -1 * 5 * Math.ceil(d.cur_edge/3)
                        }
                        
                        var dx = d.target.x - d.source.x,
                            dy = d.target.y - d.source.y
                        
                        var source_angle = Math.atan(dy/dx);
                        
                        if (Math.abs(source_angle) > (3.14/4))
                        {
                            return ("M"+(d.source.x+cor_dir)+","+d.source.y+"L"+(d.target.x + cor_dir)+","+d.target.y);
                        }
                        else
                        {
                            return ("M"+d.source.x+","+(d.source.y+cor_dir)+"L"+d.target.x+","+(d.target.y+cor_dir));
                        }
                        
                        
                        //var dx = d.target.x - d.source.x,
                        //    dy = d.target.y - d.source.y
                        //    dr = Math.sqrt(dx*dx+dy*dy);
                        //
                        //return("M"+d.source.x+","+d.source.y+"A"+dr + "," + dr + " 0 0,1" + d.target.x + "," + d.target.y)
                    }
              });
    
}

var updateTextLink = function() {
                    
                    this.attr("x1", function(d) {
                            return d.source.x;
                    }).attr("y1", function(d) {
                            return d.source.y;
                    }).attr("x2", function(d) {
                            return d.target.x;
                    }).attr("y2", function(d) {
                            return d.target.y;
                    });
            }


//var rectUpdateNode = function() {
//                    this.attr("transform", function(d) {
//                            return "translate(" + d.x + "," + d.y + ")";
//                    });
//
//            }
   
var updateNode = function() {
    
                    this.attr("transform", function(d) {
                            //right_dist = w-d.x
                            //console.log(w);
                            //if ((sel_val.x + sel_val.size) >  w) || ((sel_val.x - sel_val.size) < 0) || ((sel_val.y + sel_val.size) > h) || ((sel_val.y - sel_val.size) < 0))
                            //was just d.x and d.y, example from http://mbostock.github.io/d3/talk/20110921/bounding.html
                            return "translate(" +  d.x + "," + d.y + ")";
                    });

            }
            
      
function remove_menu_items(menu_id)
    {
        //console.log(menu_id)
        var old_div = document.getElementById(menu_id);
        while (old_div.firstChild) {
            old_div.removeChild(old_div.firstChild);
        }
        
    }

function make_text_div (text, is_bold, color, width)
{
    var temp_div = document.createElement('div');
        temp_div.innerHTML = text;
        temp_div.style.width = width;
        temp_div.style.color = color
        
        if (is_bold)
        {
            temp_div.style.fontWeight = "bold";
        }
        
    
    return(temp_div);
}

//function node_click_highlight (node_obj, vis, selected_node, shiftKey)
//{
//    var use_data = d3.select(node_obj).datum();
//    
//    //if (selected_node.keys().length == 0 && group_members.keys().length == 0)
//    //{
//    //    initial_node_type = use_data.attributes.node_type;
//    //}
//    //
//    if (selected_node.has(use_data.id) == false && (selected_node.keys().length == 0 || shiftKey == true))
//    {
//        selected_node.set(use_data.id, node_obj);
//    }
//    else if (selected_node.has(use_data.id) == false && selected_node.keys().length > 0)
//    {
//        selected_node = d3.map();
//        selected_node.set(use_data.id, node_obj);
//    }
//    
//    //trying to think about group-node moving would work...
//    //vis.selectAll("g.node").selectAll("g.nodeDrag").enter().insert('svg:g').attr("class", "nodeDrag");
//    
//    //vis.selectAll("g.node").filter(function(d,i)
//    //                               {
//    //                                    console.log(d.id)
//    //                                    return (selected_node.has(d.id)==false)
//    //                               }).insert('svg:g', ':first-child').attr("class", "nodeDrag");
//    
//    //.insert('svg:g').attr("class", "nodeDrag")
//    
//    //vis.selectAll("g.node")
//    //vis.selectAll("g.labelNode").data(force2.nodes(), function(d) {return d.id;}).enter().append('svg:g').attr("class", "labelNode");
//    //example of exiting...
//    // var anchorLink = vis.selectAll("line.labelEdges").data(labelEdges);
//    //anchorLink.exit().remove();
//    //console.log(force2.nodes().length);
//    
//    //vis.selectAll("line.link").classed("Unselected", true);
//    //vis.selectAll("g.node").each(function(d,i)
//    //{
//    //    var main_circle = d3.select(this).selectAll("circle");
//    //    
//    //    if (selected_node.has(d.id))
//    //    {
//    //        main_circle.classed("Unselected", false);
//    //        //return false;//1;
//    //    }
//    //    else
//    //    {
//    //        main_circle.classed("Unselected", true);
//    //        //return true;//back_opacity;
//    //    }
//    //});
//    //
//    //vis.selectAll("g.labelNode").selectAll("text").classed("Unselected", function(d,i)
//    //{
//    //    if (selected_node.has(d.node.id))
//    //    {
//    //        return false;//1;
//    //    }
//    //    else
//    //    {
//    //        return true;//back_opacity;
//    //    }
//    //});
//    
//}
//


 function get_node_size (cur_nodes)
{
    
    
    for (i=0; i < cur_nodes.length;i++)
    {
        if (cur_nodes[i].attributes.node_type == "MetaNode")
        {
            cur_nodes[i]['size'] = Math.min(50, cur_nodes[i].children.length*10);
        }
        else
        {
            cur_nodes[i]['size'] = 20;
        }
        
    }
   
    
    //var node_scores = [];
    //
    ////figure out the number of categories
    //
    //for (i=0; i < cur_nodes.length;i++)
    //{
    //    if (d3.map(cur_nodes[i]).has("node_score"))
    //    {
    //        node_scores.push(cur_nodes[i].node_score);
    //    }
    //}
    //
    //node_scores = node_scores.sort(function(a,b){return a-b});
    //
    //var cur_range = d3.scale.ordinal();
    //
    //cur_range.rangePoints([10, 10+(1.5*node_scores.length)]);
    //cur_range.domain(node_scores);
    //
    //for (i=0; i < cur_nodes.length;i++)
    //{
    //    if (d3.map(cur_nodes[i]).has("node_score"))
    //    {
    //        cur_nodes[i]['size'] = cur_range(cur_nodes[i].node_score);
    //    }
    //}
    //
    return(cur_nodes);
}

function get_line_size (cur_edge)
{
    var max_size = 5;
    
    var cur_range = d3.scale.ordinal()
    
    cur_range.rangePoints([1,max_size], 0)
    
    var range_list = []
    
    for (i=1; i <= max_size; i++)
    {
        range_list.push(i);
    }
    
    cur_range.domain(range_list)
    
    var use_size = 1;
    
    if ("score" in cur_edge)
    {
        use_size = Math.round(cur_edge.score/200);
    }
    else if ("rank" in cur_edge)
    {
        use_size = (max_size - Math.min(max_size, cur_edge.rank))+1;
    }
    else
    {
        console.log("WARNING: Unexpected score/rank value in cur_edge")
    }
    
    cur_edge['size'] = cur_range(use_size);
   // console.log(JSON.stringify(cur_edge));
    return(cur_edge);
}

function rev_type_translation(value)
{
    var new_obj = {}
    
    var ret_val = value;
    
    for (i in node_transl)
    {
        if (value == node_transl[i])
        {
            ret_val = i;    
        }
    }
    
    return(ret_val);
    
}

function type_translation(value, transl_type, result_type)
{
         if(typeof(result_type) === 'undefined') result_type = 'class'
      
         if(typeof(transl_type)==='undefined') transl_type = node_transl;
         
         if (result_type == 'class')
         {
            if (value in transl_type == false || transl_type[value].class == null)
            {
               //set to one of the 'default' colors
               
               var remaining_defaults = d3.set();
               
               default_css.forEach(function(x)
                                   {
                                       if (used_defaults.has(x) == false)
                                       {
                                          remaining_defaults.add(x);
                                       }
                                   })
               
              if (remaining_defaults.size == 0)
              {
                  //start over
                  var use_val = default_css.values()[0];
                  used_defaults = d3.set([use_val]);
                  
                  
              }else{
                //choose the first available one
                  var use_val = remaining_defaults.values()[0];
                  used_defaults.add(use_val);
               
              }
              
              if (value in transl_type == false)
              {
                  transl_type[value] = {name:value, class:use_val}
              }else{
                  transl_type[value].class = use_val;
              }
              
              return(use_val);
              
            }else{
               return (transl_type[value].class)
            }
            
         }else{
            if (value in transl_type)
            {
                return (transl_type[value].name);
            }
            else
            {
                return (value);
            }
         } 
}
    
function make_node_map (keep_nodes)
{
    var keep_node_map = d3.map();
    
    for (i=0; i < keep_nodes.length;i++)
    {
        if (keep_node_map.has(keep_nodes[i].id))
        {
            console.log("ERROR: Found duplicate node name in make_node_map");
        }
        else
        {
            keep_node_map.set(keep_nodes[i].id, keep_nodes[i]);
        }
    }
    
    return(keep_node_map);
}

//function make_link_map (keep_links)
//{
//    var keep_link_map = d3.map();
//    
//    for (i = 0; i < keep_links.length;i++)
//    {
//        if (keep_link_map.has(keep_links[i].source.id) == false)
//        {
//            var temp_source_map = d3.map();
//            temp_source_map.set(keep_links[i].target.id, keep_links[i]);
//            keep_link_map.set(keep_links[i].source.id, temp_source_map);
//        }
//        else
//        {
//            if (keep_link_map.get(keep_links[i].source.id).has(keep_links[i].target.id) == false)
//            {
//                keep_link_map.get(keep_links[i].source.id).set(keep_links[i].target.id, keep_links[i]);
//            }
//        }
//        
//        if (keep_link_map.has(keep_links[i].target.id) == false)
//        {
//            var temp_source_map = d3.map();
//            temp_source_map.set(keep_links[i].source.id, keep_links[i]);
//            keep_link_map.set(keep_links[i].target.id, temp_source_map);
//        }
//        else
//        {
//            if (keep_link_map.get(keep_links[i].target.id).has(keep_links[i].source.id) == false)
//            {
//                keep_link_map.get(keep_links[i].target.id).set(keep_links[i].source.id, keep_links[i]);
//            }
//        }
//    }
//    
//    return(keep_link_map);
//}

function add_existing_links_to_graph(graph_obj, new_graph)
{
    var node_pos = d3.map();
    var new_nodes = [];
    var keep_node_map = make_node_map(graph_obj.nodes);
    
    new_graph.nodes.forEach(function(d,i,ar)
                            {
                                node_pos.set(d.id, i);
                                
                                if (keep_node_map.has(d.id))
                                {
                                    new_nodes.push(keep_node_map.get(d.id));
                                }
                                else
                                {
                                    new_nodes.push(d);
                                }
                            })
    
    graph_obj.links.forEach(function(d,i,ar)
                            {
                                var temp_link = d;
                                
                                temp_link.source = node_pos.get(d.source.id);
                                temp_link.target = node_pos.get(d.target.id);
                                new_graph.links.push(temp_link);
                            })
    
    
    new_graph.nodes = new_nodes;
    return(new_graph)
}

function process_nodes_to_copy(node_list)
{
    var ret_list = [];
    
    node_list.forEach(function(d,i,ar)
                      {
                            if (d.attributes.node_type != 'MetaNode')
                            {
                                ret_list.push({id:d.id, node_type:d.attributes.node_type});
                            }
                            else
                            {
                                //in this case iterate through the MetaNode children and directly append them to the list
                                d.children.forEach(function(d2,i2,ar2)
                                                   {
                                                        ret_list.push({id:d2.id, node_type:d2.attributes.node_type});
                                                   });
                            }
                      })
    
    return(ret_list);
}

function make_unique_meta_id(vis_list)
{
    var meta_node_total = 1;
    
    for(panel in vis_list)
    {
        if ('nodes' in vis_list[panel].graph && 'links' in vis_list[panel].graph)
        {
            vis_list[panel].graph.nodes.forEach(function(d,i)
                                       {
                                            if (d.attributes.node_type == "MetaNode")
                                            {
                                                vis_list[panel].graph.nodes[i].id = 'meta_node_' + meta_node_total;
                                                meta_node_total += 1;
                                            }
                                       });
        }
        
    }
    
    return(vis_list);
}

function copy_nodes(query_nodes, subj_nodes, g_ind)
{
    var target = document.getElementById('network_spinner');
    var spinner = new Spinner(spinner_opts).spin(target);
    spinner.spin(target);
    
    var use_queries = process_nodes_to_copy(query_nodes);
    var query_set = d3.set(use_queries.map(function(d) {return d.id;}));
    
    var cur_subj =  process_nodes_to_copy(subj_nodes);
    
    var query_names = query_nodes.map(function(d)
                                      {
                                          if (d.attributes.node_type == "MetaNode")
                                          {
                                             if (d.children.length > 3)
                                             {
                                                return(d.children[0].display_name + '...' + d.children[d.children.length-1].display_name);
                                             }else{
                                                return(d.children.map(function(e) {e.display_name}).join(','));
                                             }
                                          }else{
                                             return(d.display_name);
                                          }
                                          
                                      });
    
    if ('previous_queries' in all_vis[g_ind])
    {
        //if not already there, add previous queries to query nodes
        all_vis[g_ind].previous_queries.forEach(function(d,i,ar)
                                                {
                                                    if (query_set.has(d.id) == false)
                                                    {
                                                        use_queries.push(d);
                                                    }
                                                    
                                                });
        
    }
    
    //remake query_set
    query_set = d3.set(use_queries.map(function(d) {return d.id}));
    
    //globally subtract the query nodes from the subject nodes
    var use_subj = [];
    
    cur_subj.forEach(function(d,i,ar)
                     {
                        if (query_set.has(d.id)==false)
                        {
                            use_subj.push(d);
                        }
                     });
    
    var post_node = {query:JSON.stringify(use_queries), subj:JSON.stringify(use_subj)};
    
    //the nodes will be "copied" to the new panel and new edges will be added to them
   $.post("/HitWalker2/copy_nodes/", post_node, function(data, status, xhr)
      {
           if (status == "success")
           {
               var new_graph = JSON.parse(data)
               
               console.log(history_dict[g_ind])
               
               if (history_dict[g_ind].action == "blank_panel" && Object.keys(history_dict[g_ind].values).length == 0)
               {
                  var new_panel_id = g_ind;
                  
                  history_dict[new_panel_id] = {action:'blank_panel', status:'active', values:{prev:undefined, addition:query_names}, title:query_names.join(",")};
                  
               }else{
                  var new_panel_id = make_unique_id();
                  history_dict[new_panel_id] = {action:'addition', status:'active', values:{prev:g_ind, addition:query_names}};
                  
               }
               
                all_vis[new_panel_id] = {graph:new_graph, previous_queries:use_queries};
               update_image(all_vis);
               
               //console.log(JSON.stringify(query_names));
                
               
               
               spinner.stop();
               
               d3.selectAll("g.node").each(function(d,i)
               {
                 d.fixed=false;
               });
                
                //if (force2 != undefined)
                //{
                //    //s d.fixed = true; // of course set the node to fixed so the force doesn't include the node in its auto positioning stuff
                //     node_tick(force2, node, back_node, anchorNode, anchorLink, link);
                //}
              
                //force.resume();
           }
           else
           {
               console.log(status)
           }
       
           
      }, "text");
}

//depricated in favor of direct code in update_image
//function()
//
//function add_to_graph(graph_obj, new_graph)
//        {
//            //console.log(JSON.stringify(new_graph));
//            if ('nodes' in new_graph && 'links' in new_graph)
//            {
//                console.log("new", JSON.stringify(d3.nest().key(function(d) d.attributes.type).rollup(function(d) d.length).map(new_graph.links)));
//                console.log("old", JSON.stringify(d3.nest().key(function(d) d.attributes.type).rollup(function(d) d.length).map(graph_obj.links)));
//                
//                var keep_node_map = make_node_map(graph_obj.nodes);
//                var keep_nodes = [];
//                var keep_links = [];
//                
//                //add the new nodes to the keep_node_map
//                
//                for (i = 0; i < new_graph.nodes.length;i++)
//                {
//                    if (keep_node_map.has(new_graph.nodes[i].id) == false)
//                    {
//                        keep_node_map.set(new_graph.nodes[i].id, new_graph.nodes[i])
//                    }
//                }
//                
//                for (i=0; i < graph_obj.links.length;i++)
//                {
//                    if (keep_node_map.has(graph_obj.links[i].source.id) && keep_node_map.has(graph_obj.links[i].target.id))
//                    {
//                        keep_links.push(graph_obj.links[i]);
//                    }
//                }
//                
//                console.log(JSON.stringify(keep_node_map.keys()));
//                
//               for (i = 0; i < new_graph.links.length;i++)
//               {
//                    console.log('new_link',JSON.stringify(new_graph.links[i]))
//                    if (keep_node_map.has(new_graph.links[i].source.id) && keep_node_map.has(new_graph.links[i].target.id))
//                    {
//                        console.log('in')
//                        keep_links.push(new_graph.links[i]);
//                    }
//               }
//                
//                keep_nodes = keep_node_map.values()
//                
//                console.log("combined", JSON.stringify(d3.nest().key(function(d) d.attributes.type).rollup(function(d) d.length).map(keep_links)))
//                
//                return ({graph:{nodes:keep_nodes, links:keep_links}});
//            }
//            else
//            {
//                return ({graph:graph_obj});   
//            }
//            
//        }

function make_graph_legend(vis, node_list, edge_list)
{
    vis.selectAll("g.legend").data([]).exit().remove();
    vis.selectAll("g.link_legend").data([]).exit().remove();
    
    var legend_hier = d3.layout.tree();
    
    var temp_nodes = {attributes:{node_type:'MetaNode'}, children:JSON.parse(JSON.stringify(node_list))};
    
    var legend_nodes = legend_hier.nodes(temp_nodes);
    
    var legend_links = legend_hier.links(legend_nodes);
    //jQuery.extend(true, {}, 
    
    var temp_nest = d3.nest().key(function(d) {return d.source.attributes.node_type}).key(function(d) {return d.target.attributes.node_type}).rollup(function(d) {return d.length}).map(legend_links, d3.map);
    
    var rev_nest = d3.nest().key(function(d) {return d.target.attributes.node_type}).key(function(d) {return d.source.attributes.node_type}).rollup(function(d) {return d.length}).map(legend_links, d3.map);
    
    //need to check whether any (which) of the root keys are in the values.  If so organize the hierarchies such that  
    var root_keys = temp_nest.get("MetaNode").keys();
    
    var legend_map = [];
    
    function recurse_keys (keys, cur_list)
    {
         keys.forEach(function(d)
                        {
                            var temp_obj = {name:d, type:"Node", children:[]};
                            
                           if (temp_nest.has(d))
                           {
                                var temp_keys = temp_nest.get(d).keys();
                                
                                recurse_keys(temp_keys, temp_obj.children);
                           }
                           
                           cur_list.push(temp_obj);
                           
                        });
    }
    
    root_keys.forEach(function(d)
                      {
                            //if d has more than just root as its child, it should be excluded as a duplicate
                            //also remove any MetaNodes from consideration
                            if (rev_nest.get(d).keys().length == 1 && d != "MetaNode")
                            {
                                var cur_obj = {name:d, type:"Node", children:[]};
                        
                                if (temp_nest.has(d))
                                {
                                    var cur_keys =  temp_nest.get(d).keys();
                                    
                                    recurse_keys(cur_keys, cur_obj.children);
                                }
                                
                                legend_map.push(cur_obj);
                            }
                            
                      });
    
    //roughly compute the amount of space necessary--also propagate it to the link legends as well
    var child_num = 0;
    
    function child_counter (list)
    {
        var counter = 0;
       list.forEach(function(d,i, ar)
       {
            if (d.children.length > 0)
            {
                counter += d.children.length;
                counter += child_counter(d.children);
            }
          
       })
        
        return(counter)
    }
    
    //add in the links too, at the level of the nodes
    var link_map = d3.set();
    //edge_list.forEach(function(d,i,ar)
    //{
    //    var transl_d = type_translation(d, edge_transl)
    //    if (link_map.has(transl_d) == false)
    //    {
    //        link_map.add(transl_d);
    //        legend_map.push({name:transl_d, type:"Link", children:[]});
    //    }
    //})
    
    edge_list.forEach(function(d,i,ar)
      {
         if (link_map.has(d) == false){
            link_map.add(d);
            legend_map.push({name:d, type:"Link", children:[]});
         }
      })
    
    child_num = child_counter(legend_map);
    
    child_num += link_map.values().length;
    
    var legend_tree = d3.layout.tree().size([child_num*15,100]);
    
    var legend_nodes = legend_tree.nodes({name:'root', children:legend_map})
                        .filter(function(d)
                        {
                           return(d.name != 'root'); 
                        });
    
    var legend_links = legend_tree.links(legend_nodes);
    //remove the root nodes
    
    var diagonal = d3.svg.diagonal()
        .projection(function(d) { return [d.y, d.x]; });
    
    var legend= vis.append("g").attr("class", "legend");//.attr("transform", "translate(5,0)");
    
    legend.selectAll("path.legend_link")
                .data(legend_links)
                .enter().append("path")
                .attr("class", "legend_link")
                .attr("d", diagonal);
    
    //need a manual offset here as for some reason it draws them as overlapping...
    var legend_labs = legend.selectAll("g")
        .data(legend_nodes)
        .enter().append("g")
        .attr("transform", function(d){
            
            if (d.type == "Node")
            {
                return "translate(" + d.y + "," + d.x + ")";
            }
            else
            {
                return "translate(" + d.y + "," + (d.x+30) + ")";
            }
            
            });
    
    legend_labs.append("text")
            .attr("dx", function(d) { return d.children.length > 0 ? -8 : 8; })
            .attr("dy", function(d) { return (5) })
            .attr("text-anchor", function(d) { return d.children.length > 0 ? "end" : "start"; })
            .style("font-size", "10px")
            .text(function(d) {
               if (d.type == "Node")
               {
                  return (type_translation(d.name, node_transl, 'name').replace("_", " "));
               }else{
                  return (type_translation(d.name, edge_transl, 'name').replace("_", " "));
                }});
    
    
    node_labs = legend_labs.filter(function(d,i) {return(d.type == "Node")});
    link_labs = legend_labs.filter(function(d,i) {return(d.type == "Link")});
    
    node_labs.append("circle").attr("r", 5).attr("class", function(d)
                                                    {
                                                       return(type_translation(d.name)); 
                                                    });
    link_labs.append("line").attr("class", function(d)
                            {
                               return(type_translation(d.name, edge_transl, 'class')); 
                            })
                            .attr("x1", -10)
                            .attr("x2", 5)
                            .style("stroke-width", 3);
}
function deemphasize_nodes(vis, selected_node, is_hover)
{
    if (selected_node.keys().length == 0)
    {   
        vis.selectAll("g.node").selectAll("circle[class]").classed("Unselected", false);
        vis.selectAll("path.link").classed("Unselected", false);
        vis.selectAll("g.labelNode").selectAll("text").classed("Unselected", false);
    }else{
        
        if (vis.attr('id') == selected_panel || is_hover)
        {
            vis.selectAll("g.node").each(function(d,i)
            {
                var main_circle = d3.select(this).selectAll("circle[class]");
                
                if (selected_node.has(d.id))
                {
                    main_circle.classed("Unselected", false);
                }
                else
                {
                    main_circle.classed("Unselected", true);
                }
            });
            //vis.selectAll("line.link").style("opacity", back_opacity);
            vis.selectAll("path.link").classed("Unselected", true);
            //vis.selectAll("g.labelNode").selectAll("text").style("opacity", function(d,i)
            vis.selectAll("g.labelNode").selectAll("text").classed("Unselected", function(d,i)
                                                              {
                                                                if (selected_node.has(d.node.id))
                                                                {
                                                                    return false;//1;
                                                                }
                                                                else
                                                                {
                                                                    return true;//back_opacity;
                                                                }
                                                                });
                }
    }
}

//demphasize all nodes except those that are chosen in the current panel
function unselect_all(panel_id)
{
   
    for(i in all_vis)
    {
        if (i != panel_id && i != selected_panel)
        {
            var temp_panel = d3.select("#"+i);
            temp_panel.selectAll("g.node").selectAll("circle[class]").classed("Unselected", true);
            temp_panel.selectAll("path.link").classed("Unselected", true);
            temp_panel.selectAll("g.labelNode").selectAll("text").classed("Unselected", true);
        }
    }
    
   
}

function emphasize_nodes(cur_obj, cur_name, selected_node, vis, is_click)
{
    var cur_panel = vis.attr('id');
    
    if (is_click)
    {
        var use_data = d3.select(cur_obj).datum();
        
        
        if (selected_node.has(use_data.id) == false && (selected_node.keys().length == 0 || shiftKey == true) && (selected_panel == cur_panel || selected_panel == null))
        {
            selected_node.set(use_data.id, cur_obj);
            selected_panel = cur_panel;
        }
        else if ((selected_node.has(use_data.id) == false && selected_node.keys().length > 0) || (selected_panel != cur_panel))
        {
            selected_node = d3.map();
            selected_node.set(use_data.id, cur_obj);
            selected_panel = cur_panel;
        }
        
    }
    
    var cur_ids = d3.map();
    
    if (selected_node.keys().length == 0)
    {
        cur_ids.set(cur_name, true);
    }
    else
    {
        cur_ids = selected_node;
    }
    
    unselect_all(cur_panel);
    
    deemphasize_nodes(vis, cur_ids, is_click==false);
    
    
    ////vis.selectAll("g.node").selectAll("circle").style("opacity", function(d,i)
    // vis.selectAll("g.node").each(function(d,i)
    //{
    //    var main_circle = d3.select(this).selectAll("circle[class]");
    //    
    //    if (cur_ids.has(d.id))
    //    {
    //        main_circle.classed("Unselected", false);
    //    }
    //    else
    //    {
    //        main_circle.classed("Unselected", true);
    //    }
    //});
    //////vis.selectAll("line.link").style("opacity", back_opacity);
    //vis.selectAll("path.link").classed("Unselected", true);
    //////vis.selectAll("g.labelNode").selectAll("text").style("opacity", function(d,i)
    //vis.selectAll("g.labelNode").selectAll("text").classed("Unselected", function(d,i)
    //                                                  {
    //                                                    if (cur_ids.has(d.node.id))
    //                                                    {
    //                                                        return false;//1;
    //                                                    }
    //                                                    else
    //                                                    {
    //                                                        return true;//back_opacity;
    //                                                    }
    //                                                    });
    
    return (selected_node);
}


function add_to_image(parsed_data, action, values)
{
    var new_panel_id = make_unique_id();
    
    history_dict[new_panel_id] = {action:action, title:parsed_data.title, status:'active', values:values};
    
    all_vis[new_panel_id] = {type:action, graph:parsed_data.graph};
    
    update_image(all_vis);
    
   d3.selectAll("g.node").each(function(d,i)
                                {
                                  d.fixed=false;
                                });
}

function post_to_fullfill (obj)
{
    var target = document.getElementById('network_spinner');
    var spinner = new Spinner(spinner_opts).spin(target);
    spinner.spin(target);
    
    var post_node = {choice:d3.select(obj).attr("data-value"), nodes:''};
    
    var node_dict = {};
    
    selected_node.forEach(function(key, value)
                          {
                                var cur_val = d3.select(value).datum();
                                var cur_val_list = [];
                                
                                //deal with the case that cur_val refers to a metanode
                                if (cur_val.attributes.node_type == 'MetaNode')
                                {
                                    cur_val.children.forEach(function(d){cur_val_list.push(d);});
                                }
                                else{
                                    cur_val_list.push(cur_val);
                                }
                                
                                cur_val_list.forEach(function(d)
                                                     {
                                                        if (d.attributes.node_type in node_dict)
                                                        {
                                                            node_dict[d.attributes.node_type].push({id:d.id, display_name:d.display_name});
                                                        }
                                                        else
                                                        {
                                                            node_dict[d.attributes.node_type] = [{id:d.id, display_name:d.display_name}]   
                                                        }
                                                     });
                          });
    
    post_node.nodes = JSON.stringify(node_dict);
    
    $.post("/HitWalker2/fullfill_node_query/", post_node, function(data, status, xhr)
       {
            if (status == "success")
            {
                //now add the nodes, etc to the overall graph
                
                //if no data was found, add another popover pointing to obj stating such
                
                parsed_data = JSON.parse(data);
                
                if (parsed_data.graph.nodes.length == 0)
                {
                  
                  $(obj).siblings('a.list-group-item').each(function()
                                                            {
                                                               $(this).popover('destroy');
                                                            });
                  
                  var pop_content = '<p class="text-danger">Sorry, no results were found for this query.</p>';
                    
                  $(obj).popover({container:"body", html:true, placement: 'right', title:'', content:pop_content, trigger:"manual"});
                    
                  $(obj).popover('show');
                }
                else{
                    add_to_image(parsed_data, 'query', selected_node.keys());
                
                    if (popover_ref != null)
                    {
                        delete_popover();
                    }
                }
                
                spinner.stop()
               
            }
            else
            {
                console.log(status)
            }
        
       }, "text");
}

//var set_groups = function()
//{
//    //console.log(JSON.stringify(selected_node));
//    
//    //add a new node to graph_obj.nodes
//    //remove the nodes in selected_node and place them as children in the new node
//    //remove the connections to the selected_nodes
//    //add the connections to the new group node
//    
//    //TODO:  Need to treat length == 1 as a special case, seems to work for regular gene with no innards, not so for the other cases
//    
//    //seperate functionality for exploding a group
//    
//    if (selected_node.keys().length > 0)
//        {
//            
//            var keep_links = [];
//            var keep_nodes = [];
//            var temp_nodes = d3.set();
//            
//            //need to make size based on the sizes of the nodes in selected_node
//            
//            var new_node = {size:0, id:"User_Group_"+cur_group, display_name:"User_Group_"+cur_group, children:[], attributes:{node_type:"User_Group", add_string:""}};
//            
//            for (i=0; i < graph_obj.nodes.length;i++)
//            {
//                if (selected_node.has(graph_obj.nodes[i].id) == true)
//                {
//                    if (graph_obj.nodes[i].attributes.node_type != "User_Group")
//                    {
//                        new_node.children.push(graph_obj.nodes[i]);
//                        new_node.size += graph_obj.nodes[i].size;
//                    }
//                    else
//                    {
//                        for(j=0; j < graph_obj.nodes[i].children.length; j++)
//                        {
//                            new_node.children.push(graph_obj.nodes[i].children[j]);
//                            new_node.size += graph_obj.nodes[i].children[j].size;
//                        }
//                        
//                    }
//                }
//                else
//                {
//                    keep_nodes.push(graph_obj.nodes[i])
//                }
//            }
//            
//            new_node.size = Math.min(30, new_node.size);
//            
//            keep_nodes.push(new_node);
//            
//            for (i=0; i < graph_obj.links.length;i++)
//            {
//                var temp_link;
//                
//                if (selected_node.has(graph_obj.links[i].source.id) == true && selected_node.has(graph_obj.links[i].target.id) == false)
//                {
//                    temp_link = graph_obj.links[i];
//                    temp_link.source = new_node;
//                    keep_links.push(temp_link);
//                    
//                }
//                else if (selected_node.has(graph_obj.links[i].target.id) == true && selected_node.has(graph_obj.links[i].source.id) == false)
//                {
//                    temp_link = graph_obj.links[i];
//                    temp_link.target = new_node;
//                    keep_links.push(temp_link);
//                }
//                else if (selected_node.has(graph_obj.links[i].target.id) == false && selected_node.has(graph_obj.links[i].source.id) == false)
//                {
//                    keep_links.push(graph_obj.links[i])
//                }
//            }
//            
//            
//            cur_group+=1;
//            graph_obj.nodes = keep_nodes;
//            graph_obj.links = keep_links;
//
//            selected_node = d3.map();
//
//            update_graph();
//
//            graph_stack.push(JSON.stringify(graph_obj));
//            //stack_pos = graph_stack.length - 1;
//        }
//        }

function node_tick(force2, node, back_node, anchorNode, anchorLink, link, is_dragged)
{
     force2.start();
    
    if (is_dragged == false)
    { 
        //was just d.x and d.y, example from http://mbostock.github.io/d3/talk/20110921/bounding.html for setting up the constraints
        node.each(function(d,i)
                  {
                    d.x=Math.max(d.size, Math.min(w - d.size, d.x));
                    d.y=Math.max(d.size, Math.min(h - d.size, d.y));
                  });
        
        back_node.each(function(d,i)
                       {
                            d.x=Math.max(d.size, Math.min(w - d.size, d.x));
                            d.y=Math.max(d.size, Math.min(h - d.size, d.y));
                       });
    }
    
    
    node.call(updateNode);
    back_node.call(updateNode);
    
    //subAnchorNode2.each(function(d, i) {
    anchorNode.each(function(d, i) {
            
            
            
            if(i % 2 == 0) {
                    
                    if (is_dragged)
                    {
                        d.x = d.node.x;
                        d.y = d.node.y;
                    }
                    else
                    {
                        d.x=Math.max(d.node.size, Math.min(w - d.node.size, d.node.x));
                        d.y=Math.max(d.node.size, Math.min(h - d.node.size, d.node.y));
                    }
                    
                    
            } else {
                     var text_width = 50;
                     var text_height = 5;
                     
                     d.x=Math.max(text_width, Math.min(w - text_width, d.x));
                     d.y=Math.max(text_height, Math.min(h - text_height, d.y));
                     
                    //console.log( this.childNodes[1])
                    //var b = this.childNodes[1].getBBox();
                    //hack as kept getting ns_error_failure upon adding additional graph to screen
                    var b = {width:50}
                    var diffX = d.x - d.node.x;
                    var diffY = d.y - d.node.y;
                    
                    var dist = Math.sqrt(diffX * diffX + diffY * diffY);
            
                    var shiftX = b.width * (diffX - dist) / (dist * 2);
                    shiftX = Math.max(-b.width, Math.min(0, shiftX));
                    var shiftY = 0;
                    
                    var rect_offset = parseFloat(this.childNodes[1].getAttribute("height"));
                    
                    this.childNodes[1].setAttribute("transform", "translate(" + shiftX + "," + (shiftY-rect_offset) + ")");
                    this.childNodes[2].setAttribute("transform", "translate(" + shiftX + "," + shiftY + ")");
                    //this.setAttribute("transform", "translate(" + shiftX + "," + shiftY + ")");
            }
            
    });
        
    //subAnchorNode2.call(updateNode);
    anchorNode.call(updateNode);
    //
    link.call(updateLink);
    anchorLink.call(updateTextLink);
    
                    
}

//
//maybe modify the implementation for this.
//http://bl.ocks.org/mbostock/4218871
//http://bl.ocks.org/mbostock/3231298
//use in_current_box to determine where the user wants to place the nodes
function in_current_box(cur_g, w, h)
{
    var cur_mouse = d3.mouse(cur_g.node());
    if ((cur_mouse[0] >= 0 && cur_mouse[0] <= w) && (cur_mouse[1] >= 0 && cur_mouse[1] <= h))
    {
        return (true);
    }
    else
    {
        return (false);
    }
}

function update_graph(vis, graph_obj,w,h, shiftKey)
        {
           
            var node;
            var back_node;
            var link;
            var anchorNode;
            var graph_stack = [];
            //var back_opacity = .3;
            //var max_query_size = 50;
            var cur_check_count = [];
           
            var stack_pos = -1;
            var initial_node_type = null;
            var group_members = d3.map();
            var cur_group = 0;
            var group_rect;
            //var group_colors = ['black', 'limegreen'];
            var node_drag_group;
            
            var x = d3.scale.linear()
        .domain([0, w])
        .range([0, w]);
        
        var y = d3.scale.linear()
            .domain([0, h])
            .range([h, 0]);
        
       var brush = vis.append("g")
        .datum(function() { return {selected: false, previouslySelected: false}; })
        .attr("class", "brush")
        .call(d3.svg.brush()
          .x(d3.scale.identity().domain([0, w]))
          .y(d3.scale.identity().domain([0, h]))
          .on("brushstart", function(d) {
            vis.selectAll("g.node").each(function(d) { d.previouslySelected = shiftKey && d.selected; });
          })
          
          .on("brush", function() {
            var extent = d3.event.target.extent();
            
            //ie. click on a panel but not on a node
            if (extent[0][0] == extent[1][0] && extent[0][1] == extent[1][1] && current_object===false)
            {
                //d3.selectAll('.popover').remove();
                if (popover_ref != null)
                {
                    
                    $('[data-toggle="popover"]').popover('destroy');
                    $(popover_ref).popover('destroy');
                    popover_ref = null;
                    d3.selectAll('.popover').remove();
                }
                
                selected_node = d3.map();
                selected_panel = null;
                deemphasize_nodes(vis, selected_node, false)
                unselect_all(vis.attr('id'));
                
                //vis.selectAll("g.node").selectAll("circle[class]").classed("Unselected", false)
                //vis.selectAll("path.link").classed("Unselected", false)
                //vis.selectAll("g.labelNode").selectAll("text").classed("Unselected", false)
                
                
                //
                //d3.select("#collapseFour").attr("class", "panel-collapse collapse").attr("style", "height:0px;")
                //d3.select("#collapseThree").attr("class", "panel-collapse collapse").attr("style", "height:0px;")
                
                force.stop();
                force2.stop();
                
            }
            else
            {
                vis.selectAll("path.link").classed("Unselected", true)
                vis.selectAll("g.node").each(function(d,i)
                {
                    d.selected = d.previouslySelected ^
                      (extent[0][0] <= d.x && d.x < extent[1][0]
                      && extent[0][1] <= d.y && d.y < extent[1][1]);
                    
                    if (d.selected)
                    {
                        if (selected_panel != vis.attr("id"))
                        {
                            selected_node = d3.map();
                            selected_panel = vis.attr("id");
                        }
                        
                        if (selected_node.has(d.id) == false)
                        {
                            selected_node.set(d.id, this)
                        }
                        
                        return d3.select(this).selectAll("circle[class]").classed("Unselected", false);//1
                       
                    }
                    else
                    {
                        selected_node.remove(d.id)
                        return d3.select(this).selectAll("circle[class]").classed("Unselected", true);//back_opacity;
                    }
                });
                
                vis.selectAll("g.labelNode").selectAll("text").classed("Unselected", function(d,i)
                {
                    var text_select = d.previouslySelected ^
                        (extent[0][0] <= d.x && d.x < extent[1][0]
                        && extent[0][1] <= d.y && d.y < extent[1][1]);
                        
                    var node_select = d.previouslySelected ^
                        (extent[0][0] <= d.node.x && d.node.x < extent[1][0]
                        && extent[0][1] <= d.node.y && d.node.y < extent[1][1]);
                    
                    if (node_select)
                    {   
                        return false;//1;
                    }
                    else
                    {
                        return true;//back_opacity;
                    }
                });
            }
          })
          .on("brushend", function() {
            
            //group_pop_func();
            
            d3.event.target.clear();
            d3.select(this).call(d3.event.target);
        }));
        
            //var use_charge = -3000 * 10*(1/graph_obj.nodes.length);
            
            var force = d3.layout.force()
             .gravity(1)
             .linkDistance(function(d)
                           {
                                 if (d.size)
                                 {
                                     return(d.size*30);
                                 }else{
                                     return(50);
                                 }
                           })
             //.linkStrength(10)
             .friction(.7)
             .charge(function(d)
                     {
                         if (d.attributes.node_type == "MetaNode")
                         {
                             return -3000 * 2;  
                         }else{
                              if (d.weight == 0)
                              {
                                 return(-100);
                              }else{
                                 return (-3000);
                              }
                         }
                     })
             .size([w, h]);
        
            var force2 = d3.layout.force()
                    .gravity(.1)
                    .linkDistance(0)
                    //.linkStrength(8)
                    .friction(.7)
                    .charge(-100)
                    .size([w, h]);
        
        //drag mechanism from    http://bl.ocks.org/norrs/2883411     
            var node_drag = d3.behavior.drag()
                .on("dragstart", dragstart)
                .on("drag", dragmove)
                .on("dragend", dragend);
            
            //limiting drags to left-click from http://stackoverflow.com/questions/17109474/how-to-prevent-d3-from-triggering-drag-on-right-click
            function dragstart(d, i) {
                if (d3.event.sourceEvent.which == 1)
                {
                     // stops the force auto positioning before you start dragging
                    force.stop()
                  
                    //var use_data = d3.select(p_node_obj).datum();
                    //
                    var use_data = d3.select(this).datum();
                    selected_node = emphasize_nodes(this, use_data.id, selected_node,  vis, true);
                    //if (selected_node.keys().length == 0)
                    //{
                    //    var use_data = d3.select(this).datum();
                    //    selected_node.set(use_data.id, this);
                    //}
                    
                    shouldDrag = true;
                }
              
            }
        
            function dragmove(d, i) {
            
                if (shouldDrag)
                {
                    selected_node.forEach(function(key, value)
                                      {
                                        var sel_val = d3.select(value);
                                        
                                        sel_val.attr("px", function(d) {d.px+= d3.event.dx});
                                        sel_val.attr("py", function(d) {d.py += d3.event.dy});
                                        sel_val.attr("x", function(d) {d.x += d3.event.dx});
                                        sel_val.attr("y", function(d) {d.y += d3.event.dy});
                                      });
                    
                   
                    node_tick(force2, node, back_node, anchorNode, anchorLink, link, true); // this is the key to make it work together with updating both px,py,x,y on d !
                }
            }
            
            
            
            function dragend(d, i) {
                
                if (shouldDrag)
                {
                
                    selected_node.forEach(function(key, value)
                                         {
                                           d3.select(value).datum().fixed=true;
                                         });
                   
                   //s d.fixed = true; // of course set the node to fixed so the force doesn't include the node in its auto positioning stuff
                    node_tick(force2, node, back_node, anchorNode, anchorLink, link,false);
                    force.resume();
                    //console.log(JSON.stringify());
                    shouldDrag = false;
                
                    var all_gs = d3.selectAll("g.g1");
                    //if the node(s) are dropped in another box than their own then copy them to the appropriate graph
                    if (all_gs[0].length > 1)
                    {
                        if (in_current_box(vis, w, h) == false)
                        {
                            var in_g = undefined;
                            var g_ind = "";
                            //figure out which box it was
                            all_gs.each(function(d,i)
                                        {
                                            if (in_current_box(d3.select(this), w, h))
                                            {
                                                in_g = this;
                                                g_ind = d3.select(this).attr("id");
                                            }
                                        });
                            
                            console.log(in_g)
                            
                            if (in_g != undefined)
                            {
                              
                                var subj_nodes = [];
                                var query_nodes = [];
                                
                                g_nodes = d3.select(in_g).selectAll("g.node");
                                g_nodes.each(function(d,i)
                                             {
                                                subj_nodes.push(d);
                                             })
                                
                                 selected_node.forEach(function(key, value)
                                                          {
                                                                var val_dta = d3.select(value).datum();
                                                                query_nodes.push(val_dta);   
                                                          });
                                    
                                    copy_nodes(query_nodes, subj_nodes, g_ind);
                                //
                                //if (subj_nodes.length > 0)
                                //{
                                //
                                //    selected_node.forEach(function(key, value)
                                //                          {
                                //                                var val_dta = d3.select(value).datum();
                                //                                query_nodes.push(val_dta);   
                                //                          });
                                //    
                                //    copy_nodes(query_nodes, subj_nodes, g_ind);
                                //}
                                //else
                                //{
                                // 
                                //    selected_node.forEach(function(key, value)
                                //         {
                                //           d3.select(value).datum().fixed=false;
                                //         });
                                //
                                //    node_tick(force2, node, back_node, anchorNode, anchorLink, link, false);
                                //    force.resume();
                                //}
                            }
                            else
                            {
                              
                                selected_node.forEach(function(key, value)
                                         {
                                           d3.select(value).datum().fixed=false;
                                         });
                   
                                node_tick(force2, node, back_node, anchorNode, anchorLink, link, false);
                                force.resume();
                            }
                           
                        }
                    }
                    else
                    {
                        if (in_current_box(vis, w, h) == false)
                        {
                            selected_node.forEach(function(key, value)
                                         {
                                           d3.select(value).datum().fixed=false;
                                         });
                   
                            node_tick(force2, node, back_node, anchorNode, anchorLink, link, false);
                            force.resume();
                        }
                        
                    }
                    
                        
                }
            }
            //current drag mechanism
            //function dragstart(d) {
            //    d.fixed = true;
            //    d3.select(this).classed("fixed", true);
            //  }
            //       
            //var drag = force.drag()
            //    .on("dragstart", dragstart);
            //node is not defined yet otherwise... node.call(drag);
            
            //remove_menu_items("graph_att_div");
            //
            //var graph_atts = document.getElementById("graph_att_div");
            //graph_atts.appendChild(make_text_div("Total Node Count: " + graph_obj.nodes.length, false, "black", "220px"));
            //graph_atts.appendChild(make_text_div("Total Edge Count: " + graph_obj.links.length, false, "black", "220px"));
            //graph_atts.appendChild(make_text_div("Node Types:", false, "black", "220px"));
            
            
            //remove the labels if they exist...is there a better way?
            vis.selectAll("g.labelNode").data([]).exit().remove();
            vis.selectAll("path.link").data([]).exit().remove();
            vis.selectAll("g.node").data([]).exit().remove();
            vis.selectAll("g.back_node").data([]).exit().remove();
            
            graph_obj.nodes = get_node_size(graph_obj.nodes);
            
            //initializing the node positions
            for(i=0;i < graph_obj.nodes.length;i++)
            {
                if (d3.map(graph_obj.nodes[i]).has('x') == false)
                {
                    graph_obj.nodes[i].x=Math.floor(Math.random() * w)
                    graph_obj.nodes[i].y=Math.floor(Math.random() * h)
                }
            }
            
            force.nodes(graph_obj.nodes)
                .links(graph_obj.links)
                .start();
            
            var labelNodes = [];
            var labelEdges = [];
            
            //for (i = labelEdges.length; i < graph_obj.nodes.length; i++)
            for (i = 0; i < graph_obj.nodes.length; i++)
            {
                //twice for some reason...
                labelNodes.push({'id': (i + ".1"), 'node':graph_obj.nodes[i]});
                labelNodes.push({'id': (i + ".2"), 'node':graph_obj.nodes[i]});
                labelEdges.push({'source':i*2, 'target':i*2+1});
            }
            
            force2.nodes(labelNodes)
                .links(labelEdges)
                .start();
            
            
            
            var edge_map = {}
            
            //add to each link the total number of links per edge and its current position
            graph_obj.links.forEach(function(el, ind, ar)
                                    {
                                        var cur_ind = el.source.id+'.'+el.target.id;
                                        if (cur_ind in edge_map)
                                        {
                                            edge_map[cur_ind] += 1;
                                        }
                                        else
                                        {
                                            edge_map[cur_ind] = 1;
                                        }
                                        
                                        graph_obj.links[ind].cur_edge = edge_map[cur_ind];
                                    })
            
            graph_obj.links.forEach(function(el, ind, ar)
                                    {
                                        var cur_ind = el.source.id+'.'+el.target.id;
                                        
                                        graph_obj.links[ind].edge_total = edge_map[cur_ind];
                                        
                                    })
            
            //link = vis.selectAll("line.link")
            //.data(graph_obj.links).enter().append("svg:line");
            
            link = vis.selectAll("path.link").data(graph_obj.links).enter().append("path");
            
            link.attr("class", function(d,i)
                  {
                        return ("link " + type_translation(d.attributes.type, edge_transl));
                  })
            .style("stroke-width", function(d,i)
                   {
                        return(get_line_size(d.attributes).size);
                   });
            
            //back_node makes it so that there is the illusion that the edges don't go to the center of the circle...
            //as these nodes are never made transparent.
            back_node = vis.selectAll("g.back_node")
            .data(force.nodes(), function(d){return d.id;})
            back_node.enter().append("svg:g").attr("class", "back_node");
            back_node.append("svg:circle")
            .attr("r", function(d) { return d.size|| 10; })
            .attr("fill", "#FFFFFF")
            
            node = vis.selectAll("g.node")
            .data(force.nodes(), function(d){return d.id;})
            node.enter().append("svg:g").attr("class", "node");
            node.append("svg:circle")
            .attr("r", function(d) { return d.size|| 10; })
            .attr("class", function(d,i){return d.attributes.node_type;});
           
            node.each(function(p_node,i)
                      {
                        
                        var cur_size = [(p_node.size*2), (p_node.size*2)];
                        
                        var temp_p_node = JSON.parse(JSON.stringify(p_node));
                        
                        var node_hier_filt = jQuery.extend(true, {}, temp_p_node);
                        
                        var init_pack = d3.layout.pack().size(cur_size).value(function(d){return (1)}).nodes(node_hier_filt);
                        
                        var pack_scale = {};
                        
                        //this just enforces that the leaf nodes are somewhat comparable in size...
                        
                        for (i = 0; i < init_pack.length;i++)
                        {
                            if (init_pack[i].size)
                            {
                                pack_scale[init_pack[i].display_name] = init_pack[i].size;
                                
                            }
                            else if (init_pack[i].parent && init_pack[i].parent.size)
                            {
                                
                                pack_scale[init_pack[i].display_name] = init_pack[i].parent.size/init_pack[i].parent.children.length
                                
                                
                            }
                        }
                        
                        var use_pack = d3.layout.pack().size(cur_size).value(function(d)
                                                                             {
                                                                                //console.log(d.display_name)
                                                                                if (d.display_name in pack_scale)
                                                                                {
                                                                                    return(pack_scale[d.display_name]);
                                                                                }
                                                                                else
                                                                                {
                                                                                    return(1);
                                                                                }
                                                                                
                                                                             
                                                                             }).nodes(node_hier_filt);
                        
                        //remove the toplevel object...
                        
                        
                        filt_pack = use_pack.filter(function(d)
                                        {
                                             return((d.id != p_node.id) && (d.depth <2));
                                        });
                        
                        var p_node_obj = this;
                        
                        function recurse_child (cur_node, cur_g)
                        {
                            if (cur_node.children && cur_node.children.length > 0)
                            {
                                //for simplicity of determing position relative to the current g, just repack relative to the temp_pack keeping all children the same radius
                                
                                var temp_pack = d3.layout.pack().size([cur_node.r*2, cur_node.r*2]).radius(cur_node.children[0].r).nodes(cur_node);
                                
                                temp_pack.forEach(function(cur_d,ind)
                                                  {
                                                        if (cur_d.depth == 0)
                                                        {
                                                            if (cur_d.children.length == 1)
                                                            {
                                                                var g_leftoff = cur_g.selectAll("g").data(cur_d.children, function(d) {return(d.id)}).enter().append("g")
                                                                .attr("transform", function(d) {return("scale(.75) translate(" + ( d.x-cur_node.r) + "," + (d.y-cur_node.r) + ")");});
                                                                g_leftoff.append("circle")
                                                                        .attr("class", function(d) {return(type_translation(d.attributes.node_type));})
                                                                        .attr("r", function(d) {return(d.r);});
                                                                    
                                                                recurse_child(cur_d.children[0], g_leftoff);
                                                            }
                                                            else
                                                            {
                                                                for(i=0; i < cur_d.children.length;i++)
                                                                {
                                                                    var g_leftoff = cur_g.selectAll("g").data([cur_d.children[i]], function(d) {return(d.id)}).enter().append("g")
                                                                    .attr("transform", function(d) {return("translate(" + ( d.x-cur_node.r) + "," + (d.y-cur_node.r) + ")");});
                                                                    
                                                                    g_leftoff.append("circle")
                                                                        .attr("class", function(d) {return(type_translation(d.attributes.node_type));})
                                                                        .attr("r", function(d) {return(d.r);});
                                                                    
                                                                    recurse_child(cur_d.children[i], g_leftoff);
                                                                }
                                                            }
                                                            
                                                        }
                                                        
                                                  })
                          
                            }
                            
                        }
                        
                        if (filt_pack.length == 1)
                        {
                            var cur_g = d3.select(p_node_obj).selectAll("g").data(filt_pack, function(d) {return(d.id)}).enter().append("g").attr("transform", function(d){ return("scale(.75) translate(" + (d.x-p_node.size) + "," + (d.y-p_node.size) + ")")})
                                                cur_g.append("circle")
                                                                .attr("class", function(d) {return(type_translation(d.attributes.node_type))})
                                                                 .attr("r", function(d) {return(d.r)});
                            recurse_child(filt_pack[0], cur_g);
                        }
                        else
                        {
                            filt_pack.forEach(function(use_d,i)
                                          {
                                                //console.log(use_d.id)
                                                var cur_g = d3.select(p_node_obj).selectAll("g").data([use_d], function(d) {return(d.id)}).enter().append("g").attr("transform", function(d){ return("translate(" + (d.x-p_node.size) + "," + (d.y-p_node.size) + ")")})
                                                cur_g.append("circle")
                                                                .attr("class", function(d) {return(type_translation(d.attributes.node_type))})
                                                                 .attr("r", function(d) {return(d.r)})
                                                recurse_child(use_d, cur_g);
                                          });
                        }
                    
                    var outer_node = d3.select(p_node_obj).append("circle").attr("r", p_node.size).attr("opacity", 0);
                    
                
                    //from the blog:http://blog.safaribooksonline.com/2014/03/10/creating-right-click-contextual-popup-d3/
                    //outer_node.on("click", function()
                    //              {
                    //                    var use_data = d3.select(p_node_obj).datum();
                    //                    selected_node = emphasize_nodes(p_node_obj, use_data.id, selected_node, vis, true);
                    //              })
                    
                    outer_node.on("contextmenu", function()//was onclick
                                  {
                                    d3.event.stopPropagation();
                                    
                                    d3.event.preventDefault();
                                    if (contextMenuShowing) {
                                        d3.select(".popup").remove();
                                    } else {
                                        
                                        if (popover_ref != null)
                                        {
                                            $('[data-toggle="popover"]').popover('destroy');
                                            $(popover_ref).popover('destroy');
                                            popover_ref = null;
                                            //to take care of the leftover popover which doesn't get destroyed
                                            d3.selectAll('.popover').remove();
                                        }
                                        
                                        if (selected_node.keys().length == 0 || panel_context == 'image')
                                        {
                                          
                                            var use_data = d3.select(p_node_obj).datum();
                                            var post_node = {'nodes':JSON.stringify([use_data]), 'context':panel_context}
                                            
                                            $.post("/HitWalker2/node_query/", post_node, function(data, status, xhr)
                                               {
                                                    if (status == "success")
                                                    {
                                                        
                                                        //<div class="modal-header">
                                                        //placement code from: https://github.com/twbs/bootstrap/issues/1833
                                                        $(p_node_obj).popover({content:JSON.parse(data).content, title:use_data.attributes.node_type + ": " + use_data.display_name, container:"body", html:true, trigger:"manual",
                                                                              placement: 'right'});
                                                        $(p_node_obj).popover('show');
                                                        popover_ref = p_node_obj;
                                                        
                                                        
                                                        //necessary to initialize popovers for inclusion with html elements...http://stackoverflow.com/questions/18410922/bootstrap-3-0-popovers-and-tooltips
                                                        
                                                        $('[data-toggle="popover"]').popover({container:"body"});
                                                        
                                                        
                                                    }
                                                    else
                                                    {
                                                        console.log(status)
                                                    }
                                                
                                                    
                                               }, "text");
                                        }
                                        else
                                        {
                                            var post_list = []
                                            
                                            selected_node.forEach(function(key, value)
                                                                  {
                                                                        post_list.push(d3.select(value).datum());
                                                                  });
                                            
                                            var post_node = {'nodes':JSON.stringify(post_list)}
                                            
                                            $.post("/HitWalker2/multi_node_query/", post_node, function(data, status, xhr)
                                               {
                                                    if (status == "success")
                                                    {
                                                        //add a few controls to the content
                                                        
                                                        var post_content = JSON.parse(data).content;
                                                        
                                                        $(p_node_obj).popover({content:post_content, title:"Queries", container:"body", html:true, trigger:"manual"});
                                                        $(p_node_obj).popover('show');
                                                        popover_ref = p_node_obj;
                                                    }
                                                    else
                                                    {
                                                        console.log(status)
                                                    }
                                               }, "text");
                                        }
                                        
                                    }
                                    
                                  })
                    
                    outer_node.on("mouseover", function()
                                  {
                                    if (panel_context == 'panel')
                                    {
                                       current_object = true;
                                       selected_node = emphasize_nodes(p_node_obj, d3.select(this).datum().id, selected_node,  vis, false);
                                    }
                                    
                                  })
                    
                    outer_node.on("mouseout", function()
                                  {
                                    current_object = false;
                                    deemphasize_nodes(vis, selected_node, false);
                                    //$(p_node_obj).popover('hide');
                                  })
                    
                        
                      });
            //left off here...
            //vis.call(d3.behavior.zoom().x(x).y(y).scaleExtent([1, 1]).on("zoom", zoom));
            node.call(node_drag);//or just node.call(drag);
            
            //function zoom() {
            //    
            //    vis.selectAll("g.back_node").attr("transform", transform);
            //    vis.selectAll("g.node").attr("transform", transform);
            //    vis.selectAll("g.labelNode").attr("transform", transform);
            //    vis.selectAll("line.link").attr("x1", function(d) {
            //                                return x(d.source.x);
            //                        }).attr("y1", function(d) {
            //                                return y(d.source.y);
            //                        }).attr("x2", function(d) {
            //                                return x(d.target.x);
            //                        }).attr("y2", function(d) {
            //                                return y(d.target.y);
            //                        });
            //    
            //}

            function transform(d) {
              return "translate(" + x(d.x) + "," + y(d.y) + ")";
            }
            
            var anchorLink = vis.selectAll("line.labelEdges").data(labelEdges);
            anchorLink.exit().remove();
            //console.log(force2.nodes().length);
            
            anchorNode = vis.selectAll("g.labelNode").data(force2.nodes(), function(d) {return d.id;}).enter().append('svg:g').attr("class", "labelNode");
            anchorNode.append("svg:circle").attr("r", 0).style("fill", "#FFF");
            anchorNode.append("svg:rect").attr("height",8).style("visibility", function(d,i){
               return i % 2 == 0 ? "hidden" : "visible"
            }).style("fill", "#FFF");
            anchorNode.append("svg:text").text(function(d, i) {
                                return i % 2 == 0 ? "" : d.node.display_name //+ d.node.attributes.add_string
                        }).style("fill", "black")
                           .style("font-family", "Arial").style("font-size", "10px")
                           .style("pointer-events", "none");
            
            anchorNode.each(function(d,i)
                            {
                              var text_bounds = this.childNodes[2].getBBox();
                              d3.select(this.childNodes[1]).attr("width", text_bounds.width)
                                                            .attr("height", 8);
                              
                            })
            //anchorNode.selectAll("rect").attr("width", function(d,i){
            //   console.log(this);
            //   ;
            //}).
            
            //if (graph_obj.nodes.length > 20)
            //{
            //   anchorNode.selectAll("text").style("visibility", "hidden");
            //}
            
            //may make more sense to disable the below and turn the pointer-events off as above...
            //anchorNode.on("click", function()
            //              {
            //                emphasize_nodes(this, d3.select(this).datum().node.id, selected_node, vis, true);
            //                
            //                //var cur_name = d3.select(this).datum().node.id;
            //                ////was var temp_sel = 
            //                ////vis.selectAll("g.node").selectAll("circle").style("opacity", function(d,i)
            //                //vis.selectAll("g.node").each(function(d,i)
            //                //                            {
            //                //                                 if (d.id == cur_name)
            //                //                                 {
            //                //                                     node_click_highlight(this, vis);
            //                //                                     
            //                //                                     //group_pop_func();
            //                //                                 }
            //                //                                 
            //                //                            });
            //                //
            //                ////vis.selectAll("line.link").style("opacity", back_opacity);
            //                //vis.selectAll("line.link").classed("Unselected", true);
            //                ////vis.selectAll("g.labelNode").selectAll("text").style("opacity", function(d,i)
            //                //vis.selectAll("g.labelNode").selectAll("text").classed("Unselected", function(d,i)
            //                //                                                  {
            //                //                                                    if (d.node.id == cur_name)
            //                //                                                    {
            //                //                                                        return false;//1;
            //                //                                                    }
            //                //                                                    else
            //                //                                                    {
            //                //                                                        return true;//back_opacity;
            //                //                                                    }
            //                //                                                    });
            //                //
            //              })
            //
            //        .on("mouseover", function()
            //            {
            //                current_object = true;
            //                emphasize_nodes(this, d3.select(this).datum().node.id, selected_node, vis,false, shiftKey);
            //
            //            })
            //        .on("mouseout", function()
            //            {
            //                var cur_name = d3.select(this).datum().node.id;
            //                
            //                deemphasize_nodes(vis, selected_node);
            //               
            //                current_object = false;
            //            });
            
            //var subAnchorNode2 = vis.selectAll("g.labelNode.circle.text")
            //var subAnchorNode = vis.selectAll("g.labelNode.circle")//.data(force2.nodes(), function(d) {return d.id;});
            //var anchor_drag = force2.drag()
            //.on("dragstart", function(d)
            //    {
            //        d.fixed = true;
            //        d3.select(this).classed("fixed", true);
            //    });
            //anchorNode.call(anchor_drag);
            
            
            force.on("tick",function()
                     {
                        //console.log("hello")
                        node_tick(force2, node, back_node, anchorNode, anchorLink, link, false)
                     });
            
            //for (var i = 0; i < 1000; ++i)
            //{
            //   force.tick();
            //   force2.tick();
            //}
            //force.stop();
            //force2.stop();
            
            //allows cur_pos for each node to be updated according to what is available onscreen...
            //var temp_link_map = make_link_map(graph_obj.links);
            //
            //for (i = 0; i < graph_obj.nodes.length;i++)
            //{
            //    var path_names = d3.map();
            //    
            //    //reset to zero ahead of time as below increments
            //    for (j in graph_obj.nodes[i].attributes.paths)
            //    {
            //        graph_obj.nodes[i].attributes.paths[j].cur_pos = 0;
            //        
            //        if (path_names.has(graph_obj.nodes[i].attributes.paths[j].name))
            //        {
            //            console.log("ISSUE: duplicate name found");
            //        }
            //        
            //        path_names.set(graph_obj.nodes[i].attributes.paths[j].name, j);
            //    }
            //    
            //    var cur_links = temp_link_map.get(graph_obj.nodes[i].id);
            //    
            //    //console.log(JSON.stringify(path_names))
            //    //console.log(JSON.stringify(cur_links))
            //    if (cur_links != undefined)
            //    {
            //        //figure out what the classes of cur_links is and increment cur_pos attribute appropriately.
            //        cur_links.forEach(function(key, value)
            //        {
            //            if (value.source.id == key && path_names.has(value.source.attributes.node_type))
            //            {
            //                graph_obj.nodes[i].attributes.paths[path_names.get(value.source.attributes.node_type)].cur_pos += 1;
            //            }
            //            else if (value.target.id == key && path_names.has(value.target.attributes.node_type))
            //            {
            //                graph_obj.nodes[i].attributes.paths[path_names.get(value.target.attributes.node_type)].cur_pos += 1;
            //            }
            //        });
            //    }
            //    
            //}
            //
            //if (selected_node.keys().length == 1 && group_members.keys().length == 0)
            //{
            //    node_click(selected_node.values()[0]);
            //    
            //    var cur_data = d3.select(selected_node.values()[0]).datum();
            //
            //    var d_func = function(d,i)
            //    {
            //        if (d.id == cur_data.id)
            //        {
            //            return false;//1;
            //        }
            //        else
            //        {
            //            return true;//back_opacity;
            //        }
            //    }
            //    
            //    //vis.selectAll("line.link").style("opacity", back_opacity)
            //    vis.selectAll("line.link").classed("Unselected", true);
            //    //vis.selectAll("g.node").selectAll("circle").style("opacity", function(d,i)
            //    vis.selectAll("g.node").selectAll("circle").classed("Unselected", function(d,i)
            //    {
            //        if (d.id == cur_data.id)
            //        {
            //            return false;//1;
            //        }
            //        else
            //        {
            //            return true;//back_opacity;
            //        }
            //    });
            //    
            //    //vis.selectAll("g.labelNode").selectAll("text").style("opacity", function(d,i)
            //    vis.selectAll("g.labelNode").selectAll("text").classed("Unselected", function(d,i)
            //    {
            //        if (d.node.id == cur_data.id)
            //        {
            //            return false;//1;
            //        }
            //        else
            //        {
            //            return true;//back_opacity;
            //        }
            //    });
            //    
            //}
            
        };


//possibly useful functions

  //d3.select("#next_button").on("click", function()
  //       {
  //          move_in_stack(1);
  //       })
  //      
  //      d3.select("#prev_button").on("click",function()
  //          {
  //             move_in_stack(-1);
  //          })
  //      
  //      d3.select("#del_nc_button").on("click", function()
  //          {
  //              if (selected_node.keys().length > 0)
  //              {
  //                  var keep_links = [];
  //                  var keep_nodes = [];
  //                  var temp_nodes = d3.set();
  //                  
  //                  //keep the nodes that either are selected or are connected to a selected node, first by examining links
  //                  for (i=0; i < graph_obj.links.length;i++)
  //                  {
  //                      if (selected_node.has(graph_obj.links[i].source.id) == true || selected_node.has(graph_obj.links[i].target.id) == true)
  //                      {
  //                          keep_links.push(graph_obj.links[i])
  //                          if (temp_nodes.has(graph_obj.links[i].source.id) == false)
  //                          {
  //                              temp_nodes.add(graph_obj.links[i].source.id);
  //                          }
  //                          
  //                          if (temp_nodes.has(graph_obj.links[i].target.id) == false)
  //                          {
  //                              temp_nodes.add(graph_obj.links[i].target.id);
  //                          }
  //                      }
  //                  }
  //                  
  //                  for (i=0; i < graph_obj.nodes.length;i++)
  //                  {
  //                      if (temp_nodes.has(graph_obj.nodes[i].id) == true)
  //                      {
  //                          keep_nodes.push(graph_obj.nodes[i])
  //                      }
  //                  }
  //                  
  //                  graph_obj.nodes = keep_nodes;
  //                  graph_obj.links = keep_links;
  //                  
  //                  //remove_menu_items("node_att_div");
  //                  //remove_menu_items("filter_div");
  //                  //remove_menu_items("text_div");
  //                  
  //                  selected_node = d3.map();
  //                  
  //                  update_graph();
  //                  
  //                  graph_stack.push(JSON.stringify(graph_obj));
  //                  stack_pos = graph_stack.length - 1;
  //              }
  //          })
  //         
  //          d3.select("#del_sel_button").on("click", function()
  //          {
  //              if (selected_node.keys().length > 0)
  //              {  
  //                  var keep_links = [];
  //                  var keep_nodes = [];
  //                  
  //                  //delete edges if one of the nodes is selected
  //                  for (i=0; i < graph_obj.links.length;i++)
  //                  {
  //                      if (selected_node.has(graph_obj.links[i].source.id) == false && selected_node.has(graph_obj.links[i].target.id) == false)
  //                      {
  //                          keep_links.push(graph_obj.links[i])
  //                      }
  //                  }
  //                  
  //                  for (i=0; i < graph_obj.nodes.length;i++)
  //                  {
  //                      if (selected_node.has(graph_obj.nodes[i].id) == false)
  //                      {
  //                          keep_nodes.push(graph_obj.nodes[i])
  //                      }
  //                  }
  //                  
  //                  graph_obj.nodes = keep_nodes;
  //                  graph_obj.links = keep_links;
  //                  
  //                  //remove_menu_items("node_att_div");
  //                  //remove_menu_items("filter_div");
  //                  //remove_menu_items("text_div");
  //                  
  //                  selected_node = d3.map();
  //                  
  //                  update_graph();
  //                  
  //                  graph_stack.push(JSON.stringify(graph_obj));
  //                  stack_pos = graph_stack.length - 1;
  //                  
  //              }
  //          })

//
//

//

//
//
//var group_pop_func_real = function()
//{
//    if (selected_node.keys().length >= 1)
//        {
//            d3.select('#collapseThree').attr("class", "panel-collapse in").attr("style", "height:auto;")
//            
//            if (selected_node.keys().length == 1 && group_members.keys().length == 0)
//            {
//                node_click(selected_node.values()[0]);
//            }
//            
//            //remove_menu_items("group_div");
//            if (cur_group == 0)
//            {
//                var new_div = document.getElementById("group_div");
//                
//                var group_label_div = document.createElement("div");
//                group_label_div.id = "group_label";
//                //group_label_div.style.width = "220px";
//                
//                var group_choice_div = document.createElement("div");
//                group_choice_div.id = "group_choice";
//                //group_choice_div.style.width = "220px";
//                
//                var button_div = document.createElement("div");
//                button_div.className = "btn-group btn-group-justified";
//                
//                new_div.appendChild(group_label_div);
//                new_div.appendChild(group_choice_div);
//                new_div.appendChild(button_div);
//                
//                var unset_group_button = document.createElement('a');
//                unset_group_button.className = 'btn btn-default';
//                unset_group_button.innerHTML = 'Remove All';
//                unset_group_button.onclick = function()
//                {
//                    //remove the rectangles, clear group_members and initial_node_type;
//                    vis.selectAll("g.node_rect").data([]).exit().remove();
//                    group_members = d3.map();
//                    initial_node_type = null;
//                    cur_group = 0;
//                    //remove_menu_items("group_div");
//                }
//                
//                var set_group_button = document.createElement('a');
//                set_group_button.className = 'btn btn-default'
//                set_group_button.innerHTML = 'Set Group';
//                set_group_button.onclick = function()
//                                         {
//                                            set_groups();
//                                         }
//                
//                button_div.appendChild(set_group_button);
//                button_div.appendChild(unset_group_button);
//            }
//            
//            
//        }
//}



//probably depricated functions
//
//  var set_groups_old = function()
//        {
//            if (cur_group < 2)
//            {
//                //remove_menu_items("text_div");
//                
//                cur_group += 1
//                var valid_count = 0;
//                selected_node.forEach(function(key, value)
//                {
//                    var cur_data = d3.select(value).datum();
//                    if (cur_data.attributes.node_type == initial_node_type)
//                    {
//                        if (group_members.has(cur_data.id) == false)
//                        {
//                            group_members.set(cur_data.id, cur_group);
//                            valid_count += 1;
//                        }
//                        else
//                        {
//                            console.log("ISSUE: duplicate group found")
//                        }
//                    }
//                    else
//                    {
//                        //maybe de-highlight them and remove them from selected_nodes
//                    }
//                  
//                });
//                
//                if (valid_count > 0)
//                {
//                    var group_label_div = document.getElementById("group_label");
//                    
//                    d3.select("#group_label").append("div").attr("style", "color:" + group_colors[cur_group - 1]).text("Group " + cur_group + " Type: " + initial_node_type);
//                
//                    if (cur_group > 1)
//                    {
//                        vis.selectAll("g.node_rect").data([]).exit().remove();
//                    }
//                    
//                    
//                    group_rect = vis.selectAll("g.node_rect")
//                    .data(force.nodes(), function(d){return d.id;}).enter().insert("svg:g", "g.node").attr("class", "node_rect");
//                    group_rect.append("svg:rect").filter(function(d,i)
//                                                         {
//                                                            return group_members.has(d.id)
//                                                         })
//                    .attr("width", function(d) { return d.size*2 })
//                    .attr("height", function(d) { return d.size*2 })
//                    .attr("transform", function(d) {
//                            return "translate(" + -d.size + "," + -d.size + ")";
//                    })
//                    .attr("fill", "#FFFFFF");
//                    
//                     var get_rect_color = function(node_name)
//                        {
//                            var group_num = group_members.get(node_name);
//                            
//                            return group_colors[group_num - 1]
//                        }
//                    
//                    group_rect.attr("stroke", function(d,i)
//                                    {
//                                        return get_rect_color(d.id);
//                                    });
//                    
//                    
//                    group_rect.call(rectUpdateNode);
//                
//                var temp_select = document.getElementById("group_query");
//                
//                if (temp_select == null)
//                {
//                    temp_select = document.createElement('select');
//                    temp_select.className = "form-control"
//                    temp_select.id = "group_query";
//                    
//                    var all_query = document.createElement('option');
//                    all_query.innerHTML = 'Common Neighbors';
//                    all_query.value = "common";
//                    temp_select.appendChild(all_query);
//                    
//                    var group_exec_button = document.createElement('button');
//                    group_exec_button.className = "btn btn-default"
//                    group_exec_button.type = 'button';
//                    group_exec_button.innerHTML = 'Submit';
//                    group_exec_button.onclick = function()
//                    {
//                       
//                        var other_node_ids = [];
//                        
//                        for (i=0; i < graph_obj.nodes.length;i++)
//                        {
//                            if (group_members.has(graph_obj.nodes[i].id) == false)
//                            {
//                                other_node_ids.push(graph_obj.nodes[i].db_id);
//                            }
//                            
//                        }
//                        
//                        group_dict = {};
//                        var temp_node_map = make_node_map(graph_obj.nodes);
//                        
//                        group_members.forEach(function(key, value)
//                                              {
//                                                    if (value in group_dict)
//                                                    {
//                                                        group_dict[value].push(temp_node_map.get(key).db_id);
//                                                    }
//                                                    else
//                                                    {
//                                                        group_dict[value] = [temp_node_map.get(key).db_id];
//                                                    }
//                                              });
//                       
//                        var check_ar = document.getElementsByName("group_paths");
//                        
//                        var path_names = [];
//                        
//                        for(i=0; i < check_ar.length;i++)
//                        {
//                            if (check_ar[i].checked == true)
//                            {
//                                path_names.push(check_ar[i].value);
//                            }
//                        }
//                        
//                        var query_select = document.getElementById("group_query");
//                        
//                        var query_type = query_select.options[query_select.selectedIndex].value;
//                        
//                        Dajaxice.network.get_group_query(path_callback2, {'group_dict':group_dict, 'other_node_ids':other_node_ids, 'path_names':path_names, 'query_type':query_type, 'path_filters':{}})
//                    }
//                    
//                    //this assumes that the node picked has the necessary info, which it should...
//                    //may need a better way of doing this comprehensively
//                    var node_map = make_node_map(graph_obj.nodes)
//                    var temp_node = node_map.get(group_members.keys()[0])
//                    var path_trans = d3.map()
//                    for(i in temp_node.attributes.paths)
//                    {
//                        if (path_trans.has(temp_node.attributes.paths[i].link) == false)
//                        {
//                            path_trans.set(temp_node.attributes.paths[i].link, i)
//                        }
//                    }
//                    
//                    var valid_path_list = avail_group_paths[initial_node_type];
//                    var group_choice_div = document.getElementById("group_choice");
//                    
//                    for(i=0; i < valid_path_list.length; i++)
//                    {
//                        var temp_check_div = document.createElement('div');
//                        temp_check_div.className = "checkbox";
//                        
//                        if (color_map.has(path_trans.get(valid_path_list[i])))
//                        {
//                            use_color = color_map.get(path_trans.get(valid_path_list[i]));
//                        }
//                        else
//                        {
//                            use_color = "black"
//                        }
//                        
//                        temp_check_div.style.color = use_color;
//                        
//                        var temp_check_list = document.createElement('input');
//                        //temp_check_list.style.cssFloat="left"
//                        temp_check_list.type = "checkbox";
//                        temp_check_list.name = "group_paths";
//                        temp_check_list.value = valid_path_list[i];
//                        temp_check_list.style.width="10px"
//                        temp_check_list.checked=true;
//                        
//                        var temp_check_label = document.createElement('label');
//                        temp_check_label.innerHTML = path_trans.get(valid_path_list[i]);
//                        
//                        temp_check_label.appendChild(temp_check_list);
//                        temp_check_div.appendChild(temp_check_label);
//                        group_choice_div.appendChild(temp_check_div);
//                        //group_choice_div.appendChild(make_text_div(path_trans.get(valid_path_list[i]), false, use_color, "200px"));
//                    }
//                    
//                    group_choice_div.appendChild(temp_select);
//                    
//                    group_choice_div.appendChild(group_exec_button);
//                }
//                else
//                {
//                    for (i=1; i <= cur_group;i++)
//                    {
//                        var all_query = document.createElement('option');
//                        all_query.innerHTML = 'Unique to Group ' + i;
//                        all_query.value = "group_" + i;
//                        temp_select.appendChild(all_query);
//                    }
//                }
//                }
//                else
//                {
//                    cur_group -= 1;
//                }
//            }
//        }
//        
//        var group_pop_func = function()
//        {
//            
//        }
//        
//
// function path_callback2(data)
//        {
//            
//            var full_node_map = make_node_map(graph_obj.nodes);
//            var keep_node_map = d3.map();
//            var keep_nodes = [];
//            var keep_links = [];
//            
//            //add the new nodes to the keep_node_map
//            
//            for (i = 0; i < data.existing_node_list.length;i++)
//            {
//                if (full_node_map.has(data.existing_node_list[i]))
//                {
//                    keep_node_map.set(data.existing_node_list[i], full_node_map.get(data.existing_node_list[i]))
//                }
//                else
//                {
//                    console.log("ISSUE: Cannot find key: " +  data.existing_node_list[i])
//                }
//            }
//            
//            //if (keep_node_map.has(data.cur_node_name))
//            //{
//            //    var temp_node = keep_node_map.get(data.cur_node_name)
//            //    var temp_check_map = get_check_map(data.cur_node_name);
//            //    var temp_path_count = temp_node.attributes.paths;
//            //    for (i in data.cur_path_count)
//            //    {
//            //        temp_path_count[temp_check_map.get(i)] = data.cur_path_count[i]
//            //    }
//            //    temp_node.attributes.paths = temp_path_count;
//            //    keep_node_map.set(data.cur_node_name, temp_node)
//            //}
//            //else
//            //{
//            //    console.log("ISSUE: Cannot find key: " + data.cur_node_name)
//            //}
//            
//            for (i=0; i < graph_obj.links.length;i++)
//            {
//                if (keep_node_map.has(graph_obj.links[i].source.id) && keep_node_map.has(graph_obj.links[i].target.id))
//                {
//                    keep_links.push(graph_obj.links[i])
//                }
//            }
//            
//            for (i = 0; i < data.new_node_list.length; i++)
//            {
//                if (keep_node_map.has(data.new_node_list[i].id) == false)
//                {
//                    keep_node_map.set(data.new_node_list[i].id, data.new_node_list[i])
//                }
//                else
//                {
//                    //this is probably ok, at least for the keep_existing = false as those nodes may be retrieved again.
//                    console.log("ISSUE: Duplicate node error " + data.new_node_list[i].id)
//                }
//            }
//            
//            for (i = 0; i < data.link_list.length; i++)
//            {
//                keep_links.push({'source':keep_node_map.get(data.link_list[i].source), 'target':keep_node_map.get(data.link_list[i].target)})
//            }
//            
//            keep_nodes = keep_node_map.values()
//            
//            graph_obj.nodes = keep_nodes;
//            graph_obj.links = keep_links;
//            
//            update_graph();
//         
//            graph_stack.push(JSON.stringify(graph_obj));
//            stack_pos = graph_stack.length - 1;   
//        }
//
//function expand_graph(keep_first, node_name)
//        {
//            //get the checkboxes that were selected
//            var cur_checks = document.getElementsByName("exp_checks_" + node_name);
//            
//            var path_list = [];
//            
//            for(i=0; i < cur_checks.length;i++)
//            {
//                if (cur_checks[i].checked == true)
//                {
//                    path_list.push(cur_checks[i].value)
//                }
//            }
//            
//            if (path_list.length > 0 && selected_node.keys().length == 1)
//            {
//                get_path_nodes(d3.select(selected_node.values()[0]).datum(), path_list,keep_first)
//            }
//            
//        }
//        
//        function get_path_nodes (cur_data, path_list, keep_first)
//        {
//                
//                var other_node_ids = []
//                
//                var link_map = make_link_map (graph_obj.links)
//                var cur_db_links = link_map.get(cur_data.id)
//                
//                var drop_elem = document.getElementById("keep_existing");
//                var keep_existing_res = drop_elem.options[drop_elem.selectedIndex].value;
//                
//                for (i=0; i < graph_obj.nodes.length;i++)
//                {
//                    if (graph_obj.nodes[i].db_id != cur_data.db_id)
//                    {
//                        if ((keep_existing_res == "false" && cur_db_links.has(graph_obj.nodes[i].id)) || (keep_existing_res == "true"))
//                        {
//                            other_node_ids.push(graph_obj.nodes[i].db_id)
//                        }
//                    }
//                    
//                }
//                
//                var path_filters = {};
//                var path_counts = {};
//                var any_path_count_gt_0 = false;
//                
//                var check_map = get_check_map(cur_data.id);
//                for (i=0; i < path_list.length;i++)
//                {
//                    console.log(path_list[i])
//                    var temp_att = cur_data.attributes.paths[check_map.get(path_list[i])];
//                    console.log(console.log(JSON.stringify(temp_att)))
//                    console.log(path_list[i])
//                    path_counts[path_list[i]] = temp_att;
//                    //path_filters[path_list[i]] = get_edge_filters(path_list[i]);
//                    
//                    if (temp_att.cur_pos < temp_att.count)
//                    {
//                        any_path_count_gt_0 = true;
//                    }
//                }
//                
//                if (any_path_count_gt_0 == true)
//                {
//                    Dajaxice.network.get_path_nodes(path_callback2, {'cur_node_name':cur_data.id, 'cur_node_id':cur_data.db_id, 'other_node_ids':other_node_ids,
//                                                'path_names':path_list, 'path_filters':path_filters, 'keep_first':keep_first, 'path_counts':path_counts})
//                }
//        }
//        
// function get_check_map(node_name)
//        {
//            var all_checks = document.getElementsByName("exp_checks_" + node_name);
//            
//            var name_map = d3.map();
//            
//            for (i=0; i < all_checks.length; i++)
//            {
//                var cur_check = all_checks[i].id.split("#");
//                
//                if (name_map.has(cur_check[1]) == false)
//                {
//                    name_map.set(cur_check[1], cur_check[0])
//                }
//                else
//                {
//                    console.log("ISSUE: found a duplicate key/value for checkboxes");
//                }
//            }
//            
//            return(name_map);
//        }
//
//function node_click(node_obj)
//        {
//            
//        }
//        
//        function node_click_real(node_obj)
//        {
//            
//            cur_check_count = [];
//            
//            var cur_data = d3.select(node_obj).datum();
//            
//            //selected_node.push(node_obj);
//            //here, we perform an ajax call to determine the available paths...
//            
//            var data_map = d3.map(cur_data)
//           
//            //remove_menu_items("text_div");
//            
//            var new_div = document.getElementById("text_div");
//            
//            var node_atts = d3.map(data_map.get("attributes"));
//            
//            var node_type = node_atts.get("node_type");
//
//            if (node_type != undefined)
//            {
//                d3.select(new_div).append("div").text("Node Type: " + node_type)
//                //new_div.appendChild(make_text_div("Node Type: " + node_type, false, "black", "220px"));
//            }
//            
//            var node_paths = d3.map(node_atts.get('paths'))
//            
//            node_paths.forEach(function(key, value)
//                            {
//                                //add checkboxes before each text element
//                                var check_div = document.createElement('div');
//                                check_div.className="checkbox";
//                                
//                                var use_color = "black";
//                                    
//                                if (color_map.has(value.name))
//                                {
//                                    use_color = color_map.get(value.name);
//                                }
//                                
//                                d3.select(check_div).attr("style", "color:"+ use_color)
//                                
//                                var check_label = document.createElement('label');
//                                check_label.innerHTML = value.name + " Count: " + value.cur_pos + "/" + value.count;
//                                
//                                var temp_check = document.createElement('input');
//                                temp_check.type = "checkbox";
//                                temp_check.name = "exp_checks_" + data_map.get('id');
//                                temp_check.value = value.link;
//                                temp_check.style.cssFloat="left"
//                                temp_check.id = key + "#" + value.link
//                                temp_check.onclick=function() make_filters_checkbox(this, value.count-value.cur_pos);
//                                check_label.appendChild(temp_check);
//                                check_div.appendChild(check_label);
//                                
//                                new_div.appendChild(check_div);
//                            })
//            
//            var temp_select = document.createElement('select');
//            temp_select.id = "keep_existing";
//            temp_select.className="form-control"
//            
//            var all_con = document.createElement('option');
//            all_con.innerHTML = 'All Connections';
//            all_con.value = "true";
//            temp_select.appendChild(all_con);
//            
//            var dir_con = document.createElement('option');
//            dir_con.innerHTML = 'Direct Connections';
//            dir_con.value = "false";
//            temp_select.appendChild(dir_con);
//            
//            new_div.appendChild(temp_select);
//            
//            var button_div = document.createElement("div");
//            button_div.className = "btn-group btn-group-justified";
//            
//            //put a button here called expand
//            var temp_button = document.createElement('a');
//            temp_button.className = "btn btn-default";
//            temp_button.id = 'all_neighbors';
//            temp_button.onclick = function() expand_graph(false, data_map.get('id'));
//            temp_button.disabled = true;
//            temp_button.innerHTML = "All";
//            button_div.appendChild(temp_button);
//            
//            var temp_button2 = document.createElement('a');
//            temp_button2.className = "btn btn-default";
//            temp_button2.id = 'next_neighbor';
//            temp_button2.onclick = function() expand_graph(true, data_map.get('id'));
//            temp_button2.disabled = true;
//            temp_button2.innerHTML = "Next";
//            button_div.appendChild(temp_button2);
//            
//            new_div.appendChild(button_div);
//            
//            d3.select('#collapseFour').attr("class", "panel-collapse in").attr("style", "height:auto;")
//            
//            //remove_menu_items("node_att_div");
//            //remove_menu_items("filter_div");
//            
//            //var att_div = document.getElementById("node_att_div");
//            //
//            //var node_meta = d3.map(node_atts.get("meta"));
//            //
//            //node_meta.forEach(function(key, value)
//            //{
//            //    att_div.appendChild(make_text_div(key + ":" + value, false, "black", "220px"));
//            //})
//        }
//        
//function make_filters_checkbox(checks, remaining_count)
//        {
//            var check_name = checks.id.split("#")[0]
//            
//            if (checks.checked == true)
//            {
//                cur_check_count.push({'name':check_name, 'count':remaining_count});
//            }
//            else
//            {
//                var rm_pos = [];
//                for(i=0; i < cur_check_count.length;i++)
//                {
//                    if (cur_check_count[i].name == check_name)
//                    {
//                        rm_pos.push(i);
//                    }
//                }
//                
//                if (rm_pos.length > 1)
//                {
//                    console.log("ISSUE: rm_pos > 1");
//                }
//                else if (rm_pos.length == 1)
//                {
//                    cur_check_count.splice(rm_pos[0], 1);
//                }
//                
//                //remove the current filter
//                //remove_menu_items(check_name + "_filter_div");
//            }
//            
//            var all_ns = document.getElementById('all_neighbors');
//            var next_ns = document.getElementById('next_neighbor');
//            
//            if (cur_check_count.length > 0)
//            {
//                var any_gt_zero = false;
//                var any_gt_max = false;
//                //make sure to enable All Neighbors or Next Neighbor depending on whether ANY of the queries would return greater than max_query_size elements
//                for (i=0; i < cur_check_count.length;i++)
//                {
//                   if (cur_check_count[i].count > max_query_size)
//                   {
//                        any_gt_max = true;
//                   }
//                   
//                   if (cur_check_count[i].count > 0)
//                   {
//                        any_gt_zero = true;
//                   }
//                }
//                
//                if (any_gt_zero == true && any_gt_max == true)
//                {
//                    all_ns.disabled=true;
//                    next_ns.disabled=false;
//                }
//                else if (any_gt_zero == true && any_gt_max == false)
//                {
//                    all_ns.disabled=false;
//                    next_ns.disabled=false;
//                }
//                else
//                {
//                    all_ns.disabled=true;
//                    next_ns.disabled=true;
//                }
//            }
//            else
//            {
//                all_ns.disabled=true;
//                next_ns.disabled=true;
//            }
//            
//            
//        }