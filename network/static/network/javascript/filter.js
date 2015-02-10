
    var remove_menu_items = function(menu_id)
        {
            //console.log(menu_id)
            var old_div = document.getElementById(menu_id);
            while (old_div.firstChild) {
                old_div.removeChild(old_div.firstChild);
            }
            
        }
        
    var ensure_ok_values = function(obj)
    {
        var modal_id = $(obj).attr("data-modal");
        
        var modal_alerts = $("#"+$(obj).attr("data-modal")+" > div.modal-dialog > div.modal-content > div.modal-body > div.alert")
        
        if ($(modal_alerts).length == 0)
        {
            $('#'+modal_id).modal('hide');
        }else{
            console.log("hello")
        }
    }
    
    var add_message = function(message_text)
    {
        remove_menu_items("message_text_div");
        d3.select("#message_text_div").append("p").append("small").text(message_text);
    }
    
    var reset_def_filters = function()
    {
        filter_defaults = reset_filts;
        remove_all();
        make_default_filters();
        
        for (i in html_defaults)
        {
            document.getElementById(i).value = html_defaults[i].default;
        }
        
       add_message("Parameters have been reset to default values");
        
    }
    
    var save_clicked = function()
    {
        d3.select('#save_select').selectAll("option")
            .data(param_names.values()).enter()
            .append("option").on("click", function(d,i)
                                {
                                   change_text(d);  
                                }).text(function(d,i) {return(d);});
    }
    
    var load_clicked = function()
    {
         d3.select('#load_select').selectAll("option")
            .data(param_names.values()).enter()
            .append("option").text(function(d,i) {return(d);});
    }
    
    var change_text = function(text_val)
    {
        d3.select("#save_name").attr("value", text_val);
    }
    
    var load_ok_clicked = function(prog_type)
    {
        var load_obj = document.getElementById("load_select");
        var sel_name = load_obj.options[load_obj.selectedIndex];
        
        if (sel_name != undefined)
        {
            
              $.post(prog_type+"/HitWalker2/load_parameters/", {load_name:sel_name.innerHTML},
                                        function(data, status, xhr)
                                           {
                                                if (status == "success")
                                                {
                                                    var use_data = JSON.parse(data)
                                                    
                                                    filter_defaults = JSON.parse(use_data.filters);
                                                    
                                                    remove_all();
                                                    make_default_filters();
                                                    
                                                    var parsed_params = JSON.parse(use_data.parameters);
                                                    
                                                    for (i in parsed_params)
                                                    {
                                                        document.getElementById(i).value = parsed_params[i]
                                                    }
                                                    
                                                    bad_value.push(false);
                                                    d3.select("#param_load_alert").remove();
                                                    
                                                    remove_menu_items("message_text_div");
                                                    add_message("Successfully loaded '"+ sel_name.innerHTML + "'");
                                            
                                                    $('#load_modal').modal('hide');
                                                }
                                                else
                                                {
                                                    console.log(status)
                                                }
                                            
                                                
                                           }, "text")
        }
        else{
            if (document.getElementById("param_load_alert") == null)
            {
                d3.select("#param_load_div").append("div").attr("class", "alert alert-danger").attr("id", "param_load_alert").text("Please enter a valid name"); 
            }
            
            bad_value.push(true);
        }
        
    }
    
    var save_ok_clicked = function(prog_type)
    {
        var save_name = document.getElementById("save_name");
        
        if (save_name.value != "")
        {
            var filter_hs = {};
            var param_hs = {};
            var param_ar = ['res_prob', 'zscore', 'gene_score', 'max_iter', 'conv_thresh', 'string_conf'];
            
            d3.selectAll("[name^=parent]").each(function(d,i)
                                                        {
                                                            if (this.name != undefined)
                                                            {
                                                                filter_hs[this.name] = this.value;
                                                            }
                                                           
                                                        });
            
            var var_logic = make_button_info();
            
            for (i in param_ar)
            {
                param_hs[param_ar[i]] = document.getElementById(param_ar[i]).value
            }
            
            $.post(prog_type+"/HitWalker2/save_parameters/", {save_name:save_name.value, filter_hs:JSON.stringify(filter_hs), param_hs:JSON.stringify(param_hs), var_logic:JSON.stringify(var_logic)},
                                        function(data, status, xhr)
                                           {
                                                if (status == "success")
                                                {
                                                    var use_data = JSON.parse(data);
                                                    //looks successful, add the new name to the selection box options
                                                    param_names.add(use_data.save_name);
                                                    bad_value.push(false);
                                                    d3.select("#param_save_alert").remove();
                                                    $('#save_modal').modal('hide');
                                                    add_message("Successfully Saved '"+ use_data.save_name + "'");
                                                }
                                                else
                                                {
                                                    console.log(status)
                                                }
                                            
                                                
                                           }, "text")
        }
        else
        {
            if (document.getElementById("param_save_alert") == null)
            {
                d3.select("#param_save_div").append("div").attr("class", "alert alert-danger").attr("id", "param_save_alert").text("Please enter a valid name"); 
            }
            
            bad_value.push(true);
        }
    }
    
    var make_button_info = function()
    {
        var logic_buttons = document.getElementsByName("logical_button");
        
        button_vals = {'logical_button':{}, 'comp_button':{}}
        
        for (i=0; i < logic_buttons.length;i++)
        {
            button_vals["logical_button"][logic_buttons[i].id.split("_")[2]] = d3.select(logic_buttons[i]).text();
        }
        
        var comp_buttons = document.getElementsByName("comp_button");
        
        for (i=0; i < comp_buttons.length; i++)
        {
            button_vals["comp_button"][comp_buttons[i].id.split("_")[2]] = d3.select(comp_buttons[i]).text();
        }
        
       // console.log(JSON.stringify(button_vals));
        
        return(button_vals);
    }
    
    var check_all_fields = function()
    {
        
        if (document.getElementsByName("sample_alias")[0].value == "")
        {
            if (document.getElementById("sample_alias_alert") == null)
            {
                d3.select("#sample_alias_form_group").append("div").attr("class", "alert alert-danger").attr("id", "sample_alias_alert").text("Invalid sample specified"); 
            }
            
            bad_value.push(true);
        }
        else
        {
            d3.select("#sample_alias_alert").remove();
            bad_value.push(false);
        }
        
        $("input:not([name=sample_alias])").change();
        
        var any_bad = false;
        
        for (i=0; i < bad_value.length;i++)
        {
            if (bad_value[i] == true)
            {
                any_bad = true;
            }
        }
        
        bad_value = [];
        
        return any_bad == false;
    }
    
    var remove_all = function()
    {
        
        d3.selectAll("[name=parent_group_parent]").remove();
        
    }
    
    var make_default_filters = function()
    {
        filter_num = 0;
        cur_group = 0;
        parent_group = 0;
        bad_value = [];
        
        for(i=0; i < filter_defaults.length;i++)
        {
            for(j=0; j < filter_defaults[i].length; j++)
            {
                for (k=0; k < filter_defaults[i][j].length; k++)
                {
                    var use_defaults = ""
                    if ('default' in filter_defaults[i][j][k] == false)
                    {
                        use_defaults = filter_list.get('fields')[filter_defaults[i][j][k].field].default
                    }
                    else{
                        use_defaults = filter_num,filter_defaults[i][j][k].default
                    }
                    
                    basic_filter(filter_defaults[i][j][k].field,
                                 use_defaults,
                                 filter_defaults[i][j][k].comparison,
                                 filter_defaults[i][j][k].logical)
                    filter_num+=1;
                }
            
                if (j < (filter_defaults[i].length-1))
                {
                     select_group("subgroup");
                }
               
            }
            
            if (i < (filter_defaults.length-1))
            {
                select_group("group");
            }
            
        }
    }
    
    var select_group = function(type)
    {
        
        if (type == "group")
        {
            parent_group += 1;
            cur_group+=1;
        }
        else if (type == "subgroup")
        {
            cur_group+=1;
        }
        else
        {
            console.log("ISSUE: Unexpected type found")
        }
    }
    
    var make_filter_menu = function()
    {
        var filter_box = document.getElementById("filter_button");
        
        if (document.getElementById("filter_button_ul") != null)
        {
            var filter_ul = d3.select("#filter_button_ul")
        }
        else
        {
            var filter_ul = d3.select(filter_box).append("ul").attr("class", "dropdown-menu").attr("role", "menu").attr("id", "filter_button_ul");
        }
        
        var filter_map = []
        
        d3.map(filter_list.get("fields")).forEach(function(key, value)
                       {
                            filter_map.push([{'id':key, 'name':value.name}]);
                       })
        
        var filter_li = filter_ul.selectAll("li").data(filter_map).enter().append("li");
        
        var filter_a = filter_li.selectAll("a").data(function(d,i){return d;}).enter().append("a").text(function(d,i) {return d.name;});
        
        filter_a.attr("style", "cursor:pointer");
        filter_a.on("click", function(d,i)
                    {
                        basic_filter(d.id, filter_num, d3.map(filter_list.get("fields")).get(d.id).range[0], '&gt;', 'AND');
                        filter_num += 1;
                    });
    }
    
    
    var make_logical_button = function(filter_num, default_text, default_logical)
    {
        var before_logical_div =  d3.select("#group_" + cur_group).append("div").attr("class", "form-group").attr("id", "logical_div_" + filter_num)
                .attr("name", "var_logical_divs");
        var logical_div = before_logical_div.append("div").attr("class", "btn-group col-md-5");
            
        logical_div.append("button").attr("type", "button").attr("class", "btn btn-default dropdown-toggle").attr("data-toggle", "dropdown").attr("id", "logical_button_" + filter_num)
            .text(default_text).attr("name", "logical_button")
            .append("span").attr("class", "caret");
            
        var temp_log_li = logical_div.append("ul").attr("class", "dropdown-menu").attr("role", "menu")
            .selectAll("li").data([["AND"],["OR"]]).enter().append("li")
            .selectAll("a").data(function(d,i){return d;}).enter().append("a").text(function(d,i) {return d;})
            .attr("style", "cursor:pointer");
            
            temp_log_li.on("click", function(d,i)
                           {
                                d3.select("#logical_button_" + filter_num).text(d).append("span").attr("class", "caret");
                           });
    }
    
    var check_groups = function()
    {
        if (document.getElementById("parent_group_"+parent_group) == null)
        {
             var temp_group_div = d3.select("#grouped_fields").insert("div", "#filter_main")
                .attr("id", "parent_group_"+parent_group+"_parent").attr("class", "panel panel-default").attr("name", "parent_group_parent")
                .append("div").attr("id", "parent_group_"+parent_group).attr("name", "parent_group_div").attr("class", "panel-body");
        }
        
        if (document.getElementById("group_"+cur_group) == null)
        {
            var temp_group_div = d3.select("#parent_group_"+parent_group).append("div")
                .attr("id", "group_"+cur_group+"_parent").attr("class", "panel panel-default")
                .append("div").attr("id", "group_"+cur_group).attr("name", "group_div").attr("class", "panel-body");
        }
    }
    
    var basic_filter = function(filt_id, filter_num, default_value, comp_value, default_logical)
    {   
        //still need to figure out how to uniquely identify elements of the filters...
        
        var cur_filts = d3.map(filter_list.get("fields"));
        
        console.log(cur_filts)
        
        check_groups();
        
        if (document.getElementsByName("var_main_divs").length > 0)
        {
            make_logical_button(filter_num, default_logical);
           
           //implies that a new group has been added, add the button to a new div ahead of the new one by moving to the next group 
            if (d3.select("#group_"+ cur_group).node().childNodes.length == 1)
            {
                if (d3.select("#parent_group_"+ parent_group).node().childNodes.length == 1)
                {
                    d3.select("#parent_group_"+parent_group+"_parent").classed("panel panel-default", false);
                    d3.select("#parent_group_"+parent_group).classed("panel-body", false);
                    
                    parent_group += 1;
                }
                
                d3.select("#group_"+cur_group+"_parent").classed("panel panel-default", false);
                d3.select("#group_"+cur_group).classed("panel-body", false);
                
                cur_group += 1;
                
                check_groups();
            }
            
        }
        
        var before_div = d3.select("#group_"+cur_group).append("div").attr("class", "form-group").attr("id", "main_div_"+filter_num).attr("name", "var_main_divs");
        
        before_div.append("label").attr("class", "col-md-4 control-label").attr("for", "var_input_" + filter_num).text(cur_filts.get(filt_id).name);
        
        var inp_div = before_div.append("div");
        
        if (cur_filts.get(filt_id).type == "numeric")
        {
            inp_div.attr("class", "input-group col-md-5");
            var temp_button_div = inp_div.append("div").attr("class", "input-group-btn")
            temp_button_div.append("button").attr("class", "btn btn-default dropdown-toggle").attr("data-toggle", "dropdown")
                .attr("id", "var_button_" + filter_num).attr("name", "comp_button").text(comp_dict.get(comp_value));
            var temp_li = temp_button_div.append("ul").attr("class", "dropdown-menu").selectAll("li").data([["<"],["="], [">"]]).enter().append("li");
            temp_li.selectAll("a").data(function(d,i){return d;}).enter().append("a").text(function(d,i) {return d;});
            temp_li.attr("style", "cursor:pointer");
            
            temp_li.on("click", function(d,i)
                       {
                            d3.select("#var_button_" + filter_num).text(d);
                       })
            
            inp_div.append("input").attr("class", "form-control").attr("id", "parent." + parent_group + ".group." + cur_group + ".type." + filt_id + ".num." + filter_num)
                .attr("value", default_value)
                .on("change", function(d,i)
                    {
                        check_numeric_input(this, filt_id, "main_div_"+filter_num)
                    })
                .attr("name", "parent." + parent_group + ".group." + cur_group + ".type." + filt_id + ".num." + filter_num);
            
        }
        else if (cur_filts.get(filt_id).type == "character")
        {
            //inp_div.append("span").attr("class", "input-group-addon").text("=")
            inp_div.attr("class", "input-group col-md-5");
            var temp_sel = inp_div.append("select").attr("class", "form-control").attr("multiple", "multiple").attr("id", "parent." + parent_group + ".group." + cur_group + ".type." + filt_id + ".num." + filter_num)
            .attr("name", "parent." + parent_group + ".group." + cur_group + ".type." + filt_id + ".num." + filter_num)
            var cur_sel_opt = temp_sel.selectAll("option").data(cur_filts.get(filt_id).range).enter().append("option").text(function(d,i) {return d;});
            cur_sel_opt.filter(function(d,i)
                               {
                                    if (d == default_value)
                                    {
                                        return true;
                                    }
                                    else
                                    {
                                        return false;
                                    }
                               }).attr("selected", "selected")
        }
        else
        {
            console.log("ISSUE: unexpected filter type")
        }
        
        var close_box = inp_div.append("span").attr("class", "input-group-addon glyphicon glyphicon-minus-sign").attr("style", "cursor:pointer");
        
        close_box.attr("onclick", "delete_click(" + filter_num + ")")
        
    }
    
    var check_numeric_input_basic = function(cur_obj)
    {
        var valid_range = $.map($(cur_obj).attr("data-range").split("to"), function(n) {return parseFloat(n);})
        
        var modal_parent = $(cur_obj).attr("data-modal");
        
        if (jQuery.isNumeric($(cur_obj).val()) == false || parseFloat($(cur_obj).val()) < valid_range[0] || parseFloat($(cur_obj).val()) > valid_range[1])
        {
            console.log("grrr")
            if (document.getElementById($(cur_obj).attr("id") + "_alert") == null)
            {
                
                var assoc_label = $("label[for='"+$(cur_obj).attr("id")+"']").html();
                
                $("#"+modal_parent).find("div.modal-body")
                        .append("<div id="+$(cur_obj).attr("id") + "_alert class='alert alert-danger'><p>Invalid entry for "+assoc_label.replace(":", "")+" needs to be numeric between " + valid_range[0] + "-" + valid_range[1]+"</p></div>");
            }
            
            bad_value.push(true);
        }
        else
        {
            $("#" + $(cur_obj).attr("id") + "_alert").remove();
           
            bad_value.push(false);
        }
    }
    
    var check_numeric_input = function(cur_obj, filt_id, alert_div)
    {
        var exist_filt = d3.map(filter_list.get("fields")).get(filt_id);
        
        if (jQuery.isNumeric(cur_obj.value) == false || cur_obj.value < exist_filt.range[0] || cur_obj.value > exist_filt.range[1])
        {
            if (document.getElementById(alert_div + "_alert") == null)
            {
                d3.select("#" + alert_div).append("div").attr("class", "col-md-11 alert alert-danger").attr("id", alert_div + "_alert").text("Invalid entry, needs to be numeric between " + exist_filt.range[0] + "-" + exist_filt.range[1]); 
            }
            
            bad_value.push(true);
        }
        else
        {
            d3.select("#" + alert_div + "_alert").remove();
            bad_value.push(false);
        }
    }
    
    var delete_click = function(filter_num)
    {
        var cur_filter = 0;
        
        for(i=0; i < filter_defaults.length;i++)
        {
            for(j=0; j < filter_defaults[i].length; j++)
            {
                for (k=0; k < filter_defaults[i][j].length; k++)
                {
                    if (filter_num == cur_filter)
                    {
                        filter_defaults[i][j].splice(k, 1);
                    }
                    cur_filter+=1;  
                }
                
            }
            
        }
        remove_all();
        
        make_default_filters();
    }
    
    var delete_click_dep = function(filter_num)
    {
        //figure out whether the enclosing group is empty, if it is remove it as well...
        
        var cur_node = d3.select("#main_div_" + filter_num);
        var cur_parent = cur_node.node().parentNode;
        
        var previous_logic_exists = cur_node.node().previousSibling;
        var next_logic_exists = cur_node.node().nextSibling;
        
        if (previous_logic_exists == null && next_logic_exists == null)
        {
            //for some reason the prev_parent code and below doesn't work when there is only one group and one entry, the if statement addresses that
                //seems to be an additional unexpected node that is the firstChild
                //not sure if the above is true after the addition of the upper grouping layer...
                //changed prev_parent.parentNode.childNodes[0] to prev_parent.parentNode.childNodes[1] after the upper grouping layer
            
            var prev_parent = cur_parent.parentNode;
           
            if (prev_parent.parentNode.childNodes[0] != prev_parent)
            {
                 //also check to see if the previous sibling of the grandparent (group_x_parent) had a logical grandchild that needs to be removed
                
                var prev_sibl = prev_parent.previousSibling;
            
                if (prev_sibl != null)
                {
                    var prev_sibl_child = prev_sibl.firstChild.firstChild;
                    
                    if (prev_sibl_child != null && d3.select(prev_sibl_child).attr("name") == "var_logical_divs")
                    {
                        d3.select(prev_sibl).remove();
                    }
                }
            }
            else
            {
                //and also remove a logical div if it is the next sibling...
                 
                var next_sibl = prev_parent.nextSibling;
                
                if (next_sibl != null)
                {
                    var next_sibl_child = next_sibl.firstChild.firstChild;
                    
                    if (next_sibl_child != null && d3.select(next_sibl_child).attr("name") == "var_logical_divs")
                    {
                        d3.select(next_sibl).remove();
                    }
                }
            }
            
            //then remove the entire group as there is only an input div group
            d3.select(cur_parent.parentNode).remove();
            
        }
        else if (next_logic_exists == null && previous_logic_exists != null)
        {
            d3.select(previous_logic_exists).remove();
        }
        else if (next_logic_exists != null && previous_logic_exists == null)
        {
            d3.select(next_logic_exists).remove();
        }
        else
        {
            //if there are both, then remove the next div...
            d3.select(next_logic_exists).remove();
        }
        
        cur_node.remove();
        
        var all_parents = document.getElementsByName("parent_group_div")
        //still working on this...
        d3.selectAll("[name=parent_group_parent]").filter(function(d,i)
                    {
                         var prev_sibl = this.previousSibling;
                         var next_sibl = this.nextSibling;
                         
                         if (this.childNodes[0].childNodes.length == 0)
                         {
                             return true;
                         }
                         
                         if (prev_sibl.childNodes.length > 0)
                         {
                            if (prev_sibl.childNodes[0].childNodes.length == 0  && this.childNodes[0].childNodes.length == 1 && this.childNodes[0].childNodes[0].childNodes[0].childNodes[0].getAttribute("name") == "var_logical_divs")
                            {   
                                return true;
                            }
                            
                         }
                         
                         if (next_sibl.childNodes.length > 0)
                         {  
                            if (next_sibl.nextSibling.id == "filter_main" && next_sibl.childNodes[0].childNodes.length == 0 && this.childNodes[0].childNodes.length == 1 && this.childNodes[0].childNodes[0].childNodes[0].childNodes[0].getAttribute("name") == "var_logical_divs")
                            {
                                return true;
                            }
                            
                         }
                         
                         return false;
                         
                    }).remove();
    }
    
    var recurse_nodes = function (cur_obj, is_logical)
    {
        if (cur_obj.childNodes.length > 0)
        {
            return recurse_nodes(cur_obj.childNodes[0], cur_obj.name=="var_logical_divs" || is_logical)
        }
        else
        {
            return is_logical
        }
    }