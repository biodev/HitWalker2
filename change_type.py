import sys
import os
import re
import shutil

if len(sys.argv) != 2:
    raise Exception('Need to provide a single argument')

change_type = sys.argv[1]

file_change_dict = {
    'network/config.py':[{'current':'prog_type\s*=\s*""', 'sub':'prog_type = "$$type$$"'}],
}

def rep_func(x, dict_list):
    
    if isinstance(dict_list, list):
        for i in dict_list:
            x = re.sub(i['current'], i['sub'].replace('$$type$$', change_type), x)
    else:
        x = re.sub(dict_list['current'], dict_list['sub'].replace('$$type$$', change_type), x)
   
    return x

for key,value in file_change_dict.items():
    
    #adjust files
    
    r_file = open(key, 'r')
    file_lines = r_file.readlines()
    r_file.close()
    
    new_file_lines = map(rep_func,file_lines, value*len(file_lines))
    
    w_file = open(key+'.new', 'w')
    w_file.writelines(new_file_lines)
    w_file.close()
    
    os.rename(key, key+'.old')
    os.rename(key+'.new', key)
    
    #modify the static file directories
    
    new_path = '/var/www/hitwalker2_inst/'+ change_type
    
    if os.path.exists(new_path):
        shutil.rmtree(new_path)
        
    os.mkdir(new_path)
    
    os.rename('/var/www/hitwalker2_inst/static',  new_path + '/static')
    