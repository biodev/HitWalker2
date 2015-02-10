import sys
import os

if len(sys.argv) != 2:
    raise Exception('Need to provide a single argument')

change_type = sys.argv[1]

file_change_dict = {
    'HitWalker2/settings.py':[{"current":"'/var/www/hitwalker_2_inst/hit.db'", "sub":"'/var/www/hitwalker_2_inst_$$type$$/hit.db'"}],
    'HitWalker2/wsgi.py':[{"current":"'/var/www/hitwalker_2_inst/HitWalker2'", "sub":"'/var/www/hitwalker_2_inst_$$type$$/HitWalker2'"}]
}

def rep_func(x, dict_list):
    
    if isinstance(dict_list, list):
        for i in dict_list:
            x = x.replace(i['current'], i['sub'].replace('$$type$$', change_type))
    else:
        x = x.replace(dict_list['current'], dict_list['sub'].replace('$$type$$', change_type))
   
    return x

for key,value in file_change_dict.items():
    
    r_file = open(key, 'r')
    file_lines = r_file.readlines()
    r_file.close()
    
    new_file_lines = map(rep_func,file_lines, value*len(file_lines))
    
    w_file = open(key+'.new', 'w')
    w_file.writelines(new_file_lines)
    w_file.close()
    
    os.rename(key, key+'.old')
    os.rename(key+'.new', key)
    