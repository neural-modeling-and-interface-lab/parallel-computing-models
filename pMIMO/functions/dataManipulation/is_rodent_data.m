function is_hum = is_rodent_data(f)
%checks to see if file contains human data 
f_info=whos( '-file', f);
var_names_combined = [f_info.name];
%in the data that we are using so far, all of the DNMS data has 'wire' in
%the name of variables, while the human data does not!
is_hum = isempty(strfind(var_names_combined,'before')) && isempty(strfind(var_names_combined,'after'));