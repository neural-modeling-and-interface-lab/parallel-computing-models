function was_Nev = is_Converted_Nev(f)
%checks to see if file is a .mat converted from a nev file
f_info=whos( '-file', f);
var_names_combined = [f_info.name];
%MetaTags are added when NEV files are converted, therefor, if it has a MetaData field, we assume it is a converted NEV file.
was_Nev = ~isempty(strfind(var_names_combined,'MetaTags'));