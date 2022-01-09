function val = isNex(fname)
%checks to see if given file is a .nex file
[~,~,ext] = fileparts(fname);
val = strcmp(ext,'.nex');