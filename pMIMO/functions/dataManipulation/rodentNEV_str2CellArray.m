function cellArrayTemp = rodentNEV_str2CellArray(f)
% creates a cell array from a humanNEV formatted file
    f_info=whos( '-file', f);
    n_f = length(f_info);
    %cell_strs = who('n*');
    labels = {'region','LR','label','nSpikes','BRIndex','BRChannel','BRUnit','name','session','animal','convertedNEVFormat'};
    [pathstr,sess_name,~] = fileparts(f);
    dataDir = pathstr;
    i_cell = 1;
    for i=1:n_f
        vars = parseLabels(f_info(i));
        if ~isempty(vars)
        vars{9} = sess_name;
        vars{10} = 'human'; % original name: 'human' 
        vars{11} = 1;
        cellArrayTemp(i_cell)=cellOb(vars,labels, dataDir);
        i_cell=i_cell+1;
        end
    end
    function vars = parseLabels(cell_info)
        s = cell_info.name;
        CA_ind = strfind(s,'CA');
        if isempty(CA_ind)   %Here, we say that all variables that don't have 'CA' in the title or not cells, this could change in the future!!!!
            vars = [];
        else  %format is CA1_n0111_ch011_un1_L
            vars{1} = s(CA_ind:CA_ind+2);
            vars{2} = s(end);   %the last cell is
            vars{3} = s(1:9);
            vars{4} = cell_info.size(1);
            vars{5} = str2double(s(6:9));
            vars{6} = str2double(s(13:15));
            vars{7} = str2double(s(19));
            vars{8} = s;
        end
    end
end