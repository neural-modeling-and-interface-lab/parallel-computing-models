%% This script creates a cell array from a DNMS formatted file
% Author: Xiwei She
% Last modified on: Oct. 9, 2021

function cellArrayTemp = DNMS_str2CellArray(f)
    f_info=whos( '-file', f);
    n_f = length(f_info);
    labels = {'id','region','wire','wire_cell','nSpikes','name','session','tend'};
    %if recording length is saved in .mat file, add it to the cell object
    if ismember('Stop',{f_info.name}) %add the end time of the trial if it is saved as Stop
        load(f,'Stop')
    else
        Stop = [];
    end
    [pathstr,sess_name,~] = fileparts(f);
    dataDir = pathstr;
    i_cell = 1;
    for i=1:n_f
        vars = parseLabels(f_info(i));
        if ~isempty(vars)
            vars{7} = sess_name;
            vars{8} = Stop;
            cellArrayTemp(i_cell)=cellOb(vars,labels, dataDir);
            i_cell = i_cell+1;
        end
    end
    function vars = parseLabels(cell_info)
        s = cell_info.name;
        wire_label_ind = strfind(s,'wire');  %%in all of the DNMS samples so far, there is the term wire in the name!
        if isempty(wire_label_ind)
            vars = [];
        else
            vars{1} = s(2:4);
            vars{2} = s(8:10);
            vars{3} = s(17);
            vars{4} = s(end);
            vars{5} = cell_info.size(1);
            vars{6} = s;
        end
    end
end