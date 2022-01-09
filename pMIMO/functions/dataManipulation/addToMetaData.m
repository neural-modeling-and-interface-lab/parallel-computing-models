function addToMetaData(fol,label,stage,fit_command,varargin)
[stageLabel, wildcard] = label2stagelabel(label,stage);
[verbose, overWrite, saveSubFolder] = process_options(varargin,'verbose',1,'overWrite',0,'saveSubFolder', fullfile('exampleData','exampleResult',fol));

%saveSubFolder = fullfile('..','results_human',fol);
fsave = fullfile(saveSubFolder,[fol 'metaData.mat']);
meta_exist = existMetaData(fol,stageLabel,'saveSubFolder',saveSubFolder);
if meta_exist && ~overWrite
    disp('Not Adding Data because label already exists!!!!')
else
    if meta_exist
        disp('Overwriting field in metadata!!!!')
    end
        
    load(fsave,'y_labels')
    if isempty(y_labels{1})
        findbyLabel = 0;
        load(fsave,'y_names')
    else
        findbyLabel = 1;
    end
    
    listing = dir([saveSubFolder filesep '*' wildcard '*.mat'] );
    non_empty = [listing.bytes]>150;
    listing = [listing(non_empty)];
    fileNames = {listing.name};
    for i=1:length(fileNames)
        fload = [saveSubFolder filesep fileNames{i}];
        t0=tic;
        load(fload);
        %find which index it belongs to!
        if findbyLabel
        save_index = find(ismember(y_labels,fit.yCell.label));
        else
            save_index = find(ismember(y_names,fit.yCell.name));
        end
        var_to_save(save_index) = fit.(fit_command);
        tE = toc(t0);
        if verbose
            disp(['Just loaded ' fload ', took ' num2str(tE) 's'])
        end
    end
    % we want a new variable in the metadata .mat with the name specified in label!!, we are doing this by making the label a struct fieldname 
    save_my_field.(stageLabel) = var_to_save;
    save(fsave,'-struct','save_my_field','-append')
end



