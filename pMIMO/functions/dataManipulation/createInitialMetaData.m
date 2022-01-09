function createInitialMetaData(fol,varargin)

saveSubFolder = process_options(varargin,'saveSubFolder', fullfile('exampleData','exampleResult',fol));
fsave = fullfile(saveSubFolder,[fol 'metaData.mat']);
if ~exist(fsave,'file')
    Rs = loadAllFittingResults(saveSubFolder,'wildcard','Path','verbose',1);
    Rs = [Rs.fit];
    for i=1:length(Rs)
        n_ySp(i) = Rs(i).yCell.nSpikes;
        y_labels{i} = Rs(i).yCell.label;
        %if isempty(y_labels{i})
            y_names{i} = Rs(i).yCell.name;
        %end
        yCellArray(i) = Rs(i).yCell;
        fitTime(i)=Rs(i).fitTime;
        NYFit(i) = Rs(i).getNYFit;
        %reMLETime(i)=Rs(i).reMLETime(i);
        %CV_criteriaTime(i)=Rs(i).CV_criteriaTime(i);
    end
    modelOptions=Rs(1).modelOptions;
    regOptions= Rs(1).regOptions;
    xCellArray = Rs(1).xGroup.cellArray;
    fit_reMLE_CV = Rs(1).fit_reMLE_CV;
    fit_CV_criteria = Rs(1).fit_CV_criteria;
    fileLabel = Rs(1).fileLabel;
    T_trunc = Rs(1).T_trunc;
    for i=1:length(xCellArray)
        n_xSp(i) = xCellArray(i).nSpikes;
        x_labels{i} = xCellArray(i).label;
        %if isempty(x_labels{i})
            x_names{i} = xCellArray(i).name;
        %end
    end
    [sorted_ySp, sorted_ySp_inds] = sort(n_ySp,'descend');
    [sorted_xSp, sorted_xSp_inds] = sort(n_xSp,'descend');
    save(fsave,'n_ySp','n_xSp','y_labels','x_labels','yCellArray','xCellArray','fitTime','modelOptions','regOptions','fit_reMLE_CV','fit_CV_criteria','fileLabel','T_trunc','NYFit','y_names','x_names');
else
    disp('Metafile Already exists!')
end







