function h = plotMetricComparison(fols,sessionLabels,metric,criterias,comparison,stages,varargin)
[FreqOrder, plHist,nBins,histLinws,hist_edges,labelY,noInf, saveSubFolder, globalYlim] = process_options(varargin,...
    'FreqOrder',1,'plHist',0,'nBins',10,'histLinws',ones(length(fols),1),'hist_edges',[],...
    'labelY',0,'noInf',1,'saveSubFolder',[], 'globalYlim', [-1 5]);
%noInf toggles whether or not we should use the bar_inf function, it seems that in Matlab 2015, there is a problem with the function
%comparison options: 'vs' 'none' OR one of the result types (e.g. 'zeroth')
%OR 'threshold' which plots bar charts with percent < threshold value which
%is default at one and should be used for KS plots!
switch metric
    case 'KS'
        label = 'KSStruct';
        fit_command = 'getKSStruct';
    case 'Loss'
        label = 'CVLossStrcut';
        fit_command = 'getCVLossStruct';
    otherwise
        disp('only Loss and KS metrics are availble')
end

for i_r = 1:length(fols)
    %load result data
    fol = fols{i_r};
    %if isempty(saveSubFolder)
    saveSubFolder = fullfile('exampleData','exampleResult',fol);
    % end
    if iscell(criterias)
        criteria=criterias{i_r};
    else
        criteria=criterias;
    end
    stage=stages(i_r);
    stageLabel = label2stagelabel(label,stage);
    meta_exist = existMetaData(fol,stageLabel,'saveSubFolder',saveSubFolder);
    if ~meta_exist
        addToMetaData(fol,label,stage,fit_command,'saveSubFolder',saveSubFolder)
    end
    f = getMetaDataFile(fol,'saveSubFolder',saveSubFolder);
    load(f,stageLabel,'n_ySp','y_labels')
    try
        load(f,stageLabel,'y_names')
    end
    metric_data =  eval(stageLabel);
    for i=1:length(metric_data)
        if ~isempty(metric_data(i).mu)
            if ~ismember(comparison,{'vs','none','threshold'})
                pl(i,i_r) = metric_data(i).mu.(criteria)-metric_data(i).mu.(comparison);
            else
                pl(i,i_r) = metric_data(i).mu.(criteria);
            end
            
        else
            pl(i,i_r) = nan;
        end
    end
end

if FreqOrder
    [~,pl_order]=sort(n_ySp,'descend');
    xL2 = ' (Freq Ordering)';
    pl=pl(pl_order,:);
else
    xL2 = ' (Natural Ordering)';
    
end


switch comparison
    case 'vs'
        data = pl(:,1)-pl(:,2);
        if ~plHist
            if ~noInf
                h = bar_inf(data);
            else
                h = bar(data);
            end
        else
            h = plotStairHist(data,hist_edges,'linws',histLinws);
        end
        tStr = [metric ', ' comparison ', ' sessionLabels{1} '-' sessionLabels{2}  ' (Neg Values mean first fit is better)'];
        
    case 'none'
        data = pl;
        if ~plHist
            if ~noInf
                h = bar_inf(data);
            else
                h = bar(data);
            end
            x = get(gca,'xlim');
            line(x,[1 1],'linestyle',':','color','k')
        else
            h = plotStairHist(data,hist_edges,'linws',histLinws);
        end
        tStr = [metric ', ' criteria];
        legend(sessionLabels)
        
    case 'threshold'  %this plots a bar graph with the number of values less than one
        data = sum(pl<1);
        h = bar(data);
        tStr = [metric ', ' comparison];
        ax = gca;
        ax.XTickLabel=sessionLabels;
        ylabel('count<1')
        
    otherwise
        data = pl;
        if ~plHist
            if ~noInf
                h = bar_inf(data);
            else
                h = bar(data);
            end
        else
            h = plotStairHist(data,hist_edges,'linws',histLinws);
        end
        tStr = [metric ', ' criteria '-' comparison];
        legend(sessionLabels)
end

if plHist
    xlabel(tStr);
    ylabel('Count');
elseif strcmp(comparison,'threshold')
    title(tStr)
else
    xlabel(['Output Cell' xL2]);
    title(tStr)
    if labelY && ~FreqOrder
        set(gca,'xtick',1:length(data))
        if ~isempty(y_labels{1})
            allYLabels = [y_labels];
        else
            allYLabels = [y_names];
        end
        if length(allYLabels)>10
            for i=1:length(allYLabels)
                allYLabels{i} = allYLabels{i}(3:5);
            end
        end
        ylim(globalYlim)
        set(gca,'xticklabel',allYLabels)
    end
end

if ~strcmp(comparison,'vs')
    
else
    
end
