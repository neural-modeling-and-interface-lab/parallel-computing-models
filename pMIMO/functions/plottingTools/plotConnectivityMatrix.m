function [h, ratio, possible, M] = plotConnectivityMatrix(fol,varargin)
    [freq_sort,criteria,plTitle,noPlot,onlyOrder,saveSubFolder]=process_options(varargin,'freq_sort',0,'criteria','CV','plTitle',[],'noPLot',0,'onlyOrder',[],'saveSubFolder', fullfile('exampleData','exampleResult',fol));
    not_loaded_color = [.5 .5 .5];

    if strcmp(criteria,'BIC')
        stage = 1;
        titleCritStr = ' (BIC)';
    elseif strcmp(criteria,'CV')
        stage=3;
        titleCritStr = '';
    else
        error('Criteria needs to be either CV or BIC!')
    end
    label = [criteria 'SparsityStruct'];
    fit_command = ['get' label];

    stageLabel = label2stagelabel(label,stage);
    meta_exist = existMetaData(fol,stageLabel,'saveSubFolder',saveSubFolder);
    if ~meta_exist
        addToMetaData(fol,label,stage,fit_command,'overWrite',1,'saveSubFolder',saveSubFolder)
    end
    f = getMetaDataFile(fol,'saveSubFolder',saveSubFolder);
    load(f,stageLabel,'modelOptions','n_ySp','xCellArray','yCellArray','y_labels')
    sparseData = eval(stageLabel);
    nF = length(n_ySp);
    nIn = length(xCellArray);

    if freq_sort
        [sorted_ySp2, sorted_ySp_inds] = sort(n_ySp,'descend');
        sig_input_cells={sparseData(sorted_ySp_inds).sig_input_cells};
        n_sig_inputs = {sparseData(sorted_ySp_inds).n_sig_inputs};
        sig_groups_1 = {sparseData(sorted_ySp_inds).sig_groups_1};
        sig_groups_2 = {sparseData(sorted_ySp_inds).sig_groups_2};
        y_labels=y_labels(sorted_ySp_inds);
        yCellArray=yCellArray(sorted_ySp_inds);
    else
        sig_input_cells={sparseData.sig_input_cells};
        n_sig_inputs = {sparseData.n_sig_inputs};
        sig_groups_1 = {sparseData.sig_groups_1};
        sig_groups_2 = {sparseData.sig_groups_2};
    end


    for i=1:length(sig_groups_1)
        n_sig_inputs1{i} = length(sig_groups_1{i});
        n_sig_inputs2{i} = length(sig_groups_2{i});
        all_from_sep = sort(unique([sig_groups_1{i}' sig_groups_2{i}']));
        if isempty(all_from_sep) && isempty(sig_input_cells{i}') % Modified by Xiwei Dec.14, 2021
            disp('No significant group in this case!')
        elseif ~isequal(all_from_sep,sig_input_cells{i}')
            error('Order 1 and 2 sig are not computed properly!')
        end
    end

    if ~isempty(onlyOrder)
        switch onlyOrder
            case 1
                n_sig_inputs = n_sig_inputs1;
                sig_input_cells = sig_groups_1;
            case 2
                n_sig_inputs = n_sig_inputs2;
                sig_input_cells = sig_groups_2;
            otherwise
                error('only order 1 or 2 should be specified!')
        end
    end
    %check to make sure s1 and s2 make up s_all

    inp_directions = {xCellArray.direction};

    M = zeros(nF,nIn);
    for i=1:length(sig_input_cells)
%         if ~isempty(n_sig_inputs{i})
        if ~isempty(sig_input_cells{i}) % Modified by Xiwei Dec.14, 2021
            M(i,sig_input_cells{i}) = 1;
        else
%             M(i,:) = -1;
            M(i,:) = 0;
            disp([num2str(i) ' is empty!'])
        end
        out_directions{i} = yCellArray(i).direction;
    end

    if ~noPlot
        h=imagesc(M);
        colormap([0 0 0; 1 1 0])
    else
        h = [];
    end

    if  ischar(out_directions{1})
        out_ant= find(ismember(out_directions,'Ant'));
        in_ant = find(ismember(inp_directions,'Ant'));
        out_post= find(ismember(out_directions,'Post'));
        in_post = find(ismember(inp_directions,'Post'));
        [ratio.total, count.total, possible.total] = connectivity_percentage(M);

        if freq_sort
            xlabel('CA3 Cell Number')
            ylabel('CA1 Cell Number')
        else
            xlabel(['CA3 Cell Number ' num2str(min(in_ant)) '-' num2str(max(in_ant)) ' Anterior'])
            ylabel(['CA1 Cell Number ' num2str(min(out_ant)) '-' num2str(max(out_ant)) ' Anterior'])
            M_ant_ant = M(out_ant,in_ant);
            [ratio.ant_ant, count.ant_ant, possible.ant_ant] = connectivity_percentage(M_ant_ant);
            M_ant_post = M(out_post,in_ant);
            [ratio.ant_post, count.ant_post, possible.ant_post] = connectivity_percentage(M_ant_post);
            M_post_ant = M(out_ant,in_post);
            [ratio.post_ant, count.post_ant, possible.post_ant] = connectivity_percentage(M_post_ant);
            M_post_post = M(out_post,in_post);
            [ratio.post_post, count.post_post, possible.post_post] = connectivity_percentage(M_post_post);
        end
    end
    % catch
    %     warning('Ratio Data Not Calculated')
    % end
    if isempty(plTitle)
        title([fol ', Out Freq Sort = ' num2str(freq_sort) titleCritStr])
    else
        title(plTitle)
    end

    try
        set(gca,'ytick',1:length(yCellArray))
%         set(gca,'yticklabel',num2str([yCellArray.BRIndex]'))
        ylabs = [];
        for yTemp = 1:length(yCellArray)
            ylabs = [ylabs; yCellArray(yTemp).name(3:5)];
        end
        set(gca,'yticklabel',ylabs) % Modified by Xiwei for TMS MIMO

        set(gca,'xtick',1:length(xCellArray)+1)
%         xlabs = num2str([xCellArray.BRIndex]');
        xlabs = []; 
        set(gca,'xticklabel',xlabs); % Modified by Xiwei for TMS MIMO
    catch
        warning('Cell Index is not labelled, only works with NEV files!')
    end

end

function [ratio, count, possible] = connectivity_percentage(M)
    possible = length(M(:));
    count = length(find(M(:)));
    ratio = count/possible;
end





