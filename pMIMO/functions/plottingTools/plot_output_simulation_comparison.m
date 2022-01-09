function [unstable_inds, unstable_labels, ax, y_real_filt, y_sim_filt] = plot_output_simulation_comparison(fs,varargin)
ax = [];
[outSort, sig_scale,max_map_scale,tmax,cellTogether,disp_color_bar,...
    extraLabel, plotInput,label_font_size,title_font_size,fSession,titleLabel,...
    sideLabel,xlims,dispLegend,correctLabel,new_figs,delete_fig,PAL, reBinSize,...
    clims, vert, labelAllEvents, labelCells,plotObservedOnly] ...
    = process_options(varargin, 'outSort',1,'sig_scale',1,'max_map_scale',.5,'tmax',100,...
    'cellTogether',0,'disp_color_bar',0,'extraLabel',1, 'plotInput',0,'label_font_size',10,...
    'title_font_size',10,'fSession',[],'titleLabel','','sideLabel',0,'xlims',[],...
    'dispLegend',1,'correctLabel',[],'new_figs',0,'delete_fig',0,'PAL',0,'reBinSize',[],...
    'clims',[],'vert',0,'labelAllEvents',0,'labelCells',0,'plotObservedOnly',0);
%default clims is empty which means that this function uses max_map_scale
%t define the limits of the output plotting.   if clims is deinfed, this
%will be used isntead.

%plotObservedOnly is a new option to only show observed Input and Output without
%plotting predicted Output.
%xGroup can be added to access x cell labels
figure()

if iscell(fs)
    Nf = length(fs);
else
    Nf=1;
end
for i_f = 1:Nf
    if iscell(fs)
        f = fs{i_f};
    else
        f = fs;
    end
    x_filt = [];
    y_real_filt=[];
    y_sim_filt = [];
    load(f)
    bin_size = .002;
    sigma = tmax*sig_scale;
    t_pl_min = t_range(1)*bin_size;
    t_pl_max = t_range(1)*bin_size+tmax;
    %xlims(i_f,:)
%     t_pl = (t_range(1):50000:t_range(2))*bin_size; % For common MIMO only
    t_pl = (t_range(1):t_range(2))*bin_size;
    
    %configure number of subplots
    in_pl=1;
    if ~sideLabel
        out_pl_real = 1:2;
        out_pl_pred = 3:4;
        N_pl = 4;
    else
        out_pl_real = 1:3;
        out_pl_pred = 4:6;
        N_pl = 6;
    end
    if plotInput
        N_pl = N_pl+3;
        in_pl = 1:3;
        out_pl_real = out_pl_real+3;
        out_pl_pred = out_pl_pred +3;
    elseif plotObservedOnly
        in_pl = 1:3;
        out_pl_real = 4:6;
        N_pl=6;
    end
    
    %plot event times
    if ~isempty(fSession)
        N_pl = N_pl+1;
        in_pl = in_pl+1;
        out_pl_real = out_pl_real+1;
        out_pl_pred = out_pl_pred +1;
        i_pl = sub2ind([Nf N_pl],i_f,1);
        %subplot(N_pl,1,1)
        if vert
            subplot(Nf,N_pl,i_pl)
        else
            subplot(N_pl,Nf,i_pl)
        end
        thisAx = gca;
            ax = [thisAx ax];
        if ~isempty(xlims), xlim(xlims(i_f,:)),end
        if ~labelAllEvents
            event_labels = {'Sample Present','Sample Response','Match Present','Match Response'};
        else
            event_labels = {'Sample Present','Sample Response','Match Present','Correct Response','Error Response','Focus','Reward'};
        end
        linsts = ':-:-:-:-';
        %     ccs = cbrewer('qual','Set1',4);
        %      ccs = ccs([1 1 2 2 3 3 4 4],:);
        ye = [1 1 0];
        wh = [1 1 1];
        bl = [0 0 1];
        re = [1 0 0];
        gr = [0 .7 0];
        cy = [0 1 1];
        mag = [1 0 1];
        ccs= [ye; ye; wh; wh; bl; bl; re; re];
        ccs= [ye; ye; wh; wh; ye; ye; wh; wh];
        %load(fSession,'DMS_Correct','DMS_Error','DMS_Match_Present','DMS_Match_Response','DMS_Sample_Response','DMS_Sample_Present');
        load(fSession);
        
        if exist('CORRECT_RESPONSE','var')
            if ~labelAllEvents
                if exist('ERROR_RESPONSE','var')
                    MATCH_RESPONSE = [ERROR_RESPONSE; CORRECT_RESPONSE];
                else
                    MATCH_RESPONSE = MATCH_RESPONSE;
                end
                if exist('MATCH_ON', 'var')
                    MATCH_ON = MATCH_ON;
                else
                    MATCH_ON = MATCH_RESPONSE - 2;
                end
                events = {'SAMPLE_ON','SAMPLE_RESPONSE','MATCH_ON','MATCH_RESPONSE'};
            else
                events = {'SAMPLE_ON','SAMPLE_RESPONSE','MATCH_ON','CORRECT_RESPONSE','ERROR_RESPONSE','FOCUS_ON','REWARD_ON'};
                ccs = [ye; ye; wh; gr; re; cy; mag];
                linsts = ':-:--::';
            end
        elseif exist('DMS_Sample_Present','var')
            if ~PAL
            events = {'DMS_Sample_Present','DMS_Sample_Response','DMS_Match_Present','DMS_Match_Response'};
            else
                events = {'DMS_Sample_Present','DMS_Sample_Response','DMS_Match_Present','DMS_Match_Response','PAL_Sample_Present','PAL_Sample_Response','PAL_Match_Present','PAL_Match_Response'};
            end
        elseif exist('Sample_On','var')
            events = {'Sample_On','Sample_Resp','Match_On','Match_Resp'};
        elseif exist('SamplePresent','var')
            events = {'SamplePresent','SampleResponse','MatchPresent','MatchResponse'};
        else
            warning('DNMS Events Not Loaded')
        end
        for i = 1:length(events)
            if exist(events{i},'var')  
            e = eval(events{i});
            cc = ccs(i,:);
            linst = linsts(i);
            h(i) = stem(e,ones(length(e),1),'color',cc,'marker','none','linewidth',0.5,'linestyle',linst);
            end
            hold on
        end
        %set(gca, 'XTickLabelMode', 'Manual')
        if ~vert  %unsrue why this is in horizontal mode only...
        set(gca, 'XTick', [])
        end
        set(gca, 'YTick', [])
        if vert
        xlabel('Time (s)','fontsize',label_font_size)
        view(90,90)  %hopefully rotates the plot!
        end
        if dispLegend
        legend(h,event_labels,'fontsize',14,'Location','east','Color','k','TextColor','w')
        legend(h,event_labels,'Location','east','Color','k','TextColor','w')
        %legend(h,event_labels,'Location','bestoutside','Color','k','TextColor','w')
        end
        xlim([t_pl_min t_pl_max])
        if ~isempty(xlims), xlim(xlims(i_f,:)),end
        title(titleLabel,'fontsize',title_font_size)
        if ~isempty(correctLabel)
            if correctLabel(i_f)==1
                title('Correct','fontsize',title_font_size,'color','g')
            else
                title('Error','fontsize',title_font_size,'color','r')
            end
        end
        set(gca,'Color','k')
    end
    
    
    
    unstable_inds = find(y_unstable);
    unstable_labels = yLabels(unstable_inds);
    [~,unsort_inds] = sort(sorted_ySp_inds);
    if ~outSort
        yLabels = yLabels(unsort_inds);
        y_real = y_real(:,unsort_inds);
        y_sim = y_sim(:,unsort_inds);
    else %sort the inputs by firing rate as well
        x_sps = sum(x_real);
        [~,x_sps_sort_inds] = sort(x_sps,'descend');
        x_real = x_real(:,[x_sps_sort_inds]);
    end
    
    %% Filter x_real spike timing
    for i=1:size(x_real,2)
        if isempty(reBinSize)
            x_filt(:,i) = gaussFilt(x_real(:,i),sigma)./bin_size;
        else
            tr = reBin(x_real(:,i),bin_size,reBinSize);
            x_filt(:,i) = gaussFilt(tr',sigma*bin_size/reBinSize)*reBinSize;  %%%this line needs to be fixed.
        end
    end
    %% cellTogether is mostly used operation
    if ~cellTogether
        %% Filter y_real and y_sim spike timing
        for i=1:nF
            if isempty(reBinSize)
                y_real_filt(:,i) = gaussFilt(y_real(:,i),sigma)./bin_size;
                y_sim_filt(:,i) = gaussFilt(y_sim(:,i),sigma)./bin_size;
            else
                tr_r = reBin(y_real(:,i),bin_size,reBinSize);
                tr_sim = reBin(y_sim(:,i),bin_size,reBinSize);
                y_real_filt(:,i) = gaussFilt(tr_r',sigma*bin_size/reBinSize)./reBinSize;
                y_sim_filt(:,i) = gaussFilt(tr_sim',sigma*bin_size/reBinSize)./reBinSize;
            end
        end
        if isempty(clims)
            clims = [min((y_real_filt(:))) max((y_real_filt(:)))*max_map_scale];
        end
        %% Plot Input if specified
        if plotInput || plotObservedOnly
            %subplot(N_pl,1,in_pl),
            i_pl = sub2ind([Nf N_pl],ones(size(in_pl))*i_f,in_pl);
            if new_figs
                figure()
            end
            subplot(N_pl,Nf,i_pl)
            thisAx = gca;
            ax = [thisAx ax];
            hh=imagesc(t_pl,1:nF,x_filt',clims);
            if delete_fig
                delete(hh)
            end
            set(gca, 'XTickLabelMode', 'Manual')
            set(gca, 'XTick', [])
            if disp_color_bar
                colorbar()
            end
            xlim([t_pl_min t_pl_max])
            if ~isempty(xlims), xlim(xlims(i_f,:)),end
            if i_f ==1
            ylabel('Input Cell','fontsize',label_font_size)
            end
            title(['Input Spatio-Temporal Pattern'],'fontsize',title_font_size)
        end
        
        %% plot Actual pattern
        %subplot(N_pl,1,out_pl_real),
        i_pl = [];
        for i_subpl = 1:length(out_pl_real)
            i_pl = [i_pl sub2ind([Nf N_pl],i_f,out_pl_real(i_subpl))];
        end
        if new_figs
                figure()
        end
        if vert
            subplot(Nf,N_pl,i_pl)
            hh=imagesc(1:nF,t_pl,y_real_filt,clims);
        else
            subplot(N_pl,Nf,i_pl)
            hh=imagesc(t_pl,1:nF,y_real_filt',clims);
        end
        thisAx = gca;
            ax = [thisAx ax];
        if delete_fig
            delete(hh)
        end
        if vert
            set(gca, 'YTickLabelMode', 'Manual')
            set(gca, 'YTick', [])
        else
            if  ~plotObservedOnly
                set(gca, 'XTickLabelMode', 'Manual')
                set(gca, 'XTick', [])
            end
            if labelCells
                set(gca, 'YTickLabelMode', 'Manual')
                set(gca, 'YTick', 1:size(y_real_filt,2))
                set(gca, 'YTickLabel', yLabels)
                
            end
        end
           
        if disp_color_bar
            colorbar()
        end
        if vert
            ylim([t_pl_min t_pl_max])
        else
            xlim([t_pl_min t_pl_max])
        end
        if ~isempty(xlims), xlim(xlims(i_f,:)),end
        if vert
            if i_f ==1, xlabel('Output Cell','fontsize',label_font_size), end
        else
        if i_f ==1, ylabel('Output Cell','fontsize',label_font_size), end
        end
        if sideLabel 
            if i_f ==1, ylabel('Observed Output','fontsize',title_font_size);
            else
                if vert
                    set(gca, 'XTick', []);
                else
                    set(gca, 'YTick', []);
                end
            end
        else
            if extraLabel
                titStr = ['Observed Output Spatio-Temporal Pattern, T=' num2str(round(t_pl_min)) '-' num2str(round(t_pl_max)) 's, sig=' num2str(sigma) ', Sorted=' num2str(outSort)];
                title(titStr,'fontsize',title_font_size)
            else
                titStr='Observed Output Spatio-Temporal Pattern';
                title(titStr,'fontsize',title_font_size)
            end
        end
        
        %% plot Predicted Pattern
        %subplot(N_pl,1,out_pl_pred),
        if ~plotObservedOnly
            i_pl = [];
            for i_subpl = 1:length(out_pl_pred)
                i_pl = [i_pl sub2ind([Nf N_pl],i_f,out_pl_pred(i_subpl))];
            end
            %i_pl = sub2ind([Nf N_pl],i_f,out_pl_pred);
            if new_figs
                figure()
            end
            if vert
                subplot(Nf,N_pl,i_pl)
                hh=imagesc(1:nF,t_pl,y_sim_filt,clims);
            else
                subplot(N_pl,Nf,i_pl)
                hh=imagesc(t_pl,1:nF,y_sim_filt',clims);
            end
            thisAx = gca;
            ax = [thisAx ax];
            if delete_fig
                delete(hh)
            end
            if disp_color_bar
                colorbar()
            end
            if sideLabel
                if i_f==1, ylabel('Predicted Output','fontsize',title_font_size)
                else
                    if vert
                        set(gca, 'XTick', []);
                    else
                        set(gca, 'YTick', []);
                    end
                end
            else
                title(['Predicted Output Spatio-Temporal Pattern'],'fontsize',title_font_size)
                if vert
                    xlabel('Output Cell','fontsize',label_font_size)
                    ylim([t_pl_min t_pl_max])
                else
                    ylabel('Output Cell','fontsize',label_font_size)
                    xlim([t_pl_min t_pl_max])
                    if labelCells
                        set(gca, 'YTickLabelMode', 'Manual')
                        set(gca, 'YTick', 1:size(y_real_filt,2))
                        set(gca, 'YTickLabel', yLabels)
                        
                    end
                end
            end
            
            if ~isempty(xlims), xlim(xlims(i_f,:)),end
            
        end
        if ~vert
            xlabel('Time (s)','fontsize',label_font_size)
        end
    else
        ax = [];
        for i=1:nF
            y_real_filt(:,i) = gaussFilt(y_real(:,i),sigma)./bin_size;
            y_sim_filt(:,i) = gaussFilt(y_sim(:,i),sigma)./bin_size;
            y_pl = [y_real_filt(:,i) y_sim_filt(:,i)];
            if new_figs
                figure()
            end
            thisAx = subplot(nF,1,i);
            ax = [ax thisAx];
            plot(t_pl,y_pl)
            box off
            if i~=nF
            set(gca, 'XTickLabelMode', 'Manual')
            set(gca, 'XTick', [])
            else
                legend('Recorded','Predicted')
            end
            %set(gca, 'YTickLabelMode', 'Manual')
            %set(gca, 'YTick', [])
            xlim([t_pl_min t_pl_max])
        end
        if new_figs
                figure()
            end
        subplot(nF,1,1), title(['Recorded Spiking on Top of Simulated, T=' num2str(round(t_pl_min)) '-' num2str(round(t_pl_max)) 's, sig=' num2str(sigma) ', Sorted=' num2str(outSort)])
    end
    
    %set(gcf,'position',[393   366   977   561])
    set(gcf,'position',[249         238        1055         663])
end
