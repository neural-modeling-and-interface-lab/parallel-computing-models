% Reg_fit controls the regularized fitting for one output neuron
%   
%   author: Dong Song
%           Brian Robinson, 2014-2016
%           Xiwei She, 2016-2022
classdef Reg_fit
    properties
        reg_fit_data  %This is what is created in Stage 1 when the regularization path is run
        xGroup
        yCell
        T
        T_trunc
        modelOptions
        regOptions
        errormessage
        fileLabel
        fitTime
        reMLETime
        CV_criteriaTime
        loading_error = 0;
        reMLE_CV_BIC_fit  %computed in stage 2:Here, reffitng with all criteria is run on the lambda value which corresponds to the best BIC from stage 1
        reMLE_CV_0th_fit %computed in stage 2: Here, reffitng with all criteria is run on a completely sparse model
        MLE_CV_fit %computed in stage 2: Here, refitting with all criteria is done on the full model
        CV_criteria %Here, all criteria is run
        fit_reMLE_CV   %triggers whether or not reMLE fitting is done
        fit_CV_criteria  %triggers whether or not CV criteria is done!
        refit_different_basis
        sessionIntOptions %options that stores how the modeled time bins should relate to epxerimental event timing
        sessionOb %Session object that stores trial event information
        maxIncludeIntervals
    end
    methods
        function obj = Reg_fit(xGroup,yCell,modelOptions,regOptions,fileLabel,varargin)
            [obj.T_trunc, obj.fit_reMLE_CV, obj.fit_CV_criteria, obj.sessionIntOptions, obj.sessionOb, obj.maxIncludeIntervals] = ...
                process_options(varargin, 'T_trunc',[],'fit_reMLE_CV',true,'fit_CV_criteria',true,'sessionIntOptions',[],'sessionOb',[],'maxIncludeIntervals',[]); %this gives us the option to run the fit on a truncated portion of the data for debugging
            if ~isempty(obj.T_trunc) && ~isempty(obj.sessionIntOptions)
                error('T_tunc OR sessionIntOptions can be specified, not both')
            end
            
            
            obj.xGroup=xGroup;
            obj.yCell=yCell;
            obj.modelOptions=modelOptions;
            obj.regOptions=regOptions;
            obj.fileLabel=fileLabel;
            
            %this line is tried initially to make sure intervals can be
            %retrieved without throwing an error
            [~, ~] = obj.getIncludeRemoveIntervals;
        end
        function [includeIntervals, removeIntervals] = getIncludeRemoveIntervals(obj)
            includeIntervals = [];
            removeIntervals = [];
            if ~isempty(obj.sessionIntOptions)
                if isempty(obj.sessionOb)
                    error('session IntOptions can only be specified if there is a corresponding Session Object')
                end
                binsize = obj.modelOptions.bin_size/1000; %note, all session binsizes are in s, while the model options bin_Sizes are in ms
                [includeIntervals, removeIntervals] = obj.sessionOb.getSessionIntervals(obj.sessionIntOptions,binsize);
                if ~isempty(obj.maxIncludeIntervals)
                    includeIntervals = includeIntervals(1:obj.maxIncludeIntervals,:);
                end
            end
            if ~isempty(obj.T_trunc) && ~isempty(obj.sessionIntOptions)
                error('T_tunc OR sessionIntOptions can be specified, not both')
            end
        end
        function [x,y]= getXY(obj,varargin)
            sess_t = process_options(varargin,'sess_t',0);  %this sess_t will return times and time bins corresponding to the original times recorded in the experiment, not in the time bins that correspond to the session data that is spliced togethr for fitting.
            
            xCells = obj.xGroup.cellArray;
            Nx = length(xCells);
            obj.T=0;
            if ~isempty(obj.T_trunc)
                if iscell(obj.T_trunc)   %this is the case when we want to concatenate many!!!
                    for i=1:length(obj.T_trunc)
                        t_min(i) = obj.T_trunc{i}(1)/1000*obj.modelOptions.bin_size;
                        t_max(i) = obj.T_trunc{i}(2)/1000*obj.modelOptions.bin_size;
                    end
                else
                    if length(obj.T_trunc)==1
                        t_min = 0;
                        t_max = obj.T_trunc/1000*obj.modelOptions.bin_size;
                    elseif length(obj.T_trunc)==2
                        t_min = obj.T_trunc(1)/1000*obj.modelOptions.bin_size;
                        t_max = obj.T_trunc(2)/1000*obj.modelOptions.bin_size;
                    end
                end
            end
            %parse spike times and find max duration!
            try
                for i=1:Nx
                    if isempty(obj.T_trunc)  || sess_t
                        x_sp{i}=xCells(i).getSpTimes('bin_size',obj.modelOptions.bin_size);
                    else
                        x_sp{i}=xCells(i).getSpTimesRangeSeconds(t_min,t_max,'bin_size',obj.modelOptions.bin_size);
                    end
%                     if ~isempty(obj.T_trunc)           %this is where we implement truncation of spike train
%                         if length(obj.T_trunc)==1
%                         x_sp{i}=x_sp{i}(x_sp{i}<=obj.T_trunc);
%                         elseif length(obj.T_trunc) ==2
%                         x_sp{i} = x_sp{i}(x_sp{i}>obj.T_trunc(1) && x_sp{i}<=obj.T_trunc(2));
%                         x_sp{i} = x_sp{i}-obj.T_trunc(1); 
%                         end
%                     end
                    this_T = ceil(max(x_sp{i}));
                    if this_T>obj.T
                        obj.T=this_T;
                    end
                end
             
                if isempty(obj.T_trunc)  || sess_t
                    y_sp = obj.yCell.getSpTimes('bin_size',obj.modelOptions.bin_size);
                else
                    y_sp=obj.yCell.getSpTimesRangeSeconds(t_min,t_max,'bin_size',obj.modelOptions.bin_size);
                end
%                 if ~isempty(obj.T_trunc)               %this is where we implement truncation of spike train
%                     y_sp=y_sp(y_sp<=obj.T_trunc);
%                 end
                this_T = ceil(max(y_sp));
                if this_T>obj.T
                    obj.T=this_T;
                end
                %create vector of binary spike times!
                x=zeros(obj.T,Nx);
                for i=1:Nx
                    % Exception found by Xiwei - 2017.7.25
                    if round(x_sp{i}(1)) < 1
                        x_sp{i}(1) = 1;
                    end
                    x(round(x_sp{i}),i)=1;
                end
                y=zeros(obj.T,1);
                y(round(y_sp))=1;
            catch err
                obj.loading_error=1;
                obj.errormessage=err;
            end
        end
        function fit = run_stage1(obj,folder)
            tstart = tic;
            %initialize and load spike times
            [x,y]= obj.getXY;
            
            %run fits
            if ~obj.loading_error
                try
                    %create Design Matrix
                    tic
                    try
                        disp(['Memory Ratio is ' num2str(memRatio) ', before design matrix'])
                    end
                    [includeIntervals, removeIntervals] = obj.getIncludeRemoveIntervals;
                    [V, g] = createDesign(x,y,obj.modelOptions.Bff,obj.modelOptions.Bfb,obj.modelOptions.Qff,obj.modelOptions.Qfb,obj.getOptff,obj.getOptfb,'includeIntervals', includeIntervals, 'removeIntervals',removeIntervals);
                    %                     V=sparse(V);
                    %                     y=sparse(y);
                    try
                        disp(['Memory Ratio is ' num2str(memRatio) ', after design matrix creation'])
                    end
                    model_options_str = ['Qff' num2str(obj.modelOptions.Qff) '_Qfb' num2str(obj.modelOptions.Qfb)];
                    disp('I just created the design matrix')
                    toc
                    
                    %do the fit
                    tic
                    y = trimSeries(y, includeIntervals,removeIntervals);
                    obj.reg_fit_data = groupLassoPath(V,y,g,obj.regOptions.link,obj.regOptions.max_it, obj.regOptions.n_lams,obj.regOptions.min_factor,obj.regOptions.coeff_tolerance,obj.regOptions.cv_K,obj.regOptions.cv_seed,obj.regOptions.mle_refitting);
                    disp('I just finished fitting')
                    toc
                catch err
                    obj.errormessage=err;
                end
            end
            obj.fitTime = toc(tstart);
            %save results
            fit=obj;
            %if nargin==2  %here we have to save things
            %   folder = varargin{1};
            yLab = obj.yCell.name(3:5);
            %                 fname = ['Reg_fit_y' yLab '_sess' obj.yCell.session '_' obj.fileLabel '_stagePath_d' num2str(now) '.mat'];
            fname = ['MIMOResult_' obj.yCell.session '_' obj.fileLabel '_stagePath_' yLab '.mat'];
            save([folder filesep fname],'fit')
            obj.saveCstim(folder,1);
            %end
        end
        function fit = run_stage2(obj,folder)
            tstart = tic;
            obj.reMLE_CV_0th_fit=obj.reMLE_CV_0th;
            obj.reMLE_CV_BIC_fit=obj.reMLE_CV_BIC;
            obj.MLE_CV_fit=obj.MLE_CV;
            obj.reMLETime = toc(tstart);
            fit = obj;
            %if nargin==2  %here we have to save things
                %folder = varargin{1};
                yLab = obj.yCell.name(3:5);
%                 fname = ['Reg_fit_y' yLab '_sess' obj.yCell.session '_' obj.fileLabel '_stageReMLE_d' num2str(now) '.mat'];
                fname = ['MIMOResult_' obj.yCell.session '_' obj.fileLabel '_stageReMLE_' yLab '.mat'];
                save([folder filesep fname],'fit')
                obj.saveCstim(folder,2);
            %end
        end
        function fit = run(obj,varargin)
            %do the first stage
            folder = varargin{1};
            obj=obj.run_stage1(folder);
            %fit = obj.run_stage1(folder);
            fit=obj;
            
            %% do the reMLE fitting
            disp('Done Stage 1! Run Stage 2 Now')
            if obj.fit_reMLE_CV
                obj= obj.run_stage2(folder);
                fit = obj;
            end
            %% do the CV criteria
            disp('Done Stage 2! Run Stage 3 Now')
            if obj.fit_CV_criteria
                tstart = tic;
                obj.CV_criteria = obj.calcCVCriteria;
                obj.CV_criteriaTime = toc(tstart);
                fit = obj;
                obj.save_stage3(folder)
%                 if nargin==2  %here we have to save things
%                     folder = varargin{1};
%                     fname = ['Reg_fit_y' obj.yCell.label '_sess' obj.yCell.session '_' obj.fileLabel '_stageCVCrit_d' num2str(now) '.mat'];
%                     save([folder filesep fname],'fit')
%                 end
            end
            
        end
        function save_stage3(obj,folder)
            fit = obj;
            yLab = obj.yCell.name(3:5);
            %             fname = ['Reg_fit_y' yLab '_sess' obj.yCell.session '_' obj.fileLabel '_stageCVCrit_d' num2str(now) '.mat'];
            fname = ['MIMOResult_' obj.yCell.session '_' obj.fileLabel '_stageCVCrit_' yLab '.mat'];
%             save([folder filesep fname],'fit', '-v7.3')
            save([folder filesep fname],'fit')
            obj.saveCstim(folder,3);
        end
        function [n_sig_inputs, best_lambda, min_BIC_lam_ind, lams, max_its, BIC, best_BIC_sig_groups, best_BIC_c,sig_input_cells,sig_groups_1,sig_groups_2] = getBICSparsity(obj)
            %disp(['REg fit data object is ' num2str(obj.reg_fit_data) ])
            f = obj.reg_fit_data;
            lams = f.lams;
            nL = length(lams);
            max_its = f.max_iterateds;
            for i=1:nL  %here we parse out the BIC calculated at the last iteration for each lambda value
                n_it=max_its(i);
                cr = f.criterias{i};
                BIC(i) = cr.BIC_count(n_it);
            end
            min_BIC_lam_ind = find(BIC==min(BIC));
            %             if isempty(min_BIC_lam_ind)  %in the case that this is a really bad fit and all of the BIC scores are Nan or inf, we just choose the first index
            %                 min_BIC_lam_ind = 1;
            %             end
            sig_groups = ~f.sparse_groups;
            n_sig_inputs_per_lambda = sum(sig_groups);
            n_sig_inputs = n_sig_inputs_per_lambda(min_BIC_lam_ind);
            best_lambda = lams(min_BIC_lam_ind);
            best_BIC_sig_groups = find(sig_groups(:,min_BIC_lam_ind));
            BIC_prime = diff(BIC);
            best_BIC_c = f.cs(:,min_BIC_lam_ind);
            if obj.modelOptions.Qff ~= 0 % Modified by Xiwei Dec.14, 2021 for Qff=0
                sig_input_cells = best_BIC_sig_groups/obj.modelOptions.Qff;  %the following three lines divide the groups by the model order to make this Qff invariant
            else
                sig_input_cells = [];
            end
            sig_input_cells = ceil(sig_input_cells);
            sig_input_cells = unique(sig_input_cells);
            if obj.modelOptions.Qff == 0 % Modified by Xiwei Dec.14, 2021 for Qff=0
                sig_groups_1 = [];
                sig_groups_2 = [];
            elseif obj.modelOptions.Qff == 1 
                sig_groups_1 = sig_input_cells;
                sig_groups_2 = [];
            else
                sig_groups_2 = best_BIC_sig_groups(mod(best_BIC_sig_groups,2)==0);   %all the even groups are the 2nd order!
                sig_groups_1 = best_BIC_sig_groups(~ismember(best_BIC_sig_groups,sig_groups_2));  %all the non-2nd order groups are the first order groups!
                sig_groups_2 = sig_groups_2/2; %convert the groups to input cell!
                sig_groups_1 = ceil(sig_groups_1/2); %convert the groups to input cell!
            end
        end
        function  BICSparsityStruct = getBICSparsityStruct(obj)
            BICSparsityStruct.N_possible_groups = size(obj.reg_fit_data.sparse_groups,1);
            [n_sig_inputs, best_lambda, min_BIC_lam_ind, lams, max_its, BIC, best_BIC_sig_groups, best_BIC_c,sig_input_cells,sig_groups_1,sig_groups_2] = obj.getBICSparsity;
            BICSparsityStruct.n_sig_inputs = n_sig_inputs;
            BICSparsityStruct.best_lambda=best_lambda;
            BICSparsityStruct.min_BIC_lam_ind=min_BIC_lam_ind;
            BICSparsityStruct.lams=lams;
            BICSparsityStruct.max_its=max_its;
            BICSparsityStruct.BIC=BIC;
            BICSparsityStruct.best_BIC_sig_groups=best_BIC_sig_groups;
            BICSparsityStruct.best_BIC_c= best_BIC_c;
            BICSparsityStruct.sig_input_cells=sig_input_cells;
            BICSparsityStruct.sig_groups_1=sig_groups_1;
            BICSparsityStruct.sig_groups_2=sig_groups_2;
        end
        function CVSparsityStruct = getCVSparsityStruct(obj)
            CVSparsityStruct.N_possible_groups = size(obj.reg_fit_data.sparse_groups,1);
            [n_sig_inputs, best_lambda, min_CV_lam_ind, lams, ~, av_Loss, best_CV_sig_groups, best_CV_c, meanKSS_at_best_CV_Loss,sig_input_cells,sig_groups_1,sig_groups_2] = obj.getCVSparsity;
            CVSparsityStruct.n_sig_inputs=n_sig_inputs;
            CVSparsityStruct.best_lambda=best_lambda;
            CVSparsityStruct.min_CV_lam_ind=min_CV_lam_ind;
            CVSparsityStruct.lams=lams;
            CVSparsityStruct.av_Loss=av_Loss;
            CVSparsityStruct.best_CV_sig_groups=best_CV_sig_groups;
            CVSparsityStruct.best_CV_c=best_CV_c;
            CVSparsityStruct.meanKSS_at_best_CV_Loss=meanKSS_at_best_CV_Loss;
            CVSparsityStruct.sig_input_cells=sig_input_cells;
            CVSparsityStruct.sig_groups_1=sig_groups_1;
            CVSparsityStruct.sig_groups_2=sig_groups_2;
        end
        function R = getBestCVFit(obj)
            f = obj.reg_fit_data;
            if ~isempty(f)
                lams = f.lams;
                nL = length(lams);
                av_Loss = nan(nL,1);
                Cr = obj.CV_criteria;
                for i=1:nL
                    av_Loss(i) = mean(Cr(i).outSampleLoss);
                end
                min_CV_lam_ind = find(av_Loss==min(av_Loss));
                min_CV_lam_ind=min_CV_lam_ind(end); 
            end
            R = obj.CV_criteria(min_CV_lam_ind);
        end
        function R = getCVSparsityOffset5(obj)
            R = obj.getCVSparsityOffset(.05);
        end
        function R = getCVSparsityOffset10(obj)
            R = obj.getCVSparsityOffset(.10);
        end
        function R = getCVSparsityOffset15(obj)
            R = obj.getCVSparsityOffset(.15);
        end
        function R = getCVSparsityOffset20(obj)
            R = obj.getCVSparsityOffset(.20);
        end
        function R = getCVSparsityOffset(obj,offset)
            f = obj.reg_fit_data;
            if ~isempty(f)
                lams = f.lams;
                nL = length(lams);
                Cr = obj.CV_criteria;
                for i=1:nL
                    av_Loss(i) = mean(Cr(i).outSampleLoss);
                end
                stdLoss = av_Loss - min(av_Loss);
                max_no_inf = max(stdLoss(~isinf(stdLoss)));
                if max_no_inf~=0
                stdLoss = stdLoss/max_no_inf;
                end
                ind_min_loss = find(stdLoss==0);
                ind_min_loss = ind_min_loss(end);  %find the index which is the sparsest!
                for i=ind_min_loss:-1:1           %move towards sparser solutions!
                    if stdLoss(i)<offset          %if the Loss here is still better than the offset threshold, choose this sparsity level!
                        ind_offset = i; 
                    else                          %however, if the loss goes above the threshold stop looking!
                        break
                    end
                end
                
                R.fit = obj.CV_criteria(ind_offset);
                R.lam_ind = ind_offset;
                
                sig_groups = ~f.sparse_groups;
                n_sig_inputs_per_lambda = sum(sig_groups);
                best_CV_sig_groups = find(sig_groups(:,ind_offset));
                if obj.modelOptions.Qff ~= 0 % Modified by Xiwei Dec.14, 2021 for Qff=0
                    sig_input_cells = best_CV_sig_groups/obj.modelOptions.Qff;  %the following three lines divide the groups by the model order to make this Qff invariant
                else
                    sig_input_cells = [];
                end
                sig_input_cells = ceil(sig_input_cells);
                R.sig_input_cells = unique(sig_input_cells);
                if obj.modelOptions.Qff == 0 % Modified by Xiwei Dec.14, 2021 for Qff=0
                    R.sig_groups_1 = [];
                    R.sig_groups_2 = [];
                elseif obj.modelOptions.Qff == 1
                    R.sig_groups_1 = sig_input_cells;
                    R.sig_groups_2 = [];
                else
                    sig_groups_2 = best_CV_sig_groups(mod(best_CV_sig_groups,2)==0);   %all the even groups are the 2nd order!
                    sig_groups_1 = best_CV_sig_groups(~ismember(best_CV_sig_groups,sig_groups_2));  %all the non-2nd order groups are the first order groups!
                    R.sig_groups_2 = sig_groups_2/2; %convert the groups to input cell!
                    R.sig_groups_1 = ceil(sig_groups_1/2); %convert the groups to input cell!
                end
                R.n_sig_inputs = n_sig_inputs_per_lambda(ind_offset);
            end
        end
        function [n_sig_inputs, best_lambda, min_CV_lam_ind, lams, empty, av_Loss, best_CV_sig_groups, best_CV_c, meanKSS_at_best_CV_Loss,sig_input_cells,sig_groups_1,sig_groups_2] = getCVSparsity(obj)
            empty = [];  %this is a dummy variable to make the output of this functino equivalent to that of getBICSparsity
            f = obj.reg_fit_data;
            if ~isempty(f)
                lams = f.lams;
                nL = length(lams);
                Cr = obj.CV_criteria;
                for i=1:nL
                    av_Loss(i) = mean(Cr(i).outSampleLoss);
                    %meanKSS(i) = mean(Cr(i).KS.KS_score);
                end
                min_CV_lam_ind = find(av_Loss==min(av_Loss));
                min_CV_lam_ind=min_CV_lam_ind(end);  %this takes the smallest lambda if there are multiple!!!
                sig_groups = ~f.sparse_groups;
                n_sig_inputs_per_lambda = sum(sig_groups);
                n_sig_inputs = n_sig_inputs_per_lambda(min_CV_lam_ind);
                best_lambda = lams(min_CV_lam_ind);
                best_CV_sig_groups = find(sig_groups(:,min_CV_lam_ind));
                best_CV_c = f.cs(:,min_CV_lam_ind);
                meanKSS_at_best_CV_Loss = mean(Cr(min_CV_lam_ind).outSampleLoss);
%                 x = obj.getXY;
%                 nIn = size(x,2);
                if obj.modelOptions.Qff ~= 0 % Modified by Xiwei Dec.14, 2021 for Qff=0
                    sig_input_cells = best_CV_sig_groups/obj.modelOptions.Qff;  %the following three lines divide the groups by the model order to make this Qff invariant
                else
                    sig_input_cells = [];
                end
                sig_input_cells = ceil(sig_input_cells);
                sig_input_cells = unique(sig_input_cells);
                if obj.modelOptions.Qff == 0 % Modified by Xiwei Dec.14, 2021 for Qff=0
                    sig_groups_1 = [];
                    sig_groups_2 = [];
                elseif obj.modelOptions.Qff ==1
                    sig_groups_1 = sig_input_cells;
                    sig_groups_2 = [];
                else
                    sig_groups_2 = best_CV_sig_groups(mod(best_CV_sig_groups,2)==0);   %all the even groups are the 2nd order!
                    sig_groups_1 = best_CV_sig_groups(~ismember(best_CV_sig_groups,sig_groups_2));  %all the non-2nd order groups are the first order groups!
                    sig_groups_2 = sig_groups_2/2; %convert the groups to input cell!
                    sig_groups_1 = ceil(sig_groups_1/2); %convert the groups to input cell!
                end
                %sig_input_cells=sig_input_cells(sig_input_cells<=nIn);   %this makes it so that it doesn't count the last input as significant!
            else
                n_sig_inputs=nan;
                best_lambda=nan;
                min_CV_lam_ind=nan;
                lams=nan;
                av_Loss=nan;
                best_CV_sig_groups=nan;
                best_CV_c=nan;
                meanKSS_at_best_CV_Loss=nan;
                sig_groups_2 = nan;
                sig_groups_1 = nan;
            end
        end
        function [zeroth_mu, zeroth_SE, full_mu, full_SE, BIC_mu, BIC_SE, CV_mu, CV_SE] = getKS_scores(obj)
            zF = obj.reMLE_CV_0th_fit;
            fF = obj.MLE_CV_fit;
            bF = obj.reMLE_CV_BIC_fit;
            if ~isempty(zF)
                [~, ~, min_CV_lam_ind, ~, ~, ~, ~, ~, ~] = obj.getCVSparsity;
                cvF = obj.CV_criteria(min_CV_lam_ind);
                zeroth_mu = mean(zF.KS.KS_score);
                zeroth_SE = std(zF.KS.KS_score)/sqrt(length(zF.KS.KS_score));
                full_mu = mean(fF.KS.KS_score);
                full_SE = std(fF.KS.KS_score)/sqrt(length(fF.KS.KS_score));
                BIC_mu = mean(bF.KS.KS_score);
                BIC_SE = std(bF.KS.KS_score)/sqrt(length(bF.KS.KS_score));
                CV_mu = mean(cvF.KS.KS_score);
                CV_SE = std(cvF.KS.KS_score)/sqrt(length(cvF.KS.KS_score));
            else
                zeroth_mu=nan;
                zeroth_SE=nan;
                full_mu=nan;
                full_SE=nan;
                BIC_mu=nan;
                BIC_SE=nan;
                CV_mu=nan;
                CV_SE=nan;
            end
        end
        function KSStruct = getKSStruct(obj)
            [zeroth_mu, zeroth_SE, full_mu, full_SE, BIC_mu, BIC_SE, CV_mu, CV_SE] = obj.getKS_scores;
            mu.zeroth = zeroth_mu;
            mu.full =full_mu;
            mu.BIC = BIC_mu;
            mu.CV = CV_mu;
            SE.zeroth = zeroth_SE;
            SE.full = full_SE;
            SE.BIC = BIC_SE;
            SE.CV = CV_SE;
            KSStruct.mu=mu;
            KSStruct.SE = SE;
        end
        function NYFit = getNYFit(obj)
            [~,y]=obj.getXY;
            [includeIntervals, removeIntervals] = obj.getIncludeRemoveIntervals;
            y = trimSeries(y, includeIntervals,removeIntervals);
            NYFit = sum(y);
        end
        function CVLossStruct = getCVLossStruct(obj)
            [zeroth_mu, zeroth_SE, full_mu, full_SE, BIC_mu, BIC_SE, CV_mu, CV_SE] = obj.getLosses;
            mu.zeroth = zeroth_mu;
            mu.full =full_mu;
            mu.BIC = BIC_mu;
            mu.CV = CV_mu;
            SE.zeroth = zeroth_SE;
            SE.full = full_SE;
            SE.BIC = BIC_SE;
            SE.CV = CV_SE;
            CVLossStruct.mu = mu;
            CVLossStruct.SE = SE;
        end
        function cv_design = is_cv_design(obj)
        %checks to see if cross validation should be done directly on the
        %design matrix
        
        %Note: the way cross validation is done can be changed by adding
        %regOptions.cv_design = 1, in order to maintain backwards
        %compatibility, we first check if the field exists
        cv_design = 0;
        if isfield(obj.regOptions,'cv_design')
            if obj.regOptions.cv_design==1
                cv_design = 1;
            end
        end
        
        
        end
        function [zeroth_mu, zeroth_SE, full_mu, full_SE, BIC_mu, BIC_SE, CV_mu, CV_SE] = getLosses(obj)
            zF = obj.reMLE_CV_0th_fit;
            fF = obj.MLE_CV_fit;
            bF = obj.reMLE_CV_BIC_fit;
            if ~isempty(zF)
                [~, ~, min_CV_lam_ind, ~, ~, ~, ~, ~, ~] = obj.getCVSparsity;
                cvF = obj.CV_criteria(min_CV_lam_ind);
                zeroth_mu = mean(zF.outSampleLoss);
                zeroth_SE = std(zF.outSampleLoss)/sqrt(length(zF.outSampleLoss));
                full_mu = mean(fF.outSampleLoss);
                full_SE = std(fF.outSampleLoss)/sqrt(length(fF.outSampleLoss));
                BIC_mu = mean(bF.outSampleLoss);
                BIC_SE = std(bF.outSampleLoss)/sqrt(length(bF.outSampleLoss));
                CV_mu = mean(cvF.outSampleLoss);
                CV_SE = std(cvF.outSampleLoss)/sqrt(length(cvF.outSampleLoss));
            else
                zeroth_mu=nan;
                zeroth_SE=nan;
                full_mu=nan;
                full_SE=nan;
                BIC_mu=nan;
                BIC_SE=nan;
                CV_mu=nan;
                CV_SE=nan;
            end
        end
        %% 
        %           function [zeroth_mu, zeroth_SE, full_mu, full_SE, BIC_mu, BIC_SE, CV_mu, CV_SE] = getKS_scores(obj)
        %               zF = obj.reMLE_CV_0th_fit;
        %               fF = obj.MLE_CV_fit;
        %               bF = obj.reMLE_CV_BIC_fit;
        %               [~, ~, min_CV_lam_ind, ~, ~, ~, ~, ~, ~] = obj.getCVSparsity;
        %               cvF = obj.CV_criteria(min_CV_lam_ind);
        %               zeroth_mu = mean(zF.KS.KS_score);
        %               zeroth_SE = std(zF.KS.KS_score)/sqrt(length(zF.KS.KS_score));
        %               full_mu = mean(fF.KS.KS_score);
        %               full_SE = std(fF.KS.KS_score)/sqrt(length(fF.KS.KS_score));
        %               BIC_mu = mean(bF.KS.KS_score);
        %               BIC_SE = std(bF.KS.KS_score)/sqrt(length(bF.KS.KS_score));
        %               CV_mu = mean(cvF.KS.KS_score);
        %               CV_SE = std(cvF.KS.KS_score)/sqrt(length(cvF.KS.KS_score));
        %           end
        %
        function [zKS, fKS, bKS, cvKS] = getKS_plot(obj)
            %this reutrns the first KS plot out of 5
            zF = obj.reMLE_CV_0th_fit;
            fF = obj.MLE_CV_fit;
            bF = obj.reMLE_CV_BIC_fit;
            [~, ~, min_CV_lam_ind, ~, ~, ~, ~, ~, ~] = obj.getCVSparsity;
            cvF = obj.CV_criteria(min_CV_lam_ind);
            zKS = getFirstKS(zF.KS);
            bKS = getFirstKS(bF.KS);
            fKS = getFirstKS(fF.KS);
            cvKS = getFirstKS(cvF.KS);
        end     
        %% Kernel reconstruction functions
        %re-MLE kernels!!
        function K = get_K_CV_reMLE(obj,varargin)
            [max_basis_M, makeResponse,r2Offset,k1Resp] = process_options(varargin,'max_basis_M',[],'makeResponse',0,'r2Offset',[4 20 100],'k1Resp',1);
            [~, ~, min_CV_lam_ind, ~, ~, ~, sparse_g, ~, ~] = obj.getCVSparsity;
            c = get_av_inSample_c(obj.CV_criteria(min_CV_lam_ind));
            c = c_sparse2full(c,sparse_g,obj.modelOptions.Bff,obj.modelOptions.Bfb,obj.modelOptions.Qff,obj.modelOptions.Qfb,obj.getOptff,obj.getOptfb,length(obj.xGroup.cellArray));
            K=K_Reconstruct(c,obj.modelOptions.Bff,obj.modelOptions.Bfb,obj.modelOptions.Qff,obj.modelOptions.Qfb,obj.getOptff,obj.getOptfb,'max_basis_M',max_basis_M,'makeResponse',makeResponse,'r2Offset',r2Offset,'k1Resp',k1Resp);
        end
        function K = get_K_BIC_reMLE(obj,varargin)
            max_basis_M = process_options(varargin,'max_basis_M',[]);
            c = get_av_inSample_c(obj.reMLE_CV_BIC_fit);
            [~,~,~,~,~,~,sparse_g,~]=obj.getBICSparsity;
            
            c = c_sparse2full(c,sparse_g,obj.modelOptions.Bff,obj.modelOptions.Bfb,obj.modelOptions.Qff,obj.modelOptions.Qfb,obj.getOptff,obj.getOptfb,length(obj.xGroup.cellArray));
            K=K_Reconstruct(c,obj.modelOptions.Bff,obj.modelOptions.Bfb,obj.modelOptions.Qff,obj.modelOptions.Qfb,obj.getOptff,obj.getOptfb,'max_basis_M',max_basis_M);
        end
        %gLasso kernels at BIC and CV indexes
        function K = get_K_CV_gLasso(obj,varargin)
            max_basis_M = process_options(varargin,'max_basis_M',[]);
            [~,~,~,~,~,~,~,c]=obj.getCVSparsity;
            K=K_Reconstruct(c,obj.modelOptions.Bff,obj.modelOptions.Bfb,obj.modelOptions.Qff,obj.modelOptions.Qfb,obj.getOptff,obj.getOptfb,'max_basis_M',max_basis_M);
        end
        function K = get_K_BIC_gLasso(obj,varargin)
            max_basis_M = process_options(varargin,'max_basis_M',[]);
            [~,~,~,~,~,~,~,c]=obj.getBICSparsity;
            K=K_Reconstruct(c,obj.modelOptions.Bff,obj.modelOptions.Bfb,obj.modelOptions.Qff,obj.modelOptions.Qfb,obj.getOptff,obj.getOptfb,'max_basis_M',max_basis_M);
        end
        %full kernel
        function K = get_K_full_MLE(obj,varargin)
            max_basis_M = process_options(varargin,'max_basis_M',[]);
            c = get_av_inSample_c(obj.MLE_CV_fit);
            K=K_Reconstruct(c,obj.modelOptions.Bff,obj.modelOptions.Bfb,obj.modelOptions.Qff,obj.modelOptions.Qfb,obj.getOptff,obj.getOptfb,'max_basis_M',max_basis_M);
        end    
        %% Getting kernel coefficient functions
        function [c, c_sp] = get_c_CV_reMLE(obj,varargin) %this returns stage3 result that had lowest CV Loss
            [fold, tilde] = process_options(varargin, 'fold',[],'tilde', 0);
            [~, ~, min_CV_lam_ind, ~, ~, ~, sparse_g, ~, ~] = obj.getCVSparsity;
            c = get_av_inSample_c(obj.CV_criteria(min_CV_lam_ind),'fold',fold,'tilde', tilde);
            c_sp = c;
            c = c_sparse2full(c,sparse_g,obj.modelOptions.Bff,obj.modelOptions.Bfb,obj.modelOptions.Qff,obj.modelOptions.Qfb,obj.getOptff,obj.getOptfb,length(obj.xGroup.cellArray));
            %K=K_Reconstruct(c,obj.modelOptions.Bfb,obj.modelOptions.Bff,obj.modelOptions.Qfb,obj.modelOptions.Qff);
        end
        function [c, c_sp] = get_c_CV_reMLE_per_Fold(obj) %this returns stage3 result that had lowest CV Loss
            [~, ~, min_CV_lam_ind, ~, ~, ~, sparse_g, ~, ~] = obj.getCVSparsity;
            c = get_av_inSample_c(obj.CV_criteria(min_CV_lam_ind));
            c_sp = c;
            c = c_sparse2full(c,sparse_g,obj.modelOptions.Bff,obj.modelOptions.Bfb,obj.modelOptions.Qff,obj.modelOptions.Qfb,obj.getOptff,obj.getOptfb,length(obj.xGroup.cellArray));
            %K=K_Reconstruct(c,obj.modelOptions.Bfb,obj.modelOptions.Bff,obj.modelOptions.Qfb,obj.modelOptions.Qff);
        end
        function [c1, c2, c1_fb, c2_fb, c0_OR_sig] = get_c_byIn_CV_reMLE(obj,varargin)
            [fold, tilde] = process_options(varargin, 'fold',[],'tilde', 0);
            %c1 is L1 x N, c2 is L2 x N
            c = obj.get_c_CV_reMLE('fold',fold,'tilde', tilde);
            L =  obj.modelOptions.Bff.getN;
            N = size(obj.xGroup.cellArray,2);
            L1 = getNgrCoeffs(1,L,[]);
            if obj.modelOptions.Qff ==2
                L2 = getNgrCoeffs(2,L,obj.getOptff);
            else
                L2 = 0;
            end
            [c1, c2] = get_c_by_input(c.k,L1,L2,N);  %c1 is L1 x N, c2 is L2 x N
            [c1_fb, c2_fb] = get_c_by_input(c.h,L1,L2,1);  %c1 is L1 x N, c2 is L2 x N
            if tilde
                c0_OR_sig = c.c_0;
            else
                c0_OR_sig = c.sig;
            end
        end
        function Q3fbK = getK_fbQ3(obj)
            c_all = obj.get_c_CV_reMLE;
            ch = c_all.h;
            Bfb = obj.modelOptions.Bfb;
            L = Bfb.getN;
            N12 = L+L*(L+1)/2;
            c = ch(N12+1:end);
            if isfield(obj.modelOptions, 'Q3fbBand')
                Q3 = obj.modelOptions.Q3fbBand;
            else
                Opt = obj.getOptfb;
                Q3=Opt.Q3;
            end
            B = Bfb.calcB;
            n_sig = 0;
            sig_ns = zeros(Q3.Nbands,1);
            for i=1:Q3.Nbands
                c_b = c((i-1)*L+1:i*L);
                K_temp= B*c_b;
                K(:,i) = K_temp;
                if norm(c_b)~=0
                    n_sig = n_sig+1;
                    sig_ns(i)=1;
                end
            end
            Q3fbK.K3=K;
            Q3fbK.n_sig = n_sig;
            Q3fbK.sig_ns = sig_ns;
            Q3fbK.K1 = B*ch(1:L);
        end
        function [c, c_sp] = get_c_BIC_reMLE(obj)  %this looks at the stage 2 best BIC results and returns the relevant coefficients
            c = get_av_inSample_c(obj.reMLE_CV_BIC_fit);
            [~,~,~,~,~,~,sparse_g,~]=obj.getBICSparsity;
            c_sp = c;
            
            c = c_sparse2full(c,sparse_g,obj.modelOptions.Bff,obj.modelOptions.Bfb,obj.modelOptions.Qff,obj.modelOptions.Qfb,obj.getOptff,obj.getOptfb,length(obj.xGroup.cellArray));
            %K=K_Reconstruct(c,obj.modelOptions.Bfb,obj.modelOptions.Bff,obj.modelOptions.Qfb,obj.modelOptions.Qff);
        end
        %gLasso kernels at BIC and CV indexes
        function c = get_c_CV_gLasso(obj) %this uses the criteria calculated in stage 3 and returns the path result in stage 1
            [~,~,~,~,~,~,~,c]=obj.getCVSparsity;
            %K=K_Reconstruct(c,obj.modelOptions.Bfb,obj.modelOptions.Bff,obj.modelOptions.Qfb,obj.modelOptions.Qff);
        end
        function c = get_c_BIC_gLasso(obj) %this uses the criteria calculate in stage 1 and returns the path result in stage 1
            [~,~,~,~,~,~,~,c]=obj.getBICSparsity;
            %K=K_Reconstruct(c,obj.modelOptions.Bfb,obj.modelOptions.Bff,obj.modelOptions.Qfb,obj.modelOptions.Qff);
        end
        %full kernel
        function c = get_c_full_MLE(obj)
            c = get_av_inSample_c(obj.MLE_CV_fit);
            %K=K_Reconstruct(c,obj.modelOptions.Bfb,obj.modelOptions.Bff,obj.modelOptions.Qfb,obj.modelOptions.Qff);
        end  
        %% Simulate output
        function [y, unstable, a, u, Pf, n] = simulateOutput(obj,t_range,seed,varargin)
            [unstable_a_thresh, full, unstable_timeout, stage2]=process_options(varargin,'unstable_a_thresh',1,'full',0,'unstable_timeout',[],'stage2',0);
            [x,~]= obj.getXY('sess_t',1);  %here, we are using the times from the original experimental recording!
            if t_range(2)<size(x,1)
                x = x(t_range(1):t_range(2),:);
            else
                x = x(t_range(1):end,:);
            end
            if stage2
                c=obj.get_c_BIC_reMLE;
            elseif ~full
                c = obj.get_c_CV_reMLE;
            else
                c=obj.get_c_full_MLE;
            end
            [y, unstable, a, u, Pf, n] = MISO_sim(x,c,obj.modelOptions.Bff,obj.modelOptions.Bfb,obj.modelOptions.Qff,obj.modelOptions.Qfb,obj.getOptff,obj.getOptfb,seed,'unstable_a_thresh',unstable_a_thresh,'unstable_timeout',unstable_timeout);
        end
        function [y, unstable, a, u, Pf, n] = simulateOutput_different_basis(obj,t_range,seed,varargin)
            unstable_a_thresh=process_options(varargin,'unstable_a_thresh',1);
            [x,~]= obj.getXY('sess_t',1);  %here, we are using the times from the original experimental recording!;
            if t_range(2)<size(x,1)
                x = x(t_range(1):t_range(2),:);
            else
                x = x(t_range(1):end,:);
            end
            c = obj.refit_different_basis.c;
            if ~isfield(obj.refit_different_basis,'Optff')
                obj.refit_different_basis.Optff=[];
            end
            if ~isfield(obj.refit_different_basis,'Optfb')
                obj.refit_different_basis.Optfb=[];
            end
            [y, unstable, a, u, Pf, n] = MISO_sim(x,c,obj.refit_different_basis.Bff, obj.refit_different_basis.Bfb,obj.modelOptions.Qff,obj.modelOptions.Qfb,obj.refit_different_basis.Optff,obj.refit_different_basis.Optfb,seed,'unstable_a_thresh',unstable_a_thresh);
        end
        function K = get_K_different_basis(obj)
            c = obj.refit_different_basis.c;
            K=K_Reconstruct(c,obj.refit_different_basis.Bff,obj.refit_different_basis.Bfb,obj.refit_different_basis.Qff,obj.refit_different_basis.Qfb,obj.getOptff,obj.getOptfb);
        end
        function refit_best_CV_different_basis(obj,Bff,Bfb,Optff,Optfb,folder,label,varargin)
            T_trunc = process_options(varargin,'T_trunc',[]);  %this is an additional T_trunc used in this function to test it out with shorter amounts of time!
            %refit with best CV sparsity lambda value
            [~, ~, ~, ~, ~, ~, sparse_g, ~, ~] = obj.getCVSparsity;
            [x,y]= obj.getXY;
            [includeIntervals, removeIntervals] = obj.getIncludeRemoveIntervals;
            if ~isempty(T_trunc)
                x = x(1:T_trunc,:);
                y = y(1:T_trunc,:);
                if ~isempty(includeIntervals) || ~isempty(removeIntervals) 
                    error('T_trunc has to be empty if includeIntervals and removeIntervals are specified for data selection')
                end
            end
            
            refit = MISO(x,y,Bff,Bfb,obj.modelOptions.Qff,obj.modelOptions.Qfb,Optff,Optfb,'g_include',sparse_g,...
                'includeIntervals',includeIntervals,'removeIntervals',removeIntervals);
            c_sp = refit.c;
            
            c=c_sparse2full(c_sp,sparse_g,Bff,Bfb,obj.modelOptions.Qff,obj.modelOptions.Qfb,Optff,Optfb,length(obj.xGroup.cellArray));
            %save everything
            obj.refit_different_basis.c=c;
            obj.refit_different_basis.refit = refit;
            obj.refit_different_basis.Bff=Bff;
            obj.refit_different_basis.Bfb=Bfb;
            obj.refit_different_basis.Optff = Optff;
            obj.refit_different_basis.Optfb = Optfb;
            fit = obj;
            fname = ['Reg_fit_y' obj.yCell.label '_sess' obj.yCell.session '_' obj.fileLabel '_stageDiffBasis' label '_d' num2str(now) '.mat'];
            save([folder filesep fname],'fit')
        end
        %% %%%%%%%%%%%%%%%This is an old versino of the function that used cross validation cirteria calculated with a serpate path for each fold.  The new updated one is CV based on a single path from all the data!
        %           function [n_sig_inputs, best_lambda, min_CV_lam_ind, lams, empty, av_Loss, best_CV_sig_groups, best_CV_c] = getCVSparsity(obj,Reg_Or_MLE)
        %               empty = [];  %this is a dummy variable to make the output of this functino equivalent to that of getBICSparsity
        %               f = obj.reg_fit_data;
        %               lams = f.lams;
        %               nL = length(lams);
        %               if strcmp(Reg_Or_MLE,'Reg')
        %                   av_Loss = f.CV.av_Loss_CV_Regs;
        %                   CV_cs = f.CV.c_av_CV_Regs;
        %               elseif strcmp(Reg_Or_MLE,'MLE')
        %                   av_Loss = f.CV.av_Loss_CV_MLEs;
        %                   CV_cs = f.CV.c_av_CV_MLEs;
        %               else
        %                   disp('Options are MLE or Reg')
        %               end
        %               min_CV_lam_ind = find(av_Loss==min(av_Loss));
        %               min_CV_lam_ind=min_CV_lam_ind(end);
        %               n_sparse_inputs_per_lambda=zeros(nL,1);
        %               n_sig_inputs_per_lambda=zeros(nL,1);
        %               sparse_gs=cell(nL,1);
        %               sig_gs=cell(nL,1);
        %               for i_lam = 1:nL
        %                   [n_sparse_inputs_per_lambda(i_lam), n_sig_inputs_per_lambda(i_lam),sparse_gs{i_lam},sig_gs{i_lam}] = sparseCount(CV_cs(:,i_lam),obj.modelOptions.Bff,obj.modelOptions.Bfb,obj.modelOptions.Qff,obj.modelOptions.Qfb,length(obj.xGroup.cellArray));
        %               end
        %               n_sig_inputs = n_sig_inputs_per_lambda(min_CV_lam_ind);
        %               best_lambda = lams(min_CV_lam_ind);
        %               best_CV_sig_groups = sig_gs{min_CV_lam_ind};
        %               best_CV_c = CV_cs(:,min_CV_lam_ind);
        %           end
        
        
        
        function [n_sparse_per_lambda, mono_decrease] = getSparsePerLambda(obj)
            f = obj.reg_fit_data;
            n_sparse_per_lambda = sum(f.sparse_groups)';
            n_sparse_prime = diff(n_sparse_per_lambda);
            increasing = n_sparse_prime > 0;
            mono_decrease=sum(increasing)==0;  %this returns whether the sparsities is monotonically decreasing (it should be)
        end
        
        function K = getK_gLasso_BIC_best(obj)
            [n_sig_inputs, best_lambda, min_BIC_lam_ind, lams, max_its, BIC] = getBICSparsity(obj);
            c_temp = obj.reg_fit_data.cs(:,min_BIC_lam_ind);
            K=K_Reconstruct(c_temp,obj.modelOptions.Bff,obj.modelOptions.Bfb,obj.modelOptions.Qff,obj.modelOptions.Qfb,obj.getOptff,obj.getOptfb);
            
        end
        
        function f = reMLE_CV_0th(obj,varargin)
            T_trunc = process_options(varargin,'T_trunc',[]);
            [x,y]= obj.getXY;
            [includeIntervals, removeIntervals] = obj.getIncludeRemoveIntervals;
            if ~isempty(T_trunc)
                x=x(1:T_trunc,:);
                y = y(1:T_trunc,:);
                if ~isempty(includeIntervals) || ~isempty(removeIntervals) 
                    error('T_trunc has to be empty if includeIntervals and removeIntervals are specified for data selection')
                end
            end
            %added code to allow new regOptions.cv_design
            if obj.is_cv_design
                f = MISO_cv_design(x,y,obj.modelOptions.Bff,obj.modelOptions.Bfb,obj.modelOptions.Qff,obj.modelOptions.Qfb,obj.getOptff,obj.getOptfb,'g_include',[],'LCD',obj.modelOptions.MLE_LCD,'LCD_max_it',obj.modelOptions.LCD_max_it,'K',obj.regOptions.K,...
                    'includeIntervals',includeIntervals,'removeIntervals',removeIntervals);
            else
                f = MISO_cv(x,y,obj.modelOptions.Bff,obj.modelOptions.Bfb,obj.modelOptions.Qff,obj.modelOptions.Qfb,obj.getOptff,obj.getOptfb,'g_include',[],'LCD',obj.modelOptions.MLE_LCD,'LCD_max_it',obj.modelOptions.LCD_max_it,'K',obj.regOptions.K,...
                    'includeIntervals',includeIntervals,'removeIntervals',removeIntervals);
            end
            obj.reMLE_CV_0th_fit = f;
        end
        function f = reMLE_CV_BIC(obj,varargin)
            T_trunc = process_options(varargin,'T_trunc',[]);
            [~, ~, ~, ~, ~, ~, best_BIC_sig_groups] = obj.getBICSparsity;
            [x,y]= obj.getXY;
            [includeIntervals, removeIntervals] = obj.getIncludeRemoveIntervals;
            if ~isempty(T_trunc)
                x=x(1:T_trunc,:);
                y = y(1:T_trunc,:);
                if ~isempty(includeIntervals) || ~isempty(removeIntervals) 
                    error('T_trunc has to be empty if includeIntervals and removeIntervals are specified for data selection')
                end
            end
            %added code to allow new regOptions.cv_design
            if obj.is_cv_design
                f = MISO_cv_design(x,y,obj.modelOptions.Bff,obj.modelOptions.Bfb,obj.modelOptions.Qff,...
                obj.modelOptions.Qfb,obj.getOptff,obj.getOptfb,'g_include',best_BIC_sig_groups,...
                'LCD',obj.modelOptions.MLE_LCD,'LCD_max_it',obj.modelOptions.LCD_max_it,'K',obj.regOptions.K,...
                    'includeIntervals',includeIntervals,'removeIntervals',removeIntervals);
            else
                f = MISO_cv(x,y,obj.modelOptions.Bff,obj.modelOptions.Bfb,obj.modelOptions.Qff,...
                obj.modelOptions.Qfb,obj.getOptff,obj.getOptfb,'g_include',best_BIC_sig_groups,...
                'LCD',obj.modelOptions.MLE_LCD,'LCD_max_it',obj.modelOptions.LCD_max_it,'K',obj.regOptions.K,...
                    'includeIntervals',includeIntervals,'removeIntervals',removeIntervals);
            end
            obj.reMLE_CV_BIC_fit = f;
        end
        function f = MLE_CV(obj,varargin)
            T_trunc = process_options(varargin,'T_trunc',[]);
            [x,y]= obj.getXY;
            [includeIntervals, removeIntervals] = obj.getIncludeRemoveIntervals;
            if ~isempty(T_trunc)
                x=x(1:T_trunc,:);
                y = y(1:T_trunc,:);
                if ~isempty(includeIntervals) || ~isempty(removeIntervals) 
                    error('T_trunc has to be empty if includeIntervals and removeIntervals are specified for data selection')
                end
            end
            %added code to allow new regOptions.cv_design
            if obj.is_cv_design
                 f = MISO_cv_design(x,y,obj.modelOptions.Bff,obj.modelOptions.Bfb,obj.modelOptions.Qff,...
                obj.modelOptions.Qfb,obj.getOptff,obj.getOptfb,'LCD',obj.modelOptions.MLE_LCD,...
                'LCD_max_it',obj.modelOptions.LCD_max_it,'K',obj.regOptions.K,...
                    'includeIntervals',includeIntervals,'removeIntervals',removeIntervals);
            else
                f = MISO_cv(x,y,obj.modelOptions.Bff,obj.modelOptions.Bfb,obj.modelOptions.Qff,...
                    obj.modelOptions.Qfb,obj.getOptff,obj.getOptfb,'LCD',obj.modelOptions.MLE_LCD,...
                    'LCD_max_it',obj.modelOptions.LCD_max_it,'K',obj.regOptions.K,...
                    'includeIntervals',includeIntervals,'removeIntervals',removeIntervals);
            end
            obj.MLE_CV_fit = f;
        end
        
        function f = calcCVCriteria(obj,varargin)
            T_trunc = process_options(varargin,'T_trunc',[]);
%             [x,y]= obj.getXY;
%             if ~isempty(T_trunc)
%                 x=x(1:T_trunc,:);
%                 y = y(1:T_trunc,:);
%             end
            spGs = obj.reg_fit_data.sparse_groups;
            if isempty(spGs)
                for i = 1:size(spGs, 2)
                    f(i) = obj.calcCVCriteria_individual_group([],'T_trunc',T_trunc);
                end
            else
                spG_Last = nan(size(spGs,1),1);
                for i=1:size(spGs,2)
                    if isequal(spGs(:,i),spG_Last)
                        f(i) = f(i-1);
                    else
                        sigG = find(~spGs(:,i));
                        %spG = find(spGs(:,i));
                        f(i) = obj.calcCVCriteria_individual_group(sigG,'T_trunc',T_trunc);
                        %f(i) = MISO_cv(x,y,obj.modelOptions.Bff,obj.modelOptions.Bfb,obj.modelOptions.Qff,obj.modelOptions.Qfb,obj.getOptff,obj.getOptfb,'g_include',sigG,'LCD',obj.modelOptions.MLE_LCD,'LCD_max_it',obj.modelOptions.LCD_max_it);
                    end
                    spG_Last=spGs(:,i);
                end
            end
        end
        
        function f = getUniqueCV(obj)
            spGs = obj.reg_fit_data.sparse_groups;
            spG_Last = nan(size(spGs,1),1);
            unique_inds = [];
            for i=1:size(spGs,2)
                if ~isequal(spGs(:,i),spG_Last)
                    unique_inds = [unique_inds i];
                end
                spG_Last=spGs(:,i);
            end
            f = obj.CV_criteria(unique_inds);
        end
        
        function N = getNUniqueCVInds(obj)
            spGs = obj.reg_fit_data.sparse_groups;
            spG_Last = nan(size(spGs,1),1);
            unique_inds = [];
            for i=1:size(spGs,2)
                if ~isequal(spGs(:,i),spG_Last)
                    unique_inds = [unique_inds i];
                end
                spG_Last=spGs(:,i);
            end
            N = length(unique_inds);
        end
        
        function W = getAnyWarning(obj)
            fs = obj.getUniqueCV;
            W.Bscls = nan(length(fs),1);
            W.illCs = nan(length(fs),1);
            W.itRs = nan(length(fs),1);
            for i=1:length(fs)
                fits=[fs(i).inSampleFits.stats];                
                Bscl = fits.bad_scaling_warned;
                illC = [fits.ill_conditioned_warned];
                itR = [fits.iter_reached];
                W.Bscls(i) = ~isempty(find(Bscl, 1));
                W.illCs(i) = ~isempty(find(illC, 1));
                W.itRs(i) = ~isempty(find(itR, 1));
            end
            W.anyBadScale = ~isempty(find(W.Bscls,1));
            W.anyIllConditioned = ~isempty(find(W.illCs,1));
            W.anyIterationReached = ~isempty(find(W.itRs,1));
            W.anyWarned = ~isempty(find([W.Bscls W.illCs W.itRs ], 1));
        end
        
        function [f_i, i_y, i_lam] = calcCVCriteria_individual_group(obj,sigG,varargin)
            [T_trunc,i_lam,i_y] = process_options(varargin,'T_trunc',[],'i_lam',[],'i_y',[]);
            %the i_lam and i_y values are used to keep track of the results
            %of this function when called by other batch processesl
            [x,y]= obj.getXY;
            [includeIntervals, removeIntervals] = obj.getIncludeRemoveIntervals;
            if ~isempty(T_trunc)
                x=x(1:T_trunc,:);
                y = y(1:T_trunc,:);
                if ~isempty(includeIntervals) || ~isempty(removeIntervals) 
                    error('T_trunc has to be empty if includeIntervals and removeIntervals are specified for data selection')
                end
            end
            %added code to allow new regOptions.cv_design
            if obj.is_cv_design
                f_i = MISO_cv_design(x,y,obj.modelOptions.Bff,obj.modelOptions.Bfb,obj.modelOptions.Qff,...
                    obj.modelOptions.Qfb,obj.getOptff,obj.getOptfb,'g_include',sigG,'LCD',...
                    obj.modelOptions.MLE_LCD,'LCD_max_it',obj.modelOptions.LCD_max_it,'K',obj.regOptions.K,...
                    'includeIntervals',includeIntervals,'removeIntervals',removeIntervals);
            else
                f_i = MISO_cv(x,y,obj.modelOptions.Bff,obj.modelOptions.Bfb,obj.modelOptions.Qff,...
                    obj.modelOptions.Qfb,obj.getOptff,obj.getOptfb,'g_include',sigG,'LCD',...
                    obj.modelOptions.MLE_LCD,'LCD_max_it',obj.modelOptions.LCD_max_it,'K',obj.regOptions.K,...
                    'includeIntervals',includeIntervals,'removeIntervals',removeIntervals);
            end
            
        end
        
        function pl_info = getSparsePathPlotInfo(obj)
            %get sparsity information
            n_sparse_per_lambda = obj.getSparsePerLambda;
            pl_info.sparse_N = to01(n_sparse_per_lambda);
            %get BIC information, this should always be available!
            %if ~isempty(obj.reMLE_CV_BIC_fit)
                [~, ~, min_BIC_lam_ind, lams, ~, BIC,~, ~] = obj.getBICSparsity;
                pl_info.BIC_N=to01(BIC);
                pl_info.bestBIC_ind=min_BIC_lam_ind;
                pl_info.bestBIC = pl_info.BIC_N(pl_info.bestBIC_ind);
                pl_info.lams = lams;
            %else
%                 pl_info.BIC_N=[];
%                 pl_info.bestBIC_ind=[];
%                 pl_info.bestBIC = [];
%                 pl_info.lams = obj.reg_fit_data.lams;
%             end
            %get CV information
            if ~isempty(obj.CV_criteria)
                [~, ~, min_CV_lam_ind, ~, ~, av_Loss, ~, ~, ~] = obj.getCVSparsity;
                pl_info.av_Loss_N=to01(av_Loss);
                pl_info.bestCVLoss_ind = min_CV_lam_ind;
                pl_info.bestCVLoss = pl_info.av_Loss_N(min_CV_lam_ind);
            else
                pl_info.av_Loss_N=[];
                pl_info.bestCVLoss_ind = [];
                pl_info.bestCVLoss = [];
            end
        end
        function Opt = getOptff(obj)
            if isfield(obj.modelOptions,'Optff')
                Opt = obj.modelOptions.Optff;
            else
                Opt = [];
            end
        end
        function Opt = getOptfb(obj)
            if isfield(obj.modelOptions,'Optfb')
                Opt = obj.modelOptions.Optfb;
            else
                Opt = [];
            end
        end
        function SparseChangeVsTime = getStage3FitTimeVsSparsityChange(obj)
            spGChange = diff(obj.reg_fit_data.sparse_groups,1,2);  %this has how each group changes for each lambda iteration
            totalSpGChange = sum(spGChange);                      %this is the total amount of changes in each lambda iteration
            spGChange_inds = [1 (1+find(totalSpGChange))];  %the first evaluation
            nPathGChange = length(spGChange_inds);         %this is the total amont of iterations that have a change
            coeff_per_lambda = sum(obj.reg_fit_data.cs~=0);
            coeff_per_glm_lambda=coeff_per_lambda(spGChange_inds);
            total_coeffs_estimated = sum(coeff_per_glm_lambda);
            
            for i = 1:length(spGChange_inds)
                i_lam = (spGChange_inds(i));
                CV_lam = obj.CV_criteria(i_lam);
                CV_stats = [CV_lam.inSampleFits.stats];
                total_lam_iters(i) = sum([CV_stats.max_iter]);
                av_lam_glm_iters(i) = mean([CV_stats.max_iter]);
            end
            SparseChangeVsTime.sum_mle_iters = sum(total_lam_iters);
            SparseChangeVsTime.sum_mle_iters_times_coeffs = total_lam_iters*coeff_per_lambda(spGChange_inds)';
            SparseChangeVsTime.av_lam_glm_iters=av_lam_glm_iters;         
            SparseChangeVsTime.coeff_per_glm_lambda=coeff_per_glm_lambda;
            [~, ~, ~, ~, ~, av_Loss, ~, ~, ~] = obj.getCVSparsity;
            av_Loss_N=to01(av_Loss);
            
            SparseChangeVsTime.av_Loss_N_eval = av_Loss_N(spGChange_inds);
            
            SparseChangeVsTime.total_coeffs_estimated=total_coeffs_estimated;
            SparseChangeVsTime.nPathGChange = nPathGChange;
            SparseChangeVsTime.CV_criteriaTime=obj.CV_criteriaTime;
        end
        function InpInf = getInpInf(obj,varargin)
            % loops thtough all selected input / yff and fits the model
            % without that selected input / yff (keeping the all of the
            % reamining selected inputs / yff).  Calculates the loss.
            T_trunc = process_options(varargin,'T_trunc',[]);
            c = obj.get_c_CV_reMLE;
            c = c2c_temp(c);
            [x,y]= obj.getXY;
            [includeIntervals, removeIntervals] = obj.getIncludeRemoveIntervals;
            if ~isempty(T_trunc)
                x = x(1:T_trunc,:);
                y = y(1:T_trunc);
                if ~isempty(includeIntervals) || ~isempty(removeIntervals) 
                    error('T_trunc has to be empty if includeIntervals and removeIntervals are specified for data selection')
                end
            end
            y_trim = trimSeries(y, includeIntervals,removeIntervals); %trimmed y by intervals needs to be used for loss calculation (full y needs to be used for design matrix creation)
            V = createDesign(x,y,obj.modelOptions.Bff,obj.modelOptions.Bfb,obj.modelOptions.Qff,obj.modelOptions.Qfb,obj.getOptff,obj.getOptfb,...
                'includeIntervals',includeIntervals,'removeIntervals',removeIntervals);
            L_full = ProbitLoss(c, V, y_trim)./length(y_trim);
            L=ones(size(x,2)+1,1)*L_full;
            x_calc = obj.getCVSparsityStruct.sig_input_cells;
            for i=x_calc'
                tstart = tic;
                x_temp = x;
                if i<=size(x,2)
                    x_temp(:,i) = zeros(length(x),1);
                    y_temp = y;
                else
                    y_temp = zeros(length(x),1);
                end
                V = createDesign(x_temp,y_temp,obj.modelOptions.Bff,obj.modelOptions.Bfb,obj.modelOptions.Qff,obj.modelOptions.Qfb,obj.getOptff,obj.getOptfb,...
                'includeIntervals',includeIntervals,'removeIntervals',removeIntervals);
                tend = toc(tstart);
                disp(['Done Design in ' num2str(i) ' out of ' num2str(size(x,2)) ' in ' num2str(tend) 's'] )
                tstart = tic;
                L(i) = ProbitLoss(c, V, y_trim)./length(y_trim);
                tend = toc(tstart);
                disp(['Done Loss in ' num2str(i) ' out of ' num2str(size(x,2)) ' in ' num2str(tend) 's'] )
            end
            inp_DI = L-L_full;
            InpInf.inp_DI = inp_DI;
            InpInf.L=L;
            InpInf.L_full = L_full;
        end
        function T = getStage2FitTime(obj)
            T.BIC =sum([obj.reMLE_CV_BIC_fit.inSampleFits.Telapsed]);
            T.zero =sum([obj.reMLE_CV_0th_fit.inSampleFits.Telapsed]);
            T.full =sum([obj.MLE_CV_fit.inSampleFits.Telapsed]);
            T.total = T.BIC+T.zero+T.full;
        end
        function T = getStage3FitTime(obj)
            for i=1:length( obj.CV_criteria)
                Tlam(i) = sum([obj.CV_criteria(i).inSampleFits.Telapsed]);
            end
            T.lam = Tlam;
            T.unique = unique(T.lam);
            T.total = sum(T.unique);
        end
%         function Cstim = getCStim(obj,sigma)
%             [~, ~, min_CV_lam_ind, ~, ~, ~, sparse_g, ~, ~] = obj.getCVSparsity;
%             c_tilde = get_av_inSample_c_tilde(obj.CV_criteria(min_CV_lam_ind));
%             c_tilde_sp = c_tilde;
%             c_tilde = c_sparse2full(c_tilde_sp,sparse_g,obj.modelOptions.Bff,obj.modelOptions.Bfb,obj.modelOptions.Qff,obj.modelOptions.Qfb,obj.getOptff,obj.getOptfb,length(obj.xGroup.cellArray));
%             Cstim = c_tilde2C_stim(c_tilde, sigma);
%         end
%         function saveCstim(obj,outputPath,sigma,prefix)
%             C=obj.getCStim(sigma);
%             NeuronIndex = parseHuman6CellName(obj.yCell.name);
%             f = [prefix num2str(NeuronIndex) '.mat'];
%             fsave = fullfile(outputPath,f);
%             modelInfo.regOptions=obj.regOptions;
%             modelInfo.modelOptions=obj.modelOptions;
%             modelInfo.T_trunc = obj.T_trunc;
%             modelInfo.xGroup = obj.xGroup;
%             modelInfo.yCell = obj.yCell;
%             save(fsave,'C','sigma','modelInfo');
%         end
        function Cstim = getCStimStage3(obj,sigma)  %this is stage 3 where we get
            [~, ~, min_CV_lam_ind, ~, ~, ~, sparse_g, ~, ~] = obj.getCVSparsity;  %this returns the sparse_g that corresponds with the best CV fit
            c_tilde_sp = get_av_inSample_c_tilde(obj.CV_criteria(min_CV_lam_ind)); %this returns the average sparse c_tilde for the best DV fit
            c_tilde = c_sparse2full(c_tilde_sp,sparse_g,obj.modelOptions.Bff,obj.modelOptions.Bfb,obj.modelOptions.Qff,obj.modelOptions.Qfb,obj.getOptff,obj.getOptfb,length(obj.xGroup.cellArray));
            Cstim = c_tilde2C_stim(c_tilde, sigma);
        end
        function Cstim = getCStimStage2(obj,sigma)  %this is stage 2 where we get the cross-validated fit at the best CV
            [~,~,~,~,~,~,sparse_g,~]=obj.getBICSparsity; %this returns the sparsity at the best BIC fit
            c_tilde_sp = get_av_inSample_c_tilde(obj.reMLE_CV_BIC_fit); %this returns the sparse tilde at the best BIC fit
            c_tilde = c_sparse2full(c_tilde_sp,sparse_g,obj.modelOptions.Bff,obj.modelOptions.Bfb,obj.modelOptions.Qff,obj.modelOptions.Qfb,obj.getOptff,obj.getOptfb,length(obj.xGroup.cellArray));
            Cstim = c_tilde2C_stim(c_tilde, sigma);
        end
        function Cstim = getCStimStage1(obj,sigma)  %this is stage 1 where we get the group lass fit with the Best BIC
            BIC_struct = obj.getBICSparsityStruct;
            c_tilde_unsegmented = BIC_struct.best_BIC_c;
            c_tilde.c_0=c_tilde_unsegmented(1);
            Ncoeffs=obj.getNcoeffs;
            c_tilde.k = c_tilde_unsegmented(2:Ncoeffs.ff+1);
            c_tilde.h = c_tilde_unsegmented(Ncoeffs.ff+2:end);
            Cstim = c_tilde2C_stim(c_tilde, sigma);
        end
        function Cstim = getCStimFull(obj,sigma) %this gets the coefficiets from the full model
            c_tilde = get_av_inSample_c_tilde(obj.MLE_CV_fit);
            Cstim = c_tilde2C_stim(c_tilde, sigma);
        end
        function Ncoeffs=getNcoeffs(obj)
            Nin=length(obj.xGroup.cellArray);
            Ncoeffs.ff = order2numc(obj.modelOptions.Bff.getN, obj.modelOptions.Qff)*Nin;
            Ncoeffs.fb = order2numc(obj.modelOptions.Bfb.getN, obj.modelOptions.Qfb);
        end
        function saveCstim(obj,folder,stage,varargin)  %the default saving options are arbitraary but correspond to legacy formats
            [prefix,sigma]=process_options(varargin,'prefix','Order2Out','sigma',10);
            %specify and create saving folder
            %outpath = fullfile(folder,['CoeffStage' num2str(stage)]);
            outpath = folder;
            if ~exist(outpath,'dir')
                mkdir(outpath)
            end
            switch stage
                case 1
                    C=obj.getCStimStage1(sigma);
                case 2
                    C=obj.getCStimStage2(sigma);
                case 3
                    C=obj.getCStimStage3(sigma);
                case 'full'
                    C=obj.getCStimFull(sigma);
            end
            switch stage
                case 1
                    prefix = ['Stage1' prefix];
                case 2
                    prefix = ['Stage2' prefix];
                case 3
                    prefix = ['Stage3' prefix];
                case 'full'
                    prefix = ['StageFull' prefix];
            end
            if ~isempty(obj.yCell.BRIndex)  %this is the case when the yCell is created from a NEV file
                f = [prefix num2str(obj.yCell.BRIndex) '.mat'];
            else
                f = [prefix obj.yCell.name '.mat'];
            end
            fsave = fullfile(outpath,f);
            modelInfo.regOptions=obj.regOptions;
            modelInfo.modelOptions=obj.modelOptions;
            modelInfo.T_trunc = obj.T_trunc;
            modelInfo.xGroup = obj.xGroup;
            modelInfo.yCell = obj.yCell;
            finishedTime=datestr(now);
            save(fsave,'C','sigma','modelInfo','stage','finishedTime');
        end
        %% Deprecated functions
    end
end