classdef pMDM
    % This runs the Memory Decoding Model with parallel computing strategies
    % on USC CARC Cluster
    
    % Author: Xiwei She, Brian Robinson
    
    properties
        DataFolder % patient folder
        Category % category for regression
        oF % output file
        
        % fitting options
        num_bagging %number of bagging splits
        resolutions % number of b splines resolutions
        d % number of b spline knots
        L % memory length on either side of the event in ms
        par %this indicates if fitting is done in a parfor loop or not
        
        % intermediate variables
        c %this is the category classification used in estimation
        SAMPLE_RESPONSE % Timestamp for matching neural signals with behaviorals
        R_exp %this stores the results of the non-control data
        tFit %elapsed fitting time
        
    end
    
    methods
        
        function obj = pMDM(DataFolder, Category, varargin)
            [obj.num_bagging, obj.resolutions, obj.d,obj.L, obj.par] = process_options(varargin,...
                'num_bagging',20,'resolutions',50:150,'d',3,'L',2000,'par',0);
            
            % Read information from recorded files
            iF1 = strcat('PATH_OF_BEHAVIORAL_DATA.mat');
            load(iF1, 'SAMPLE_RESPONSE', 'CATEGORY_INFO');
            c = CATEGORY_INFO(:, Category);
            
            %store object properties
            obj.c = c;
            obj.SAMPLE_RESPONSE = SAMPLE_RESPONSE;
            obj.DataFolder = DataFolder;
            obj.Category = Category;
            
            %specify output file and save memory decoding setup
            obj.oF = strcat('PATH_OF_OUTPUT_DATAs.mat');
            MDfit = obj;
            save(obj.oF,'MDfit')
        end
        
        function SpikeTensor = getSpikeTensor(obj)
            % Preprocessing
            iF2 = strcat('PATH_OF_NEURAL_DATA.mat');
            load(iF2, 'X', 'Y'); % X: CA3, Y:CA1
            
            SpikeTensor = Train2Tensor([X Y], obj.SAMPLE_RESPONSE, obj.L);
            
        end
        function thisSplit = runSplit(obj, SpikeTensor, varargin)
            % In each split
            
            % The resoluiton loop
            for rIndex = 1:length(obj.resolutions)
                thisResolution = obj.resolutions(rIndex);
                
                % MDM Estimation
                [modelOutputs] = MemoryDecodingModel(SpikeTensor, thisResolution, obj.d);
                
                % Save results
                thisSplit(rIndex).modelOutputs = modelOutputs; % Model outputs
                thisSplit(rIndex).tElapsed = tElapsed; % Running time
            end
        end
        function R = runAllSplits(obj, varargin)
            SpikeTensor=obj.getSpikeTensor;
            
            % The split loop
            parfor split = 1:obj.num_bagging
                R_sepSplit{split} = runSplit(SpikeTensor);
            end
            
            
            % Put R back in the original format
            for split = 1:obj.num_bagging
                for res = 1:length(obj.resolutions)
                    R(split,res) = R_sepSplit{split}(res);
                end
            end
            
        end
        function MDfit = run(obj,varargin)
            if obj.par
                poolOb = obj.setupPool;
            end
            tStart = tic;
            obj.R_exp = obj.runAllSplits;
            
            obj.tFit = toc(tStart);
            MDfit = obj;
            save(obj.oF,'MDfit','-v7.3')
            if obj.par
                poolOb.delete;
            end
        end
        function poolOb = setupPool(obj)
            if ~isunix
                poolOb = parpool;
            else  %here, the default is to put a seperate node per trial, with the number of workers equal to the number of trials
                nWorkers = obj.num_bagging;
                
                clusProf = get_SLURM_cluster('/home/rcf-proj/USER_ACCOUNT/matlab_storage','/usr/usc/matlab/R2020a/SlurmIntegrationScripts','--time=50:00:00 --partition USERNAME');
                
                poolOb = parpool(clusProf,nWorkers);
            end
        end
    end
    
end

