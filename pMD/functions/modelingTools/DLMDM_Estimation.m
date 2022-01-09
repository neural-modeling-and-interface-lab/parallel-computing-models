classdef DLMDM_Estimation
    % This runs the double-layer memory decoding for a given subject, 
    % category, and fitting options
    
    properties
        oF % output file
        
        % fitting options
        num_split %number of bagging splits
        m_all % array with number of b splines
        d % order of b spline knots
        L % memory length (decoding window)
        par % this indicates if fitting is done in parallel or not
        
        % intermediate variables
        target %this is the category classification used in estimation
        R_exp %this stores the results of the non-control data
        tFit % elapsed fitting time
        R_first % FirstLayer MD Results
        R_second % SecondLayer MD Results
        
        % Nest CV variables
        TrainingSet_SpikeTensor
        TrainingSet_target
        TestingSet_SpikeTensor
        TestingSet_target
    end
    
    methods
        function obj = DLMDM_Estimation(currentFold, runCase, varargin)
            [obj.num_split, obj.m_all, obj.d,obj.L, obj.par] = process_options(varargin,...
                'num_split',20, 'm_all',50:150, 'd',3, 'L',2000, 'par',0);
            warning off;
            % Load Data
            iF1 = strcat('exampleData\exampleInput\NestedCVDLMDM_', runCase, '_fold_', mat2str(currentFold),'.mat');
            load(iF1, 'target', 'TrainingSet_SpikeTensor', 'TrainingSet_target', 'TestingSet_SpikeTensor', 'TestingSet_target');
            
            %store object properties
            obj.target = target;
            obj.TrainingSet_SpikeTensor = TrainingSet_SpikeTensor;
            obj.TrainingSet_target = TrainingSet_target;
            obj.TestingSet_SpikeTensor = TestingSet_SpikeTensor;
            obj.TestingSet_target = TestingSet_target;
            
            %specify output file and save memory decoding setup
            obj.oF = strcat('exampleData\exampleResult\DLMDMEstimation_', runCase,'_fold_',mat2str(currentFold),'.mat');
            DLMDfit = obj;
            save(obj.oF,'DLMDfit')
        end

        function [firstR, secondR] = runSplit(obj, ti, varargin)
            % Training input and output
            SpikeTensor = obj.TrainingSet_SpikeTensor;
            thisTrial_c = obj.TrainingSet_target;
            
            % Testing input and output
            testingTensor = obj.TestingSet_SpikeTensor;
            testingTarget = obj.TestingSet_target;
            
            % get all parameters from object used in the fit of each trial
            firstR = struct;
            observes_secondLayer = zeros(length(thisTrial_c), length(obj.m_all));
            observes_secondLayer_testing = zeros(length(testingTarget), length(obj.m_all));
            
            % First Layer MD training
            disp(['============== First-layer MDM start on node ', mat2str(ti),' =============='])
            for mi = 1:length(obj.m_all)
                m = obj.m_all(mi);
                P = SpikeTensor2BSplineFeatureMatrix(SpikeTensor, m, obj.d);
                
                % Model Fit
                [FL_inside_ConfusionMatrix, FL_inside_predictions, FL_inside_Coefficients, FL_inside_FitInfo, FL_inside_CrossValSet_0, FL_inside_CrossValSet_1, FL_inside_probabilities, FL_inside_Deviance] = ModelFit(P, thisTrial_c);
                
                % Save first layer intermediate results
                firstR(mi).Resolution = mi; % current b-spline resolution
                firstR(mi).FL_inside_ConfusionMatrix = FL_inside_ConfusionMatrix; % confusion matrix
                firstR(mi).FL_inside_MCC = mcc(FL_inside_ConfusionMatrix); % MCC
                firstR(mi).FL_inside_predictions = FL_inside_predictions; % prediction
                firstR(mi).FL_inside_Coefficients = FL_inside_Coefficients; % Fiting Results of each fold - Coefficients
                firstR(mi).FL_inside_FitInfo = FL_inside_FitInfo; % Fiting Results of each fold - Fitting Information
                firstR(mi).FL_inside_CrossValSet_0 = FL_inside_CrossValSet_0; % Patition setting for label 0
                firstR(mi).FL_inside_CrossValSet_1 = FL_inside_CrossValSet_1; % Patition setting for label 1
                firstR(mi).FL_inside_probabilities = FL_inside_probabilities; % predicted probability
                firstR(mi).FL_inside_Deviance = FL_inside_Deviance;
                firstR(mi).FL_inside_Target = thisTrial_c;
                
                % Pick the global lambda
                globalMinDeviance = min(FL_inside_Deviance);
                globalIndex = find(FL_inside_Deviance == globalMinDeviance);
                if length(globalIndex) > 1
                    globalIndex = globalIndex(1);
                end
                
                % First Layer Nested testing
                P_test = SpikeTensor2BSplineFeatureMatrix(testingTensor, m, obj.d);
                
                numFold = length(FL_inside_Coefficients);
                FL_outside_probabilities = zeros(numFold, length(testingTarget));
                FL_outside_predictions = zeros(numFold, length(testingTarget));
                FL_outside_MCC = zeros(numFold, 1);
                FL_outside_GlobalCoefficients = zeros(numFold, size(FL_inside_Coefficients{1}, 1));
                FL_outside_GlobalC0 = zeros(numFold, 1);
                for tempI = 1:numFold
                    
                    globalC0 = FL_inside_FitInfo{tempI}.Intercept(globalIndex);
                    FL_outside_GlobalC0(tempI) = globalC0;
                    FL_outside_GlobalCoefficients(tempI, :) = FL_inside_Coefficients{tempI}(:, globalIndex);
                    c_i = P_test * FL_inside_Coefficients{tempI}(:, globalIndex) + globalC0;
                    c_p = 1 ./ (1 + exp(-c_i));
                    FL_outside_probabilities(tempI, :) = c_p;
                    FL_outside_predictions(tempI, :) = double(c_p>=0.5);
                    FL_outside_CM = confusionmat(testingTarget, double(c_p>=0.5));
                    if (size(FL_outside_CM,1)==1&&size(FL_outside_CM,2)==1)
                        FL_outside_CM = [FL_outside_CM(1,1) 0;0 0];
                    end
                    FL_outside_MCC(tempI, 1) = mcc(FL_outside_CM);
                    
                end
                
                % Average the predicted probabilities of all folds
                FL_outside_probabilities_foldAveraged = mean(FL_outside_probabilities, 1);
                FL_outside_CM_foldAveraged = confusionmat(testingTarget, double(FL_outside_probabilities_foldAveraged>=0.5));
                if (size(FL_outside_CM_foldAveraged,1)==1&&size(FL_outside_CM_foldAveraged,2)==1)
                    FL_outside_CM_foldAveraged = [FL_outside_CM_foldAveraged(1,1) 0;0 0];
                end
                FL_outside_MCC_foldAveraged = mcc(FL_outside_CM_foldAveraged);
                
                % Save FL testing results
                firstR(mi).FL_outside_GlobalCoefficients = FL_outside_GlobalCoefficients;
                firstR(mi).FL_outside_GlobalC0 = FL_outside_GlobalC0;
                firstR(mi).FL_outside_probabilities = FL_outside_probabilities;
                firstR(mi).FL_outside_predictions = FL_outside_predictions;
                firstR(mi).FL_outside_MCCs = FL_outside_MCC;
                firstR(mi).FL_outside_probabilities_foldAveraged = FL_outside_probabilities_foldAveraged;
                firstR(mi).FL_outside_CM_foldAveraged = FL_outside_CM_foldAveraged;
                firstR(mi).FL_outside_MCC_foldAveraged = FL_outside_MCC_foldAveraged;
                firstR(mi).FL_outside_Target = testingTarget;

                fprintf('Done split: %d; base learner: %d \n', ti, mi);
                % Save the probability for the second layer MD
                observes_secondLayer(:, mi) = FL_inside_probabilities;
                observes_secondLayer_testing(:, mi) = FL_outside_probabilities_foldAveraged; % Now take the averaged probability
            end
            disp(['============== First-layer MDM end on node ', mat2str(ti),' =============='])
            
            % Second Layer MD training
            disp(['============== Second-layer MDM start on node ', mat2str(ti),' =============='])
            % Model Estimation
            [SL_inside_ConfusionMatrix, SL_inside_predictions, SL_inside_Coefficients, SL_inside_FitInfo, SL_inside_CrossValSet_0, SL_inside_CrossValSet_1, SL_inside_probabilities, SL_inside_Deviance] = ModelFit(observes_secondLayer, double(thisTrial_c));
            
            % Save second layer intermediate results
            secondR.SL_inside_ConfusionMatrix = SL_inside_ConfusionMatrix; % confusion matrix
            secondR.SL_inside_MCC = mcc(SL_inside_ConfusionMatrix); % MCC
            secondR.SL_inside_predictions = SL_inside_predictions; % prediction - out Sample
            secondR.SL_inside_Coefficients = SL_inside_Coefficients; % Fiting Results of each fold - Coefficients
            secondR.SL_inside_FitInfo = SL_inside_FitInfo; % Fiting Results of each fold - Fitting Information
            secondR.SL_inside_CrossValSet_0 = SL_inside_CrossValSet_0; % Patition setting for label 0
            secondR.SL_inside_CrossValSet_1 = SL_inside_CrossValSet_1; % Patition setting for label 1
            secondR.SL_inside_probabilities = SL_inside_probabilities; % prediction - probability
            secondR.SL_inside_Deviance = SL_inside_Deviance;
            
            % Second Layer MD Nested testing
            % Test with the self CV method
            % treat every model as different trals
            % Pick the global lambda
            globalMinDeviance2 = min(SL_inside_Deviance);
            globalIndex2 = find(SL_inside_Deviance == globalMinDeviance2);
            if length(globalIndex2) > 1
                globalIndex2 = globalIndex2(1);
            end
                
            numFold = length(SL_inside_Coefficients);
            SL_outside_probabilities = zeros(numFold, length(testingTarget));
            SL_outside_predictions = zeros(numFold, length(testingTarget));
            SL_outside_GlobalCoefficients = zeros(numFold, size(SL_inside_Coefficients{1}, 1));
            SL_outside_GlobalC0 = zeros(numFold, 1);
            for tempI = 1:numFold
                
                globalC0 = SL_inside_FitInfo{tempI}.Intercept(globalIndex2);
                SL_outside_GlobalC0(tempI) = globalC0;
                SL_outside_GlobalCoefficients(tempI, :) = SL_inside_Coefficients{tempI}(:, globalIndex2);
                c_i2 = observes_secondLayer_testing * SL_inside_Coefficients{tempI}(:, globalIndex2) + globalC0;
                c_p2 = 1 ./ (1 + exp(-c_i2));
                SL_outside_probabilities(tempI, :) = c_p2;
                SL_outside_predictions(tempI, :) = double(c_p2>=0.5); % Can adjust this threshold latter
            end
            
            % Calculate the average performance
            SL_outside_probabilities_modelAveraged = mean(SL_outside_probabilities, 1);
            SL_outside_predictions_modelAveraged = double(SL_outside_probabilities_modelAveraged>=0.5);
            SL_outside_CM_modelAveraged = confusionmat(testingTarget, SL_outside_predictions_modelAveraged);
            if (size(SL_outside_CM_modelAveraged,1)==1&&size(SL_outside_CM_modelAveraged,2)==1)
                SL_outside_CM_modelAveraged = [SL_outside_CM_modelAveraged(1,1) 0;0 0];
            end
            SL_outside_MCC_modelAveraged = mcc(SL_outside_CM_modelAveraged);
            
            % Save secondMD outside intermediate results
            secondR.SL_outside_probabilities = SL_outside_probabilities;
            secondR.SL_outside_predictionss = double(SL_outside_probabilities>=0.5);
            secondR.SL_outside_MCC_modelAveraged = SL_outside_MCC_modelAveraged;
            secondR.SL_outside_predictions_modelAveraged = double(SL_outside_predictions_modelAveraged);
            secondR.SL_outside_probabilities_modelAveraged = SL_outside_probabilities_modelAveraged;
            secondR.testingTarget = testingTarget;
            
            disp(['============== Second-layer MDM end on node ', mat2str(ti),' =============='])
        end
        
        function [firstR, secondR] = runAllSplits(obj,varargin)
            
            if obj.par
                parfor ti = 1:obj.num_split
                    [R_sepSplit{ti}, R2_sepSplit{ti}] = obj.runSplit(ti);
                end
            else
                R_sepSplit = cell(obj.num_split, 1);
                R2_sepSplit = cell(obj.num_split, 1);
                for ti = 1:obj.num_split
                    [R_sepSplit{ti}, R2_sepSplit{ti}] = obj.runSplit(ti);
                end
            end
            
            for ti = 1:obj.num_split
                for mi = 1:length(obj.m_all)
                    firstR(ti,mi) = R_sepSplit{ti}(mi);
                end
                secondR(ti) = R2_sepSplit{ti};
            end
            
        end
        
        function DLMDfit = run(obj,varargin)
            if obj.par
                poolOb = obj.setupPool;
            end
            tStart = tic;
            [obj.R_first, obj.R_second] = obj.runAllSplits;
            obj.tFit = toc(tStart);
            DLMDfit = obj;
            save(obj.oF, 'DLMDfit', '-v7.3')
            if obj.par
                poolOb.delete;
            end
        end
        
        % Here, the default is to put a seperate node per bagging split, 
        % with the number of workers equal to the number of splits
        function poolOb = setupPool(obj)
            if ~isunix % This will run on your local machine using each cpu core as one node
                poolOb = parpool;
                
            else  % This will run on the cluster
                nWorkers = obj.num_split;
                
                % This setting should depend on your corresponding
                % hpc/cluster settings
                clusProf = get_SLURM_cluster('/home/rcf-proj/tb/xiweishe/matlab_storage','/usr/usc/matlab/R2018a/SlurmIntegrationScripts','--time=50:00:00 --partition berger');
                poolOb = parpool(clusProf,nWorkers);
            end
        end
    end
    
end

