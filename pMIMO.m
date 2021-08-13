% This runs the MIMO Model with parallel computing strategies on USC
% CARC Cluster
    
% Author: Xiwei She, Brian Robinson

clearvars
%% initialize all paths and make sure they all get added to the cluster
neuralData = 'Neural_Data.mat';
behavioralData = 'Behavioral_Data.mat';
neuralDataFiles = {neuralData};
behavioralDataFiles = {behavioralData};

extraFiles = [behavioralDataFiles  neuralDataFiles];

%% Cluster Settings
max_nodes = 150;    %only used in balance parallel loads
new_clus = 1;

%% Model Settings
%Define model creation parameters
M=1000;
L=3;
a_ff=.97;
a_fb = .95;
modelOptions.Bff = basis('lagSeries',M,L,a_ff);
modelOptions.Bfb = basis('lagSeriesZS',M,L,a_fb);
modelOptions.Qff = 2;
modelOptions.Qfb = 2;
modelOptions.MLE_LCD = 0;
modelOptions.LCD_max_it = nan;
modelOptions.bin_size = 2;

% Regularization fitting parameters
n_lams = 30; % Number of Lambda (i.e., the lambda loop)
regOptions.link = 'probit';
regOptions.max_it = 50; % The number of max iteration 
regOptions.n_lams = n_lams; % The number of lambda

%% Regression Parameters
regOptions.min_factor = 1e-2;
regOptions.coeff_tolerance = 1e-4;
regOptions.cv_K = 0;
regOptions.cv_seed = 99;
regOptions.mle_refitting = 0;
regOptions.K=5;  % the CV loop
regOptions.cv_design = 1;  % does correct cross-validation


%% Create Cell groups and place to save everything
all_cells = cellGroup(neuralDataFiles);

%place to save
resultFolder = fullfile('..','results_human',folderN);
if ~exist(resultFolder,'dir')
    mkdir(resultFolder)
end

%filter output cells
outputCells=all_cells;
outputCells=outputCells.filterGroup('region',outRegion);
if maxOutputCells<= length(outputCells.cellArray)
    outputCells=cellGroup(outputCells.cellArray(1:maxOutputCells));
end
Nout = length(outputCells.cellArray);% number of output neurons

%filter input cells
inputCells=all_cells;
inputCells=inputCells.filterGroup('region',inRegion);
if maxInputCells <= length(inputCells.cellArray)
    inputCells = cellGroup(inputCells.cellArray(1:maxInputCells));
end

%% Claim Job Parameters
MIMOModelEstimation = cell(length(Nout));
j = cell(length(MIMOModelEstimation));

%% Define the Fits
for i = 1:Nout % the output loop
    yCell = outputCells.cellArray(i);
    MIMOModelEstimation{i} = Reg_fit(inputCells,yCell,modelOptions,regOptions);
end

%% Distribute modeling tasks to the cluster
if localRun % Run non-parallel
    Reg_fit = run_batch_fit(MIMOModelEstimation{2},resultFolder);
else % Run parallel!
    N_ret_var = 0;
    for i = 1:length(MIMOModelEstimation)
        tic
        vars = {MIMOModelEstimation{i},resultFolder};
        j{i} = batch('run_batch_fit',N_ret_var,vars,'CaptureDiary',true,'AttachedFiles',extraFiles,'AutoAttachFiles',false);   %for batch mode
        disp(['Submitted ' num2str(i) ' out of ' num2str(length(MIMOModelEstimation))])
        toc
    end
end