%% This is the main script of the MIMO model Estimation

% Please cite: 
%   Song, D., Chan, R. H., Marmarelis, V. Z., Hampson, R. E., 
%   Deadwyler, S. A., & Berger, T. W. (2007). Nonlinear dynamic modeling 
%   of spike train transformations for hippocampal-cortical prostheses. 
%   IEEE Transactions on Biomedical Engineering, 54(6), 1053-1066.

% Input: Data from MIMOPreprocessing.m
% Output: save to file

% Author: Dong Song, 
%         Brian Robinson, 2014-2016
%         Xiwei She, 2016-2022

clear;clc;

% Define case information
localRun = 1; % change to 0 if you want to use cluster
fileLabel='fit example'; % Define your fitting case name

%% initialize all paths and make sure they all get added to the cluster
p1 = fullfile('functions');
addpath(genpath(p1));

DataFolder = 'exampleData\exampleInput\';
DataName = 'exampleData'; 

dataFileName = strcat(DataFolder,DataName,'_ChannelUnit.mat');
tLoad = load([DataFolder,DataName,'_Correct_Intervals.mat']);

T_trunc = tLoad.correct_intervalsALL;
dataFiles = {dataFileName};
p1Files = dir_matlab_files(p1);
extraFiles = [p1Files  dataFiles];

%% Cluster Settings
max_nodes = 150;    %only used in balance parallel loads
new_clus = 1;

%% Data Selection Settings
%Define data parameters
inRegion = 'CA3';
outRegion = 'CA1';
maxInputCells = Inf; % May changed in the quick debug mode
maxOutputCells = Inf;% May changed in the quick debug mode
folderN = fileLabel;

%% Model Settings
%Define model creation parameters
M = 1000;
L = 3;
a_ff = .97;
a_fb = .95;
modelOptions.Bff = basis('lagSeries',M,L,a_ff);
modelOptions.Bfb = basis('lagSeriesZS',M,L,a_fb);
modelOptions.Qff = 2;
modelOptions.Qfb = 2;
modelOptions.MLE_LCD = 0;
modelOptions.LCD_max_it = nan;
modelOptions.bin_size = 2;

%% Regularization fitting parameters
n_lams = 30;
regOptions.link = 'probit';
regOptions.max_it = 50; % The number of max iteration % May changed in the quick debug mode
regOptions.n_lams = n_lams; % The number of lambda % May changed in the quick debug mode

%% Regression Parameters
regOptions.min_factor = 1e-2;
regOptions.coeff_tolerance = 1e-4;
regOptions.cv_K = 0;
regOptions.cv_seed = 99;
regOptions.mle_refitting = 0;
regOptions.K=5;  %this controls the folds of cross validation!
regOptions.cv_design = 1;  % does correct cross-validation
fit_reMLE_CV=1;
fit_CV_criteria = 1;

%% Create Cell groups and place to save everything
all_cells = cellGroup(dataFiles);

%place to save
resultFolder = fullfile('exampleData','exampleResult',folderN);
if ~exist(resultFolder,'dir')
    mkdir(resultFolder)
end

%filter output cells
outputCells=all_cells;
outputCells=outputCells.filterGroup('region',outRegion);
if maxOutputCells<= length(outputCells.cellArray)
    outputCells=cellGroup(outputCells.cellArray(1:maxOutputCells));
end
Nout = length(outputCells.cellArray);

%filter input cells
inputCells=all_cells;
inputCells=inputCells.filterGroup('region',inRegion);
if maxInputCells <= length(inputCells.cellArray)
    inputCells = cellGroup(inputCells.cellArray(1:maxInputCells));
end
warning off 'parallel:cluster:CannotSaveCorrectly'

%% Claim Job Parameters
Reg_fits = cell(length(Nout));
j = cell(length(Reg_fits));

%% Define the Fits
for i = 1:Nout
    yCell = outputCells.cellArray(i);
    Reg_fits{i} = Reg_fit(inputCells,yCell,modelOptions,regOptions,fileLabel,'T_trunc',T_trunc,'fit_reMLE_CV',fit_reMLE_CV,'fit_CV_criteria',fit_CV_criteria);
end

%% Start MIMO Estimation
if localRun==1 % Run locally
    for i = 1:length(Reg_fits)
        disp(['Current running batch fit: ', mat2str(i)])
        Reg_fit = run_batch_fit(Reg_fits{i},resultFolder);
    end
    
else % Run parallel on cluster
    
    N_ret_var = 0;
    for i = 1:length(Reg_fits)
        tic
        vars = {Reg_fits{i},resultFolder};
        % This part need to be modified according to the hpc settings
        j{i} = batch('run_batch_fit',N_ret_var,vars,'CaptureDiary',true,'AttachedFiles',extraFiles,'AutoAttachFiles',false);   %for batch mode
        disp(['Submitted ' num2str(i) ' out of ' num2str(length(Reg_fits))])
        toc
    end
end