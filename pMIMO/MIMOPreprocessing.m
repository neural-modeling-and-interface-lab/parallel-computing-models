%% This script is used to preprocess neural signals for the MIMO modeling

% Please cite: 
%   Song, D., Chan, R. H., Marmarelis, V. Z., Hampson, R. E., 
%   Deadwyler, S. A., & Berger, T. W. (2007). Nonlinear dynamic modeling 
%   of spike train transformations for hippocampal-cortical prostheses. 
%   IEEE Transactions on Biomedical Engineering, 54(6), 1053-1066.

% Input: Neural data
% Output: save to file

% Author: Dong Song, 
%         Brian Robinson, 2014-2016
%         Xiwei She, 2016-2022

clear all;clc;
addpath(genpath('function'))
iF = 'exampleData\exampleInput';

%% Preprocessing
load(strcat(iF, '\exampleData_neural.mat'));

% Bin each behavioral events into intervals
binsize = 2; %specified in ms
correct_intervalsALL = [];
for i=1:length(SUCCESS) % We only use correct trials
    
    preceding_ALL_S_PHASE = ALL_S_PHASE(ALL_S_PHASE<SUCCESS(i));
    if ~isempty(preceding_ALL_S_PHASE) % Fixed by Xiwei
        preceding_ALL_S_PHASE = preceding_ALL_S_PHASE(end);
    else
        preceding_ALL_S_PHASE = 0;
    end
    
    after_ALL_S_PHASE = ALL_S_PHASE(ALL_S_PHASE>SUCCESS(i));
    if i==length(SUCCESS) && numel(after_ALL_S_PHASE)==0 %for the last trial
        try
            after_ALL_S_PHASE = STOP;
        catch
            after_ALL_S_PHASE = SUCCESS(end);
        end
    else
        after_ALL_S_PHASE = after_ALL_S_PHASE(1);
    end
    
    after_ITI_ON_B = after_ALL_S_PHASE*1000/binsize;
    preceding_ITI_ON_B = preceding_ALL_S_PHASE*1000/binsize;
    correct_intervalsALL = [correct_intervalsALL {[preceding_ITI_ON_B, after_ITI_ON_B]}];
    intervalMaxTs(i) = after_ALL_S_PHASE;
end

savePath = strcat(iF, '\exampleData_Correct_Intervals.mat');
save(savePath,'correct_intervalsALL');

%% Match neurons to CA1 and CA3 channels
CA3_Chan = [1,4,6,8,9,13,14,15,22,25,26,29,30,31,32, 97,98,101,108,109,111];
CA1_Chan = [33,34,35,38,45,46,47,49,53,55,56,57,61,62,74,85,90,91,92];

L_channels = [1:62];
R_channels = [85:111];

OutputIndices = CA1_Chan;
source_indices = CA3_Chan;

savePath = strcat(iF, '\exampleData_ChannelUnit.mat');
save(savePath);