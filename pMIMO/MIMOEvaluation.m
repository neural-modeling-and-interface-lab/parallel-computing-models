%% This script is used to visualize the modeling performance of MIMO model

% Please cite: 
%   Song, D., Chan, R. H., Marmarelis, V. Z., Hampson, R. E., 
%   Deadwyler, S. A., & Berger, T. W. (2007). Nonlinear dynamic modeling 
%   of spike train transformations for hippocampal-cortical prostheses. 
%   IEEE Transactions on Biomedical Engineering, 54(6), 1053-1066.

% Input: Data from MIMOEstimation.m
% Output: save to file

% Author: Xiwei She, 2016-2022

clear;clc

p1 = fullfile('toolbox');
addpath(genpath(p1));

fols{1} = 'fit example'; % This is the output folder of the MIMOEstimation.m

% This is the title name, suggested keep same as above 'fols'
session_names{1} = fols{1};

%% Create Initial Meta Data For all fits!
for i=1:length(fols)
    createInitialMetaData(fols{i});
end

%% Plot Connectivity Matrix Comparisons
figure()
for i=1:length(fols)
    plotConnectivityMatrix(fols{i})
    title(session_names{i})
end

%% Plot KS Metric Comparison
FreqOrder = 0;
figure()
for i=1:length(fols)
    comparison = 'none';
    criterias = {'CV','zeroth'};
    sessionLabels = criterias;
    metric = 'KS';
    stages = [3 3];
    labelY = 1;
    plotMetricComparison(fols([i i]),sessionLabels,metric,criterias,comparison,stages,'FreqOrder',FreqOrder,'labelY',labelY);
    title(session_names{i});
    title(session_names{i});
    xlabel('Output Cell');
    ylabel('KS Score');
end

%% Plot Log Lik Metric Comparison
FreqOrder = 0;
figure()
for i=1:length(fols)
    comparison = 'vs';
    criterias = {'CV','zeroth'};
    sessionLabels = criterias;
    metric = 'Loss';
    stages = [3 3];
    labelY = 1;
    plotMetricComparison(fols([i i]),sessionLabels,metric,criterias,comparison,stages,'FreqOrder',FreqOrder,'labelY',labelY);
    title(session_names{i});
    title(session_names{i});
    xlabel('Output Cell');
    ylabel('Log Likelihood ratio');
end