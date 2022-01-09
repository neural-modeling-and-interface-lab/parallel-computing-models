function [CM_outSample, c_predict_outSample, B_outSample, FitInfo_outSample, CrossValSet_0, CrossValSet_1, c_probability_outSample, theDeviance_outSample_fold] = ModelFit(inputFeature, target)
%% This function used to separate the data set into Training and Testing Sets which contain equal number of 1 and 0 labels
% And use this method to compare the fitting parameters lambda from overall
% data set with picked lambda from each fold
% Created by Xiwei 2017-9-18

warning off;

%% Do the classification of each fold to make pure out of sample classification
lambdaPool = power(10, 1:-0.25:-5); % Define the lambda pool

zeroIndex = find(target==0); % Find those 0 labels
oneIndex = find(target==1); % Find those 1 labels

% Decide which group is the minority set
loopNum = min(length(zeroIndex), length(oneIndex));
if loopNum > 10
    loopNum = 10;
end

% Label 0 group
zeroSetC = target(zeroIndex, :);
zeroSetP = inputFeature(zeroIndex, :);

% Label 1 group
oneSetC = target(oneIndex, :);
oneSetP = inputFeature(oneIndex, :);

if ((length(zeroIndex)==1) || (length(oneIndex)==1) || (isempty(zeroIndex)) || (isempty(oneIndex)))
    CM_outSample = 0; c_predict_outSample = zeros(length(target),1); B_outSample = []; FitInfo_outSample = [];
    CrossValSet_0 = []; CrossValSet_1 = []; c_probability_outSample = zeros(length(target),1); theDeviance_outSample_fold = [];
    
    if (size(CM_outSample,1)==1 && size(CM_outSample,2)==1)
        CM_outSample = [CM_outSample(1,1) 0;0 0];
    end
    return;
else
    CrossValSet_0 = cvpartition(length(zeroSetC),'KFold', 10);
    CrossValSet_1 = cvpartition(length(oneSetC),'KFold', 10);
end

% Predict label
c_predict_outSample = zeros(length(target), 1);
c_probability_outSample = zeros(length(target), 1);

% Store all B and FitInfo from each fold
B_outSample = cell(loopNum, 1);
FitInfo_outSample = cell(loopNum, 1);

for loopFold = 1:loopNum
    % Label 0 Data Separation
    P0_train = zeroSetP(training(CrossValSet_0, loopFold),:);
    c0_train = zeroSetC(training(CrossValSet_0, loopFold),:);
    P0_test = zeroSetP(test(CrossValSet_0, loopFold),:);
    c0_test = zeroSetC(test(CrossValSet_0, loopFold),:);
    index0 = test(CrossValSet_0, loopFold)==1;
    
    % Label 1 Data Separation
    P1_train = oneSetP(training(CrossValSet_1, loopFold),:);
    c1_train = oneSetC(training(CrossValSet_1, loopFold),:);
    P1_test = oneSetP(test(CrossValSet_1, loopFold),:);
    c1_test = oneSetC(test(CrossValSet_1, loopFold),:);
    index1 = test(CrossValSet_1, loopFold)==1;
    
    % Find the original index in c for those testing label
    originalCIndex = [zeroIndex(index0); oneIndex(index1)];
    
    % Training Set
    P_train = [P0_train;P1_train];
    c_train = [c0_train;c1_train];
    
    % Testing Set
    P_test = [P0_test;P1_test];
    c_test = [c0_test;c1_test];
    
    % Lassoglm
    [B_outSample_fold,FitInfo_outSample_fold] = lassoglm(P_train,double(c_train),'binomial', 'Lambda', lambdaPool, 'MaxIter', 1e3);
    B_outSample{loopFold} = sparse(B_outSample_fold);
    FitInfo_outSample{loopFold} = FitInfo_outSample_fold;
    C_v = B_outSample_fold; % selected coefficients
    C0 = FitInfo_outSample_fold.Intercept;
    theDeviance_outSample_fold = zeros(size(C_v, 2), 1);
    c_p_outSample_fold = zeros(size(C_v, 2), size(c_test, 1));
    
    for loopLambda = 1:size(C_v, 2) % Loop through all lambda
        c_i_outSample = P_test * C_v(:, loopLambda) + C0(loopLambda);
        c_p_outSample_fold(loopLambda, :) = 1 ./ (1 + exp(-c_i_outSample));
        probability_outSample_temp = (1 - c_test) + 2 * (c_test - 0.5) .* c_p_outSample_fold(loopLambda, :)';
        theDeviance_outSample_fold(loopLambda) = -2 * sum(log(probability_outSample_temp));
    end
   
    % Then find the minimum deviance and it's parameters
    minDeviance = min(theDeviance_outSample_fold);
    indexDev = find(theDeviance_outSample_fold == minDeviance);
    if length(indexDev) > 1
        indexDev = indexDev(1);
    end
    sc_p = c_p_outSample_fold(indexDev, :);
    c_predict_fold = double(sc_p > 0.5);
    
    c_predict_outSample(originalCIndex) = c_predict_fold;
    c_probability_outSample(originalCIndex) = sc_p;
    
end

CM_outSample = confusionmat(double(target),c_predict_outSample); % confusion matrix

if (size(CM_outSample,1)==1 && size(CM_outSample,2)==1)
    CM_outSample = [CM_outSample(1,1) 0;0 0];
end