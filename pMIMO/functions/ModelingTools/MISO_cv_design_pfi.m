function f = MISO_cv_design_pfi(x,y,Bff,Bfb,Qff,Qfb,Optff,Optfb,varargin)
% performs K-fold crossvalidation for MISO estimation
%   Inputs:
%       Optional Name Value Pairs:
%           K
%           cv_seed
%           g_include
%           LCD
%           LCD_max_it
%           link
%   Output fields:
%       inSampleFits, K dimensional struct array with in-sample
%           MISO_fit_obs
%       outSampleLoss, Kx1, out of sample average negative log likelihood
%       KS, struct containing fields which are all Kx1 that summarize out
%           of sample KS test
%
% author: Brian Robinson, 2014-2016
% Modified by Xiwei She to include Permutation Feature Importance (pfi)
% Last Modified Date: Dec. 17, 2021

[K, cv_seed, g_include, LCD, LCD_max_it,link, pfi]=process_options(varargin,'K',5,'cv_seed',9013,'g_include',-1,'LCD',0,'LCD_max_it',200,'link','probit', 'pfi', 1);  %sets default fold numbers for cross validation and seed for randomly splitting indices for CV
T = length(y);
[V, ~, ff_inds, fb_inds] = createDesign(x,y,Bff,Bfb,Qff,Qfb,Optff,Optfb,'g_include',g_include);
rng(cv_seed);
cv_inds = crossvalind('Kfold',T,K); %get indices to use for cross validaion

% yVal_allFolds = cell(K, 1);
% tr_inds_allFolds = cell(K, 1);
% val_inds_allFolds = cell(K, 1);
% c_tilde_allFolds = cell(K, 1);

% Variables for the permutation feature importance check
if pif == 1
    groupSize = Bff.getN + Bff.getN * (Bff.getN + 1) / 2; % Each input corresponds to how many features
    numInput = size(V, 2) / groupSize - 1; % Exclude the output feedback feature
    permuteTime = 20;
    outSampleLoss_permute = zeros(K, numInput);
end

for i_fold = 1:K %lop through folds
    disp(['Starting Fold ' num2str(i_fold)])
    val_inds = find(cv_inds == i_fold);  %out-sample time indices
    tr_inds = find(cv_inds ~= i_fold);  %in-sample time indices
    %% Perform MISO on in-sample data
    V_tr = V(tr_inds,:);    %in-sample design matrix
    y_tr = y(tr_inds);      %in-sample training observations
    MISO_fit_ob = MISO_fromV(V_tr,y_tr,ff_inds,fb_inds,'LCD',LCD,'LCD_max_it',LCD_max_it,'link',link);
    MISO_fit_ob.Bff=Bff;
    MISO_fit_ob.Bfb=Bfb;
    cv_fit(i_fold)=MISO_fit_ob;
    
    %% calculate out of sample predicted firing probability
    V_val = V(val_inds,:);
    y_val = y(val_inds);
    c_tilde = cv_fit(i_fold).c_tilde;
    c = [c_tilde.c_0; c_tilde.k; c_tilde.h];
    P = glmval(c,V_val,link);
    
    %% calculate negative log likelihood and KS test on out-of sample
    %probability prediction
    [Z{i_fold}, b{i_fold}, b95{i_fold}, KS_score(i_fold)]= discretetimerescaling(P, y_val);
    
%     yVal_allFolds{i_fold} = y_val;
%     tr_inds_allFolds{i_fold} = tr_inds;
%     val_inds_allFolds{i_fold} = val_inds;
%     c_tilde_allFolds{i_fold} = c_tilde;
    
    outSampleLoss(i_fold) = ProbitLoss(c, V_val, y_val)./length(y_val);
    if isinf(outSampleLoss(i_fold))  %% if an infinite value is reached, we will not bother doing the rest of the folds, this saves computation time
        break
    end
    
    %% Permutation Feature Importance - Added by Xiwei Dec. 17, 2021
    if pfi == 1
        for n = 1:numInput
            permuteGroup = [ (n-1)*groupSize + 1 :  n*groupSize ];
            outSampleLoss_permute_temp = zeros(permuteTime, 1);
            for i_permute = 1:permuteTime
                V_val_permute = V_val;
                permuteSeed = randperm(size(V_val_permute, 1));
                V_val_permute(:, permuteGroup) = V_val_permute(permuteSeed, permuteGroup);
                outSampleLoss_permute_temp(i_permute) = ProbitLoss(c, V_val_permute, y_val)./length(y_val);
            end
            outSampleLoss_permute(i_fold, n) = mean(outSampleLoss_permute_temp);
        end
    else
        outSampleLoss_permute = [];
    end
    
end
%store results in struct
f.inSampleFits = cv_fit;
f.outSampleLoss = outSampleLoss;
f.KS.Z = Z;
f.KS.b = b;
f.KS.b95 = b95;
f.KS.KS_score = KS_score;
f.outSampleLoss_permute = outSampleLoss_permute;

% f.V = V;
% f.ff_inds = ff_inds;
% f.fb_inds = fb_inds;
% f.yVal_allFolds = yVal_allFolds;
% f.tr_inds_allFolds = tr_inds_allFolds;
% f.val_inds_allFolds = val_inds_allFolds;
% f.c_tilde_allFolds = c_tilde_allFolds;