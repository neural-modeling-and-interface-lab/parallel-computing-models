function [loss,P,Pf] = ProbitLoss(w,X,y,varargin)
%ProbitLoss calculates the negative loglikeihood for a probit GLM given a
% covariate matrix, estimated parameters, and observations
% Inputs:
%   w is parameters or estimated coefficients (N x 1)
%   X is covariates or design matrix TxN
%   y is observations Tx1  (1's or 0's)
% 
% Outputs:
%   loss: negative loglikelihood (not normalized by number of samples)
%   P: likelihood of each observation
%   Pf: firing probability estimated by the model
%
% author: Dong Song, modified by Brian Robinson to use 0 or 1 for y
% observations 2014-2016, uses normcdf to estimate probability instead of
% erf

[n,p] = size(X);
[offset] = process_options(varargin,'offset',zeros(n,1));  %default is to have no offset!!
if p == (length(w)-1)  %if there are no column of ones included in X for the constant offset term, this will add them!
    X = [ones(n,1) X];
end
Xw = X*w+offset;

Pf = normcdf(Xw); %calculate estimated firing probability
P = Pf.*(y==1) + (1-Pf).*(y==0); %calculate likelihood of each observation

loss = -sum(log(P)); %calculate negative log likelihood of all observations
