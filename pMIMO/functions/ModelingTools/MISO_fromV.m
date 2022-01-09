function MISO_fit_ob = MISO_fromV(V,y,ff_inds,fb_inds,varargin)
%performs MISO identifaction given a design matrix.  
%Inputs:
%   ff_inds and fb_inds designate which design matrix columns are
%       feedforwad and feedback respectively 
%        ([ff_inds;  fb_inds]'=1:size(V,2)
%Output fields:
%   glmfit returned fields
%     MISO_fit_ob.sig
%     MISO_fit_ob.dev
%     MISO_fit_ob.stats
%   LCD returned fields
%     MISO_fit_ob.convergence_max_diff
%     MISO_fit_ob.converged
%     MISO_fit_ob.max_iterated
%   Coefficients (Parameterized
%     MISO_fit_ob.c_tilde (un-normalizezd, grouped coefficients)
%     MISO_fit_ob.c  (normalized, grouped coefficients)
%     MISO_fit_ob.c_temp (un-normalized, vector of coefficients, the result of regression on the design matrix)
%   Additional Fit Information
%     MISO_fit_ob.loss (negative log likelihood (not normalized by samples)
%     MISO_fit_ob.T (number of observations)
%     MISO_fit_ob.Telapsed (time elapsed in seconds)
%     MISO_fit_ob.link, (link funciont used for fitting)
%
% author: Brian Robinson, 2014-2016
T = size(V,1);
[offset, LCD, LCD_max_it,link] = process_options(varargin,'offset',zeros(T,1),'LCD',0,'LCD_max_it',200,'link','probit');  %default is to have no offset

Tstart = tic;
if ~LCD
    [c_temp, dev, stats] = glmfit2(V,y,'binomial','link',link,'offset',offset);
    switch link
        case 'probit'
            loss = ProbitLoss(c_temp, V, y);
        otherwise
            loss = nan;
            warning('Loss Calculation not added for logit yet!')
    end
else
    lam = 0;
    method='groupLasso';
    c_init = zeros(size(V,2)+1,1);
    g = [ones(size(V,2),1)];
    coeff_tolerance = 1e-4;
    [c_temp, ~, loss, ~, ~, ~, ~, ~, ~, convergence_max_diff, converged, max_iterated] = groupGLM_indRupdate(V,y,g,c_init,method,link,lam,'max_it',LCD_max_it,'coeff_tolerance',coeff_tolerance);
end
Telapsed = toc(Tstart);

c_tilde.c_0 = c_temp(1); 
c_tilde.k = c_temp(ff_inds+1); 
c_tilde.h = c_temp(fb_inds+1);
if isempty(c_tilde.k)
    c_tilde.k=[];
end
if isempty(c_tilde.h)
    c_tilde.h=[];
end
%normalize estimated coefficients
c.sig = -1/c_tilde.c_0;
c.h = c_tilde.h*c.sig;
c.k = c_tilde.k*c.sig;
sig = c.sig;

if ~LCD
    stats = rmfield(stats,{'resid','residp','residd','resida','wts'});
    stats.covb=diag(stats.covb);
    stats.coeffcorr=diag(stats.coeffcorr);
    MISO_fit_ob.sig=sig;
    MISO_fit_ob.dev=dev;
    MISO_fit_ob.stats=stats;
else
    MISO_fit_ob.convergence_max_diff=convergence_max_diff;
    MISO_fit_ob.converged=converged;
    MISO_fit_ob.max_iterated=max_iterated;
end
MISO_fit_ob.c_tilde=c_tilde;
MISO_fit_ob.c=c;
MISO_fit_ob.loss=loss;
MISO_fit_ob.T=T;
MISO_fit_ob.c_temp = c_temp;
MISO_fit_ob.Telapsed = Telapsed;
MISO_fit_ob.link = link;