function [fit_data] = groupLassoPath(X,y,g,link,max_it, n_lams,min_factor,coeff_tolerance,cv_K,cv_seed,mle_refitting)
%this function fits lambda values along a path
%n_lams and min_factor determine the path
%max_it, coeff_tolerance determine how many iterations the algorithm does.
%cv_K specifies the number of cross validation folds.  If set to 0, cross validation is not performed, (This is a deprecated way of calculating the cross validation!)
%cv_seed specifies the random seed that determines how the indices are randomly divided for the cross validation trials
%   
%   author: Brian Robinson, 2014-2016

P = size(X,2);
method = 'groupLasso';

%initialize lambda path!
lams = setupLassoLamba(X,y,g,n_lams,min_factor,link);

%intialize cells where each lambda data will be stored
c_traces = cell(n_lams,1);
sparse_groups_traces = cell(n_lams,1);
convergence_max_diffs = nan(max_it,n_lams);
convergeds = nan(n_lams,1);
max_iterateds = nan(n_lams,1);
criterias = cell(n_lams,1);
penalties = cell(n_lams,1);
objective_traces  = cell(n_lams,1);
sparse_groups=nan(max(g),n_lams);
cs=nan(P+1,n_lams);
%initialize c values to all zero!
c = zeros(P+1,1);

if cv_K>0
    cv_N = length(y);
    rng(cv_seed);
    cv_inds = crossvalind('Kfold',cv_N,cv_K);
    %initilialize CV data storage
    av_Loss_CV_MLEs=zeros(n_lams,1);
    av_Loss_CV_Regs=zeros(n_lams,1);
    c_av_CV_MLEs= zeros(P+1,n_lams);
    c_av_CV_Regs= zeros(P+1,n_lams);
end

%loop through lambda values starting at lambda_max and going down to lambda
%min
for i_l = 1:length(lams)
    disp(['Starting eval of lambda ' num2str(i_l) ' out of ' num2str(length(lams))])
    lam=lams(i_l);
    c_init=c;
    try
        disp(['Memory Ratio is ' num2str(memRatio) ', before this lambda iteration'])
    end
    [c, M, L, c_trace, L_trace, objective_trace, sparse_groups_trace, penalty, criteria, convergence_max_diff, converged, max_iterated] = groupGLM_indRupdate(X,y,g,c_init,method,link,lam,'max_it',max_it,'coeff_tolerance',coeff_tolerance);
    c_traces{i_l} = c_trace;
    sparse_groups_traces{i_l} = sparse_groups_trace;
    final_sparse_groups = sparse_groups_trace(max_iterated,:);
    sparse_groups(:,i_l)=final_sparse_groups';
    convergence_max_diffs(:,i_l)=convergence_max_diff;
    convergeds(i_l)=converged;
    max_iterateds(i_l)=max_iterated;
    criterias{i_l}=criteria;
    penalties{i_l}=penalty;
    objective_traces{i_l}= objective_trace;
    cs(:,i_l)=c;
    if cv_K>0  %This is an old way of calculation Cross-Validation!
        [av_Loss_CV_MLE, av_Loss_CV_Reg, c_av_CV_MLE, c_av_CV_Reg] = groupGLM_CV(X,y,g,c_init,method,link,lam,cv_inds,max_it,coeff_tolerance);
        av_Loss_CV_MLEs(i_l) = av_Loss_CV_MLE;
        av_Loss_CV_Regs(i_l) = av_Loss_CV_Reg;
        c_av_CV_MLEs(:,i_l) = c_av_CV_MLE;
        c_av_CV_Regs(:,i_l) = c_av_CV_Reg;
    end
    if mle_refitting
        c_MLE(:,i_l) = mle_refit(X,y,c,link);
    end
end
fit_data.lams = lams;
fit_data.c_traces=c_traces;
fit_data.sparse_groups_traces=sparse_groups_traces;
fit_data.convergence_max_diffs=convergence_max_diffs;
fit_data.convergeds=convergeds;
fit_data.max_iterateds=max_iterateds;
fit_data.criterias=criterias;
fit_data.penalties=penalties;
fit_data.objective_traces=objective_traces;
fit_data.sparse_groups=sparse_groups;
fit_data.cs=cs;
fit_data.fit_options.tol=coeff_tolerance;
fit_data.fit_options.n_lams=n_lams;
fit_data.fit_options.max_it=max_it;
fit_data.fit_options.min_factor=min_factor;
if cv_K>0
    fit_data.CV.av_Loss_CV_MLEs=av_Loss_CV_MLEs;
    fit_data.CV.av_Loss_CV_Regs=av_Loss_CV_Regs;
    fit_data.CV.c_av_CV_MLEs=c_av_CV_MLEs;
    fit_data.CV.c_av_CV_Regs=c_av_CV_Regs;
end
if mle_refitting
    fit_data.c_MLE = c_MLE;
end