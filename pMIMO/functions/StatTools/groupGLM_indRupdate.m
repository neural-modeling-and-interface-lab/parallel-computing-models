function [c, M, L, c_trace, L_trace, objective_trace, sparse_groups_trace, penalty, criteria, convergence_max_diff, converged, max_iterated] = groupGLM_indRupdate(X,y,g,c_init,method,link,lambda,varargin)
% performs group regressino for a particular lambda value
%   
%   author: Brian Robinson, 2014-2016
try
    disp(['Memory Ratio is ' num2str(memRatio) ', at beginning of groupGLM estimation'])
end
[group_bridge_gamma, max_it, loss_tolerance, delta_stability,coeff_tolerance,verbose]=process_options(varargin,'group_bridge_gamma',.5,'max_it',20,'loss_tolerance',1e-3,'delta_stability',1e-8,'coeff_tolerance',eps,'verbose',1);
%default delta_stability of 1e-8 is chosen to match that of grpreg
T = size(X,1);
P = size(X,2);

%X is a T x P matrix, where T is the time bins, P is the number of
%covariates, the first column of X is not all ones!!!!

% y is a T x 1 matrix

% g is (P) x 1 that specifies group index for each of the covariates  (we
% assume that the first g will always be 1 for now, but we can change
% this up later,  we also assume groups go from 1: a max value;

%link is the link function used.  These modifications will only really be
%used in approximating the loss function w

%c_init is (P+1) x 1 and is the initial set of coefficients or parameters,
%c(1) is c0!!!

M = zeros(P,P);   %not calculated yet, but going to do so soon


js = 1:max(g);

[X, mu_v, std_v]= standardizeX(X);   %X is changing values here!  right now I am using the same variable name to simplify notation.
L_prev = Inf;    %the next estimate of the loss will always be better than infinity!
L = nan(P,1);    %loss function is calculated in each iteration
c = standardizeC(c_init,mu_v, std_v);
c_trace = nan(P+1,max_it+1);
c_trace(:,1)=c;
L_trace = nan(max_it,1);
criteria.df_count = nan(max_it,1);
criteria.BIC_count= nan(max_it,1);
criteria.AIC_count= nan(max_it,1);
criteria.GCV_count= nan(max_it,1);
criteria.df_rat = nan(max_it,1);
criteria.BIC_rat= nan(max_it,1);
criteria.AIC_rat= nan(max_it,1);
criteria.GCV_rat= nan(max_it,1);

penalty = nan(max_it,1);
objective_trace = nan(max_it,1);
N_sparse_groups = nan(max_it,1);
sparse_groups_trace = nan(max_it,max(g));
convergence_max_diff = nan(max_it,1);


%precomputing Xj divisions will speed up code 
Xjs=cell(1,max(g));
XXjs=cell(1,max(g));
for i=1:max(g)
    j_ind = g==i;
    Xjs{i} = X(:,j_ind);
    XXjs{i}= Xjs{i}'*Xjs{i};  %precomputing Xj'*Xj will also speed up computation!
end

X = [ones(size(X,1),1) X];  %this is pre_computed to speed to add a column of ones to X to speed up calculations!

for i_it = 1:max_it
    criteria.df_rat(i_it)=0; %initialize degree of freedom to zero!
    c_old = c;
    w = lossFunc(y,X,c,link);  %w should be a T x 1 vector, in loss function, w is calculated as square root of variance!!
    %update B0
    r = calc_res(X,c,y,link)./w;    %r should be a Tx1 vector    %division in grpreg is made to match that in grpreg code!!!
    c0 = c(1);
    c0 = w'*r/sum(w)+c0;     %this line should do the same as the previous 3, but we need to test this!
    %UPDATE RESIDUALS
    r = r-(c0-c(1))*ones(T,1);   %new minus old
    %Store new coefficient
    c(1) = c0;
    c_trace(1,i_it+1) = c0;
    
    %update Bj
    N_sparse_groups(i_it) = 0;
    for i_j=js
        %disp(['Iteration ' num2str(i_it) ', group ' num2str(i_j)])
        j_ind = find(g==i_j);
        cj = c(j_ind+1);
        [cj, r, dfg] = update_cj(Xjs{i_j},XXjs{i_j},cj,lambda,w,method,group_bridge_gamma,delta_stability,r);
        criteria.df_rat(i_it)=criteria.df_rat(i_it)+dfg;
        c(j_ind+1) = cj;
        c_trace(j_ind+1,i_it+1) = cj;
        
        sparse_groups_trace(i_it,i_j) = (sum(norm(cj))==0);  %if all coefficents are zero, the norm is 0 and it is a sparse group!
    end
    L = calcLoss(X,c,y,link);
    L_trace(i_it)=L;
    if strcmp(method,'groupLasso')
        penalty(i_it) = groupLassoPenalty(c(2:end),g);
        objective_trace(i_it) = L/(2*T)+lambda*penalty(i_it);  %The loss function is scaled by 2 times the number of samples as in the definition in Breheny!!
    elseif strcmp(method,'groupBridge')
        penalty(i_it) = groupBridgePenalty(c(2:end),g,group_bridge_gamma);
        objective_trace(i_it) = L/(2*T)+lambda*penalty(i_it);
    else
        disp('Penalty only defined for groupBridge or groupLasso!')
    end

%     if(  ( L_prev-L_trace(i_it))  < loss_tolerance)   %if improvement in loss function is less than a certain threshold, stop!
%         break
%     end
    L_prev = L;
    
    
    %Calculte the Information Criteria two different ways!!!
    criteria.df_count(i_it) = sum(c~=0);
    df=criteria.df_count(i_it);
    criteria.BIC_count(i_it) = 2*L+df*log(T);
    criteria.AIC_count(i_it) = 2*L+df*2;
    criteria.GCV_count(i_it) = 2*L/(1-(df/T))^2;
    
    df=criteria.df_rat(i_it);
    criteria.BIC_rat(i_it) = 2*L+df*log(T);
    criteria.AIC_rat(i_it) = 2*L+df*2;
    criteria.GCV_rat(i_it) = 2*L/(1-(df/T))^2;
    
  
    [converged, convergence_max_diff(i_it)] = checkConvergence(c_old,c,coeff_tolerance);
    if verbose
    disp(['Iteration ' num2str(i_it) ', Max Convergence Difference is' num2str(convergence_max_diff(i_it))]);
    end
    if converged
        if ~verbose
            disp(['Iteration ' num2str(i_it) ', Max Convergence Difference is' num2str(convergence_max_diff(i_it))]);
        end
        break;
    end
end
df_count = criteria.df_count(i_it);
criteria.analysis.BIC_count_df_portion=df_count*log(T);
criteria.analysis.AIC_count_df_portion=df_count*2;
df_rat = criteria.df_rat(i_it);
criteria.analysis.BIC_rat_df_portion=df_rat*log(T);
criteria.analysis.AIC_rat_df_portion=df_rat*2;
criteria.analysis.double_L=2*L;

max_iterated = i_it;
c = unStandardizeC(c,mu_v, std_v);
c_trace = unStandardizeC(c_trace,mu_v, std_v);

try
    disp(['Memory Ratio is ' num2str(memRatio) ', at end of groupGLM estimation'])
end

end


%% Functions that are called


function [cj, r,dfg] = update_cj(Xj,XXj,cj,lambda,w,method,group_bridge_gamma,delta_stability,r)
dfg=0;  %initializes group degrees of freedom to zero.        
        Pj = length(cj);
        T = size(Xj,1);
        lasso_thresh = 1/T * norm(Xj'*r+XXj*cj)/sqrt(Pj);          %Breheny and Huang, 2009 eq 14
        if strcmp('groupLasso',method) && (lasso_thresh<=lambda)
            %set group to zero, need to calc how this will affect the
            %residuals
            r = r+Xj*cj;
            cj=zeros(Pj,1);
            
        elseif (cj~=zeros(size(cj))) | strcmp(method,'groupLasso')   %updated so that calc is still done if cj is previously set to zero and it is the groupLasso method!
            for i_p=1:length(cj)
                
                if strcmp(method,'groupBridge')
                    lambda_j= lambda*group_bridge_gamma*Pj.^group_bridge_gamma*norm(cj,1).^(group_bridge_gamma-1);
                    
                elseif strcmp(method,'groupLasso')
                    %lasso_thresh = 1/T * norm(Xj'*r+Xj'*Xj*cj)/sqrt(Pj);          %Breheny and Huang, 2009 eq 14
                    lambda_j = lambda*sqrt(Pj)/(norm(cj)+delta_stability);
                end
                cjk = cj(i_p);
                Xjk = Xj(:,i_p);
                [cjk, dfj] = update_cjk_indR(r,Xjk,w,cjk,lambda_j,method);
                r = r-(cjk-cj(i_p))*Xjk;      %new-old
                cj(i_p) = cjk;
                dfg=dfg+dfj;  %updates degree of freedom for each coefficent

            end
        end   
end

function r = calc_res(X,c,y,link) %calculate residual!
if strcmp(link,'probit')
%r = y - normcdf(c(1) + X*c(2:end));  
r = y - normcdf(X*c);  %this function updated to take in an X matrix that has a column of ones
else
    disp('NEED TO PROGRAM RESIDUAL CALC FOR logit!!!')
end
end

function [converged, max_diff] = checkConvergence(c_old,c,coeff_tolerance)

old_zero_ind = (c_old ==0);
new_zero_ind = (c ==0);
if sum(old_zero_ind~=new_zero_ind)>0  %if the zero locations are not the same, it hasn't converged, it can converge though, if the zero locations are the same on the next iteration!   
    converged = 0;
    max_diff = 1;
else
    c=c(c~=0); %division by zero in the next step means that we can only do the comparison with nonzero elements (the nonzero elements here are the same!)
    c_old = c_old(c_old~=0);  
    max_diff = max(abs((c_old - c)./c_old));
    if max_diff > coeff_tolerance
        converged = 0;
    else
        converged = 1;
    end
end

end







