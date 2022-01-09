function lams = setupLassoLamba(X,y,g,N,min_factor,link)
%% Brian Robinson
%   
%   author: Brian Robinson, 2014-2016
lam_max = maxLassoLambda(X,y,g,link) + 1e-6;
lam_min = min_factor*lam_max;
lams=logspace(log10(lam_min),log10(lam_max),N);
lams = fliplr(lams);  %this makes the lambda values descending for use in the groupLasso path!
end

function max_lam = maxLassoLambda(X,y,g,link)
%standardize X
[X, mu_v, std_v]= standardizeX(X); 

%calculate intercept residuals
r = calc_intercept_only_residuals(y,link);
%for each group, find the maximal lambda value
%loop through each j
js = 1:max(g);
max_lam = 0;
for i_j=js
        
        j_ind = find(g==i_j);
        Xj = X(:,j_ind);
        Pj = size(Xj,2);
        T = size(X,1);
        g_max_lam = 1/T * norm(Xj'*r)/sqrt(Pj);
        if g_max_lam>max_lam
            max_lam = g_max_lam;
        end
        disp(['Group ' num2str(i_j) ', max_lam = ' num2str(g_max_lam)])
    %lasso_thresh = 1/T * norm(Xj'*r)/sqrt(Pj); 
end
end


