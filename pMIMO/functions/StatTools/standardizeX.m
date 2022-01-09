function [X0s, mu_v, std_v]= standardizeX(X)
%standardizeX standardizes a design matrix
%we standardize one column at a time in order to use less memory!
%   
%   author: Brian Robinson, 2014-2016
Ti = size(X,1);
N = size(X,2);
X0s=nan(size(X));
mu_v=nan(1,N);
std_v=nan(1,N);

for i=1:N
    X_temp = X(:,i);
    mu_v_temp = mean(X_temp);
    if norm(X_temp)==0 %if there is an empty column, we return a std dev of 1 so that it essentially just does not get modified
        X0s(:,i)=X_temp;
        mu_v(i)=mu_v_temp;
        std_v(i) = 1;
    else
        mu_M = ones(Ti,1)*mu_v_temp;
        X0 = X_temp-mu_M;
        std_v_temp= std(X0) * sqrt((Ti-1)/Ti);  %matlab's std function normalizes by N-1 instead of N, this is a correction factor
        std_M = ones(Ti,1)*std_v_temp;
        X0s_temp = X0./std_M;
        X0s(:,i)=X0s_temp;
        mu_v(i)=mu_v_temp;
        std_v(i) = std_v_temp;
    end
end

end