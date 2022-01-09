function r = calc_intercept_only_residuals(y,link)
%   
%   author: Brian Robinson, 2014-2016
%this code was made by modifying the groupGLM intercept algorithm to only
%include the very first iteration when c is initiliaized to zero!!
%  c = 0;
% for i = 1:100
%     T = length(y);
%     X = ones(T,0);
%    
%     w = lossFunc(y,X,c,link);  %w should be a T x 1 vector, in loss function, w is calculated as square root of variance!!
%     if strcmp(link,'probit')
%         %r = (y-.5)./w;    %w = -w;
%         r = calc_res(X,c,y,link)./w;
%     else
%         disp('NEED TO PROGRAM RESIDUAL CALC FOR logit!!!')
%     end
%     c0=c;
%     c0 = w'*r/sum(w)+c0;     %this line should do the same as the previous 3, but we need to test this!
%     %r = calc_res(X,c,y,link)./w;
%     r = r-(c0-c(1))*ones(T,1);   %new minus old
%     r2 = calc_res(X,c,y,link)./w;
%     c = c0;
%  
% end

%calc residual using glmfit!
T = length(y);
X = ones(T,0);
X1 = ones(T,1);
c = glmfit(X,y,'binomial',link);
w = lossFunc(y,X1,c,link);
r = calc_res(X,c,y,link)./w;
end



function r = calc_res(X,c,y,link) %calculate residual!
if strcmp(link,'probit')
    r = y - normcdf(c(1) );
else
    disp('NEED TO PROGRAM RESIDUAL CALC FOR logit!!!')
end
end

%
%     %update B0
%
%
%     r = calc_res(X,c,y,link)./w;    %r should be a Tx1 vector    %division in grpreg is made to match that in grpreg code!!!
%
%     function r = calc_res(X,c,y,link) %calculate residual!
% if strcmp(link,'probit')
% r = y - normcdf(c(1) + X*c(2:end));
% else
%     disp('NEED TO PROGRAM RESIDUAL CALC FOR logit!!!')
% end
% end
%
%
%     c0 = c(1);
% %     x0 = ones(T,1);            %x0 should be Tx1;
% %     z0 = (x0').*(w');     % z should be a 1 x T row vector!
% %     c0 = z0*r/z0*x0+c0;
%     c0 = w'*r/sum(w)+c0;     %this line should do the same as the previous 3, but we need to test this!
%     r = r-(c0-c(1))*ones(T,1);   %new minus old