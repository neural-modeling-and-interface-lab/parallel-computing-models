%% Brian Robinson
%% Draft Version 4.30.15

function [cjk, df] = update_cjk_indR(r,Xjk,w,cjk,lambda_j,method)
%   
%   author: Brian Robinson, 2014-2016

%delta_stability is a small term for making sure that the lambda in the groupLasso calculation doesn't reach inf when cj drops to zero (as recommended in Breheny). we can try a few terms to find an ideal value
%Pj = size(Xj,2);
T = size(Xjk,1);
% Old way of calculating updated cjk when it wa not an individual value!
% Zj = Xjk'.*w';
% num = 1/T*(Zj*r+diag(Zj*Xjk).*cjk);
% den = 1/T*diag(Zj*Xjk);

xw=Xjk'.*w'; %(1 x T)
xwr = xw*r;  %(1 x T) * (T x 1) = 1 x 1
xwx =  xw*Xjk;  % (1xT)*(Tx1) = 1x1
num=1/T*(xwr+xwx*cjk);
den=1/T*xwx;

if strcmp(method,'groupMCP') || strcmp(method,'groupBridge')
    %     if strcmp(method,'groupMCP')
    %         %lambda_j = lambda;
    %         disp('groupMCP not programmed yet!!!!!')
    %     elseif strcmp(method,'groupBridge')
    %         lambda_j= lambda*group_bridge_gamma*Pj.^group_bridge_gamma*norm(cj,1).^(group_bridge_gamma-1);
    %     end
    cjk = soft_threshold(num,lambda_j)./den;
    
elseif strcmp(method,'groupLasso')
    
    cjk = num./(den+lambda_j);
    
end

df = abs(cjk)/abs(xwr/xwx+cjk);  %this line calculates df by comparison to the unpenlized coefficent estimate!

end