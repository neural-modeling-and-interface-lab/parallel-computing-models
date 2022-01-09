function c_standard = standardizeC(c,mu_v, std_v)
%   
%   author: Brian Robinson, 2014-2016
c_standard = zeros(size(c));
it = size(c,2);      %this refers to the amount of iterations!
c_standard(2:end,:) = c(2:end,:).*repmat(std_v',1,it);
c_standard(1,:) = c(1,:) + (mu_v./std_v)*c_standard(2:end,:);


