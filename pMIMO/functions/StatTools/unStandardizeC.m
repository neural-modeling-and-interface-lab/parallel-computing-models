function c = unStandardizeC(c_standard,mu_v, std_v)
% c_standard is a P x it size matrix, where each row corresponds to a
% parameter estimate and each column corresponds to an estimate at a
% specific iteration

%mu_v and std_v are 1xP vectors calculated when standardizing the input matrix, X
%   
%   author: Brian Robinson, 2014-2016
c = zeros(size(c_standard));
it = size(c_standard,2);      %this refers to the amount of iterations!
c(1,:) = c_standard(1,:) - (mu_v./std_v)*c_standard(2:end,:);
c(2:end,:) = c_standard(2:end,:)./repmat(std_v',1,it);  %need to verify this works for de-standardizing and restandardizing!!!!
