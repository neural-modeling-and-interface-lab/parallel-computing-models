function [Zs, b, b95, KS_score, Tau]= discretetimerescaling(Pb, y)
%Applies time rescaling and KS test
% Inputs:
%   Pb, Tx1 estimated probability per bin
%   y, Tx1 observations per bin
% Ouputs:
%   Zs, sorted transformation to uniform distribution
%   b, cdf of reference uniform distribution
%   b95, cdf 95 confidence interval
%   KS_score 
%   Tau, rescaled intervals
%
%   author: Dong Song
T_time = find(y); % spike location;


q= -log(1-Pb);
Tau = 0; % in case there's no interval
for i = 2:length(T_time)
    r=rand(1,1);
    Tau(i-1) = sum( q((T_time(i-1)+1):(T_time(i)-1)) ) - log(1-r*(1-exp(-q(T_time(i))))); % Rescaled Tau
end

Z = 1 - exp(-Tau); % Transform to Uniform Distribution
Zs = sort(Z);

n = length(Zs);
b = 1:1:n; b=b';
b = (b-.5)/n;
b95 = horzcat((b-1.36/sqrt(n)), (b+1.36/sqrt(n)));
KS_score = max(abs(Zs'-b)) /  max(abs(b95(:,1)-b));

