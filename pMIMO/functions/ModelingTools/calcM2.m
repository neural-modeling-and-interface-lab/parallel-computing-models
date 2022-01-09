function M2= calcM2(V,L,Opt)
% calculates second order basis function design matrix given first order
% design matrix
% Inputs:
%   V, first order design matrix  (T x L)
%   L, number of columns in first order design matrix
%   Opt, default is empty. If Opt.Q2_noSq = 1, the squared second order design matrix components will not be returned
%Outputs:
%   M2, second order design matrix
%
% author: Dong Song, modified by Brian Robinson 2014-2016 to incorporate
% 2nd order options
if isfield(Opt,'Q2_noSq') %
    Q2_noSq = Opt.Q2_noSq;
else
    Q2_noSq = 0;
end
M2 = zeros(length(V), getNgrCoeffs(2,L,Opt));
ni = 0;
if ~Q2_noSq
    for i = 1:L
        for j = i:L
            ni = ni + 1;
            M2(:,ni) = V(:,i).*V(:,j);
        end
    end
else
    for i = 1:L
        for j = i:L
            if i==j
                continue
            end
            ni = ni + 1;
            M2(:,ni) = V(:,i).*V(:,j);
        end
    end
end
end