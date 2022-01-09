function k3plot = pk3_multiB(B,t,c3v,varargin)
Q2_noSq = process_options(varargin,'Q2_noSq',0);
% K3PLOT calculate the k3 matrix
%   k3plot = pk2(M, a, t, c3v), where
%   M is memory length, a is alpha, t 
%   is [min max bin], c3v is 3rd order
%   LVK coefficient vector

if ~Q2_noSq
    L = floor(sqrt(2*length(c3v)));
else
    L = ceil(sqrt(2*length(c3v)));
end

% prepare for the lagureer base functions
b=B.calcB;
% for k=1:L
% lagu(1:M,k)=lagu_re(k,M,a);
% end

% resample
ti = roundti(t);
n = 1;
for i=ti(1):ti(3):ti(2)
    db(n,:) = b(i,:);
    n = n + 1;
end

if Q2_noSq
    c3v = cNosq2cFull(c3v,L);
end

[c3half,c3full] = v2m(c3v);

% make the 3rd order kernel matrix
k3plot = zeros(n-1,n-1);
for i=1:L
    for j=1:L
        vv = c3full(i,j)*db(:,i)*(db(:,j)');
        k3plot = k3plot + vv;
    end
end
