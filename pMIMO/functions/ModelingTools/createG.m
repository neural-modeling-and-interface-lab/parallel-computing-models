function g = createG(B,Q,Nin,Opt)
%   
%   author: Brian Robinson, 2014-2016
if isfield(Opt,'Q3')
    Q3=Opt.Q3;
else
    Q3=[];
end
g_ind=1;
g=[];
L=B.getN;
if Q == 0 % Modified by Xiwei Dec.13,2021
    n_gx = 0;
elseif Q==1
    n_gx = L;
elseif Q == 2
    n_gx = [getNgrCoeffs(1,L,Opt) getNgrCoeffs(2,L,Opt)];
elseif Q ==3
    n_gx = [getNgrCoeffs(1,L,Opt) getNgrCoeffs(2,L,Opt) ones(1,Q3.Nbands)*L];
end
for i_in=1:Nin
    for i_g=1:length(n_gx)  %create group index
        g = [g; ones(n_gx(i_g),1)*g_ind];
        g_ind=g_ind+1;
    end
end
end