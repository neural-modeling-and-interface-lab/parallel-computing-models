function [MA, n_g] = X2V(X, basisOb, ModelOrder,Opt)
%X2V creates a design matrix from a single input signal given a basis
%object and a model order
%
%Inputs:
%   x is binary, (T,1) input spike train vector
%   basisOb is a object of class basis
%   ModelOrder is 1,2, or 3 (1 or 2 is recommended)
%   Opt may contain additional options for how to create the design matrix,
%   (current options that are supported are including non-causal basis
%   functions)
%Outputs:
%   MA is the design matrix that is the input filtered through the basis functions
%   n_g contains the number of basis functions for each group (e.g. if
%   ModelOrder=2 with 3 first order and 6 second order design
%   columns, n_g=[3 6];
%
%   author: Dong Song, modified by Brian Robinson to incorporate basis
%   function objects 2014-2016

if ModelOrder == 0
    MA = [];
    n_g=[];
else
    B = basisOb.calcB();
    L = basisOb.getN();  %this is the number of basis functions
    n_g = L;
    M = basisOb.getM();  %this is the memory length of the basis function
    T = length(X);
    V=zeros(T,L);
    for i=1:L
        conv_temp = conv(X,B(:,i));
        V(:,i) = conv_temp(1:T);
    end
    M1 = V;
    %% Supports additional options for creating design matrix
    includeNonCausal = 0;
    if isfield(Opt,'includeNonCausal')  %here we can specify an option that will additionally include the noncasual version of the basis function
        if Opt.includeNonCausal==1
            includeNonCausal = 1;
            Bnc = [B(end:-1:2,:); B];
            Vnc=zeros(T,L);
            for i=1:L
                conv_temp = conv(X,Bnc(:,i));
                Vnc(:,i) = conv_temp(M:T+M-1);
            end
            M1 = [Vnc M1];
            n_g = 2*L;
            if ModelOrder>1
                warning('noncausal basis functions only currently support 1st order models!')
            end
        end
    end
    if isfield(Opt,'onlyNonCausal')
        if Opt.onlyNonCausal==1
            if includeNonCausal ==1
                error('only one of the two options - includeNonCausal AND onlyNonCausal can be specified ')
            end
            Bnc = B(end:-1:1,:);
            Vnc=zeros(T,L);
            for i=1:L
                conv_temp = conv(X,Bnc(:,i));
                Vnc(:,i) = conv_temp(M:T+M-1);
            end
            M1 = Vnc;
        end
    end
    
end
    
if ModelOrder==1
    MA = M1;
     
elseif ModelOrder==2
    M2= calcM2(V,L,Opt);
    MA = horzcat(M1, M2);
    n_g = [getNgrCoeffs(1,L,Opt) getNgrCoeffs(2,L,Opt)];
    
elseif ModelOrder==3
    M2= calcM2(V,L,Opt);
    [M3 n_g]= calcM3(X,V,L,Opt,basisOb);
    MA = horzcat(M1, M2, M3);
    
else
    if ModelOrder~=0
        error('Invalid Model Order');
    end
end
end


function [M3, n_g] = calcM3(X,V,L,Opt,basisOb)
    if isfield(Opt,'Q3')
       Q3 = Opt.Q3;
    else
        Q3 = [];
    end

    if ~isempty(Q3)
        Bands = Q3.calcQ3(X);
        M3 = [];
        n_g = [getNgrCoeffs(1,L,Opt) getNgrCoeffs(2,L,Opt) ];
        for i=1:Q3.Nbands
            M_temp=X2V(Bands(:,i), basisOb, 1);
            M3 = [M3 M_temp];
            n_g = [n_g size(M_temp,2)];
        end
    else
        M3 = zeros(length(V), L*(L+1)*(L+2)/6);
        ni = 0;
        for i = 1:L
            for j = i:L
                for k = j:L
                    ni = ni + 1;
                    M3(:,ni) = (V(:,i).*V(:,j)).*V(:,k);
                end
            end
        end
        M3 = calcM3(X,V,L);
        n_g = [getNgrCoeffs(1,L,Opt) getNgrCoeffs(2,L,Opt) L*(L+1)*(L+2)/6];
    end


end