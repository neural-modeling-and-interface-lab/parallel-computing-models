%orthExp basis function class created by orthogonalizing exponential
% functions with different time constants
%   interfaces with basis
%
%   author: Brian Robinson 2014-2016
classdef orthExp

    
    properties
        M
        taus
        N
    end
    
    methods
        function obj = orthExp(M,taus)
            obj.M = M;
            obj.taus = taus;
            obj.N = length(obj.taus);
        end
        function B = calcB(obj)
            B = zeros(obj.M, length(obj.taus));
            for i=1:length(obj.taus)
                tau = obj.taus(i);
                B(:,i) = exp(-(1:obj.M)./tau);
            end
            B = gschmidt(B);          
        end
        function N = getN(obj)
            N = obj.N;
        end
        function M = getM(obj)
            M = obj.M;
        end
        function s = descrip(obj)
            s=['M=' num2str(obj.M) ', taus=' num2str(obj.taus) ', N=' num2str(obj.N) ];
        end
        function ise = eq(ob1,ob2)
            if (ob1.M==ob2.M) && (ob1.taus==ob2.taus) 
                ise = true;
            else
                ise = false;
            end
        end
    end
end

