%lagSeriesZS defines a laguerre basis function object that is zero sum
%   commonly used for fitting feedback functions to limit postive
%   feedback
%   Can be interfaced with basis
%
%   author: Brian Robinson 2014-2016
classdef lagSeriesZS

    
    properties
        M
        L
        a
    end
    
    methods
        function obj = lagSeriesZS(M,L,a)
            obj.M = M;
            obj.L = L;
            obj.a = a;
        end
        function B = calcB(obj)
            B = zeros(obj.M, obj.L+1);
            for i = 1:(obj.L+1)
                B(:,i) = lagu_re(i, obj.M, obj.a);
            end
            B_new = B(:,2:end) - repmat(B(:,1),1,obj.L); %subtracts the first basis function from the remaing basis functions
            B = gschmidt(B_new); %re-orthogonalizes
        end
        function N = getN(obj)
            N = obj.L;
        end
        function M = getM(obj)
            M = obj.M;
        end
        function s = descrip(obj)
            s=['M=' num2str(obj.M) ', L=' num2str(obj.L) ', a=' num2str(obj.a) ];
        end
        function ise = eq(ob1,ob2)
            if (ob1.L==ob2.L) && (ob1.M==ob2.M) &&(ob1.a==ob2.a) 
                ise = true;
            else
                ise = false;
            end
        end
    end
end

