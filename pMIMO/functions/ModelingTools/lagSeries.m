%lagSeries defines a laguerre basis functions object calculated in
%series
%   can be interfaced with basis
%
%   author: Brian Robinson 2014-2016
classdef lagSeries

    
    properties
        M % memory length in bins
        L % order
        a % alpha value
    end
    
    methods
        function obj = lagSeries(M,L,a)
            obj.M = M;
            obj.L = L;
            obj.a = a;
        end
        function B = calcB(obj)
            B = zeros(obj.M, obj.L);
            for i = 1:obj.L
                B(:,i) = lagu_re(i, obj.M, obj.a);
            end       
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

