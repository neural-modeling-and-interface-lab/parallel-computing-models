% immediateB defines a basis function class with an impulse response
%   interfaces with basis
%
%   author: Brian Robinson, 2014-2016
classdef immediateB

    
    properties
        M
        N
    end
    
    methods
        function obj = immediateB(M)
            obj.M = M;
            obj.N = 1;
        end
        function B = calcB(obj)
            B = zeros(obj.M,obj.N);
            B(1)=1;

        end
        function N = getN(obj)
            N = obj.N;
        end
        function M = getM(obj)
            M = obj.M;
        end
        function s = descrip(obj)
            s=[' '];
        end
        
        function ise = eq(ob1,ob2)
            
            if (ob1.M==ob2.M) 
                ise = true;
            else
                ise = false;
            end
        end
        
    end
    
end

