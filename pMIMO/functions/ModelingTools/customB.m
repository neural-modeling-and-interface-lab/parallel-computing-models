%customB defines a custom single basis function
%   Often used in comparisons to test a specific single waveform
%
%   author: Brian Robinson, 2014-2016
classdef customB

    
    properties
        B
        M
        N
    end
    
    methods
        function obj = customB(B)
            obj.B = B;
            obj.M = length(B);
            if isempty(B)   %this allows us to make a dummy custom function with an empty B
            obj.N=0;
            else
            obj.N = 1;
            end
        end
        function B = calcB(obj)
            B = obj.B;
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
            
            if (ob1.B==ob2.B)  
                ise = true;
            else
                ise = false;
            end
        end
        
    end
    
end

