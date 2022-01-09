%basis is a modular wrapper class for any type basis function object
%   basis hasa common calcB, getN, getM, and descrip method that can be
%   used with several basis function class types
%   
%   author: Brian Robinson, 2014-2016
classdef basis
    
    properties
        type
        basisOb
    end
    
    methods
        function obj = basis(type, varargin)
            obj.type=type;
            switch type
                case 'customB'
                    B=varargin{1};
                    obj.basisOb = customB(B);
                case 'immediateB'
                    M = varargin{1};
                   obj.basisOb = immediateB(M);
                case 'orthExp'
                    M = varargin{1};
                    taus = varargin{2};
                    obj.basisOb = orthExp(M,taus);
                case 'lagSeries'
                    M = varargin{1};
                    L = varargin{2};
                    a = varargin{3};
                    obj.basisOb = lagSeries(M,L,a);
                case 'lagSeriesZS'
                    M = varargin{1};
                    L = varargin{2};
                    a = varargin{3};
                    obj.basisOb = lagSeriesZS(M,L,a);
                case 'lag' % deprecated laguerre class calculation
                    M = varargin{1};
                    L = varargin{2};
                    a = varargin{3};
                    if nargin == 5
                        fr_0 = varargin{4};
                    else
                        fr_0 = 0;
                    end
                    obj.basisOb = lagB(M,L,a,fr_0);
                case 'bspl'
                    dp = varargin{1};
                    i = varargin{2};
                    L = varargin{3};
                    p = varargin{4};
                    if nargin ==6
                        rmlb=varargin{5};
                    else
                        rmlb = 1;
                    end
                    obj.basisOb = splB(dp, i, L, p,rmlb);
                case 'pieceB'
                    Ls = varargin{1};
                    obj.basisOb = pieceB(Ls);
                otherwise
                    disp('Not valid basis type choice, only lag and bspl supported') 
            end
        end
        function B = calcB(obj)
            B=obj.basisOb.calcB();
        end
        function N = getN(obj)
            N=obj.basisOb.getN();
        end
        function M = getM(obj)
            M=obj.basisOb.getM();
        end
        function s= descrip(obj)
            s=obj.basisOb.descrip();
        end
        function ise = eq(ob1,ob2)
            if strcmp(ob1.type,ob2.type)
                ise = (ob1.basisOb==ob2.basisOb);
            else
                ise=false;
            end
        end
        
    end

end

