function n = getNgrCoeffs(Q,L,Opt)
% returns number of columns in a design matrix for a model order given the number of first order columns
% Inputs: 
    %Q, model order
    %L, number of first order design columns
% Outputs:
    %n number of Q order design columns
if isfield(Opt,'Q2_noSq')
   Q2_noSq = Opt.Q2_noSq;
else
   Q2_noSq = 0;
end

switch Q
    case 0
        n=0;
    case 1
    n=L;
    case 2
        if Q2_noSq
            n=L*(L+1)/2-L;
        else
            n = L*(L+1)/2;
        end
    case 3
        disp('Calc in parent function for now!!')

end