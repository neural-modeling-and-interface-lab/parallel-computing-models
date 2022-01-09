function P = groupLassoPenalty(c,g)
%This function calculates the groupLasso penalty term
%it is assumed that c and g are the same length, and that g is all consecutive values starting from one (e.g. 1 1 1 2 2 3 3 3 3 4 4 4)
%unpenalized terms (the offset), should not be passed in the c vector!!!!!
%   
%   author: Brian Robinson, 2014-2016

P=0;
js = 1:max(g);
for i_j = js
    cj = c(g==i_j);
    K = length(cj);     %In Song 2009, this group weighting term is called ai, and removed since all groups have the same number of coefficients.  This term is included in Breheny and Huang 2009, however, so we include it here..
    P=P+sqrt(K)*norm(cj);   
end