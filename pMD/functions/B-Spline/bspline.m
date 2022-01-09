function b = bspline(dp, i, L)% dp- order, i- knots, L- length
t = linspace(0,1,L+1); % range of data
t = t(2:end);

s = linspace(0,1,i);% knot sequence
nknots = length(s); % number of knots

knots=augknt(s,dp); % knot points
colmat=spcol(knots,dp,brk2knt(t,3));
b=colmat(1:3:end,:); % B-spline basis functions

