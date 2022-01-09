function b =spl(dp,i)% dp- order, i- knots
t = linspace(0,1,2001); % range of data
t = t(2:end);

s = linspace(0,1,i);% knot sequence
nknots = length(s); % number of knots

knots=augknt(s,dp); % knot points
colmat=spcol(knots,dp,brk2knt(t,3));
b=colmat(1:3:end,:); % B-spline basis functions

