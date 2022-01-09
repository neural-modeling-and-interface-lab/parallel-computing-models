function C = c_tilde2C_stim(c_tilde, sigma)
%   
%   author: Brian Robinson, 2014-2016
C0 = c_tilde.c_0*sigma+1;
Ck = c_tilde.k*sigma;
Ch = c_tilde.h*sigma;
C = [Ck; C0; Ch];