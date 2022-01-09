function c = get_av_inSample_c_tilde(cv_fit)
%   
%   author: Brian Robinson, 2014-2016
c_all = [cv_fit.inSampleFits.c_tilde];
c.h=mean([c_all.h],2);
c.k=mean([c_all.k],2);
c.c_0=mean([c_all.c_0]);

end
