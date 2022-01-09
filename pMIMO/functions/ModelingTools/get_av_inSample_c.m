function c = get_av_inSample_c(cv_fit,varargin)
[fold, tilde] = process_options(varargin, 'fold',[],'tilde', 0);
%calculates average coefficient parameters for many in-sample fits
switch tilde
    case 0
        if isempty(fold)
            c_all = [cv_fit.inSampleFits.c];
        else
            c_all = [cv_fit.inSampleFits(fold).c];
        end
        c.sig=mean([c_all.sig]);
    case 1
        if isempty(fold)
            c_all = [cv_fit.inSampleFits.c_tilde];
        else
            c_all = [cv_fit.inSampleFits(fold).c_tilde];
        end
        c.c_0=mean([c_all.c_0]);
end
        c.h=mean([c_all.h],2);
        c.k=mean([c_all.k],2);
        
end
