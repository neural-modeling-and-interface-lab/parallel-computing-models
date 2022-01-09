function c = c_sparse2full(c_sp,sparse_g,Bff,Bfb,Qff,Qfb,Optff,Optfb,Nin)
%   
%   author: Brian Robinson, 2014-2016
%Q3fbBand = process_options(varargin,'Q3fbBand',[]);
if isfield(c_sp,'sig')
    c.sig = c_sp.sig;
    if isempty(c.sig)
        c.h=[];
        c.k=[];
        disp('No sigma value, we assume this data is empty!')
        return
    end
else  %this is the case when we are making the full data with the un-normalized coefficients
    c.c_0 = c_sp.c_0;
    if isempty(c.c_0)
        c.h=[];
        c.k=[];
        disp('No sigma value, we assume this data is empty!')
        return
    end
end

gk = createG(Bff,Qff,Nin,Optff);
if ~isempty(gk)
    gM = max(gk);
    g_ff=gk;
    g_fb=[];
    [sparse_g_vector, ~, ~,~]= createG_sparse(gk,sparse_g,g_ff,g_fb);
    for i_g=1:gM
        sparse_g_inds = find(i_g==sparse_g_vector);
        g_inds = find(i_g==gk);
        if isempty(sparse_g_inds)
            c.k(g_inds) = 0;
        else
            c.k(g_inds) = c_sp.k(sparse_g_inds);
        end
    end
    c.k=c.k';
else % Modified by Xiwei Dec.13, 2021 for Qff = 0
    gM = 0;
    c.k = [];
end

    sparse_g_fb = sparse_g-gM;
    sparse_g_fb=sparse_g_fb(sparse_g_fb>0);
    gh = createG(Bfb,Qfb,1,Optfb);
if ~isempty(gh)    
    gMh = max(gh);
    g_ff=[];
    g_fb=gh;
    [sparse_g_fb_vector, ~, ~,~]= createG_sparse(gh,sparse_g_fb,g_ff,g_fb);
    for i_g=1:gMh
        sparse_g_inds = find(i_g==sparse_g_fb_vector);
        g_inds = find(i_g==gh);
        if isempty(sparse_g_inds)
            c.h(g_inds) = 0;
        else
            c.h(g_inds) = c_sp.h(sparse_g_inds);
        end
    end
    c.h=c.h';
else % Modified by Xiwei Dec.13, 2021 for Qfb = 0
    c.h = [];
end


