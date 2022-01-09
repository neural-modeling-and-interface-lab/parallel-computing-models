function [g_sp, ff_inds, fb_inds,group_inds_to_keep]= createG_sparse(g,g_include,g_ff,g_fb)
%   
%   author: Brian Robinson, 2014-2016
n_ff_vars=0;
n_fb_vars=0;
group_inds_to_keep=[];
for i=g_include'
    group_inds_to_keep=[ group_inds_to_keep find(i==g')];  %loops through each group to keep and adds their corresponding indices
    n_ff_vars=n_ff_vars+length(find(i==g_ff));
    n_fb_vars=n_fb_vars+length(find(i==g_fb));
end
g_sp=g(group_inds_to_keep);
ff_inds = 1:n_ff_vars;
fb_inds = (n_ff_vars+1):(n_ff_vars+n_fb_vars);
