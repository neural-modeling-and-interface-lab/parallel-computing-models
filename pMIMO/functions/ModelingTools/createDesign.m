function [V, g, ff_inds, fb_inds] = createDesign(x,y,Bff,Bfb,Qff,Qfb,Optff,Optfb,varargin)
%creates MISO design matrix given, x, y, basis functions objects and model
%order
%
%Inputs:
%   x, binary TxN input matrix where T is number of bins, N is
%       number of inputs
%   y, binary Tx1 ouput matrix
%   Bff, basis object corresponding to inputs (ff refers to feedforward)
%   Bfb, basis object corresponding to output (fb refers to feedback)
%   Optff, feedforward design matrix options, default is empty
%   Optfb, feedback design matrix options, default is empty
%   Optional name value pair arguments:
%       g_include, if specified it will create a design matrix with only
%       the specified groupo indices
%Outputs:
%   V, design matrix TxnCol
%   g, group indices nColx1, creates an index for which group each design
%       matrix column belongs to.
%   ff_inds, returns the column indices that correspond to the feedforward
%       portion of the deisng matrix
%   fb_inds, returns the column indices that correspond to the feedback
%       portion of the design matrix
%
% author: Brian Robinson, 2014-2016

[g_include] = process_options(varargin,'g_include',-1);  %default value of -1 indicates that we are not changing V for any groups
%Create filtered data for x input spike trains
T=length(y);
V=zeros(T,0);    %V is the design matric used for fitting
g = [];   %each column of the design matrix corresponds to the group indicated by the row in g. 
g_ind=1;
for i=1:size(x,2)
    [V_temp, n_gx] = X2V(x(:,i), Bff, Qff,Optff);
    V = [V V_temp];
    for i_g=1:length(n_gx)  %create group index
        g = [g; ones(n_gx(i_g),1)*g_ind];
        g_ind=g_ind+1;
    end
end
ff_inds = 1:size(V,2);
g_ff=g;

%Create filtered data for output y spike train
y_f = [0; y(1:end-1)];
[V_temp, n_gy] = X2V(y_f, Bfb, Qfb,Optfb);
V = [V V_temp];
for i_g=1:length(n_gy)
    g = [g; ones(n_gy(i_g),1)*g_ind];
    g_ind=g_ind+1;
end
fb_inds = (length(ff_inds)+1):size(V,2);
g_fb=g((length(g_ff)+1):end);

%This whole section will re-recrate the desing matrix with only the group indices indicated by g_include
if isempty(g_include) || sum(g_include ~=-1)
    [g, ff_inds, fb_inds,group_inds_to_keep]= createG_sparse(g,g_include,g_ff,g_fb);
    V = V(:,group_inds_to_keep);
end