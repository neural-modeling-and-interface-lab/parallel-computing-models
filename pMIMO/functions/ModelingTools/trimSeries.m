function [y_trim, t_bins] = trimSeries(y, includeIntervals,removeIntervals)
%trims a binned series based on according to includeIntervals and
%removeIntervals
% 
% Inputs:
%   y, T x M series to be trimmed
%   includeIntervals, 2 x NIntervals matrix, each row represents an interval,
%       the 1st column is the start index, the second column is the end
%       index. The include intervals are relative to the binned y
%       series. If empty, all indices will be included
%   removeIntervals, 2 x NIntervals matrix formatted as above.  All of
%       these intervals wil be trimmed out after the includeIntervals
%       are created. 
%
% Outputs:
%   y_trim, Ttrim X 1 series after trimming (Ttrim<=T)
%   t_bins, time bin indices of y_trim relative to the original y
%       e.g. y(t_bins)=y_trim
%   
%   author: Brian Robinson, 2014-2016
t_bins = (1:size(y,1))';
if isempty(includeIntervals)
    includeIntervals = [1 max(t_bins)];
end
includeMask = getMask(includeIntervals,t_bins);
removeMask = getMask(removeIntervals,t_bins);
y_trim = y(includeMask & ~removeMask,:);
t_bins=t_bins(includeMask & ~removeMask);
