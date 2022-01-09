function mask = getMask(intervals,t)
%getMask(intervals,t) applies a mask as specified in the given intervals to to the time
%bins specified in t
%   Inputs:
%       intervals is a 2 column matrix, the first column is the start bins, the
%           second column is the end bins
%
%       t is a T x 1 vector, the time bins that we are interested in masking with the given
%           intervals.  
%   Outputs:
%       mask is a Tx1 boolean vector that is true for entries where the
%       specified intervals include the given time indices
%
%For example, if t=(100:200)' and intervals is [50 150], mask
%would return a logical [true(51,1) false(50,1)] because in that original
%t, the first 51 bins are in the specified interal
%e.g. t(1) = 100, t(51)=150 (these bins are in the specified interval)
%   t(52) = 151, t(101) = 200 (these bins are not in the specified
%   interval)
%   
%   author: Brian Robinson, 2014-2016

nInt = size(intervals,1);

%create base time mask, this creates the mask as applied to sequential
%index bins starting from one
maxBin = max([intervals(:); t]); %chooses max bin size to fit all intervals and all supplied t indices
mask_base_time = false(maxBin,1);
for i=1:nInt;
    indStart = intervals(i,1);
    indEnd = intervals(i,2);
    mask_base_time(indStart:indEnd) = true;
end

%applies interval mask to supplied t indices
mask = mask_base_time(t);