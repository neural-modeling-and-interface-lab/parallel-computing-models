function [Results, fileNames] = loadAllFittingResults(folder,varargin)
% This function loads all mat files in a folder and returns them in a cell array
% [Results, fileNames] = loadAllFittingResults(folder,varargin)
% Input Arguments
% folder: which folder to search
% Optional Input Arguments
% wildcard: only return results with this wildcard
% verbose: 0 or 1
% min_bytes: (default is 150)

[wildcard, verbose, min_bytes] = process_options(varargin,'wildcard','','verbose',1,'min_bytes',150);
%get all folder conent names that end in mat
listing = dir([folder filesep '*' wildcard '*.mat'] );
non_empty = [listing.bytes]>min_bytes;
listing = [listing(non_empty)];
fileNames = {listing.name};
if isempty(fileNames)
    Results = [];
end
for i=1:length(fileNames)
    fload = [folder filesep fileNames{i}];
    t0=tic;
    try
        Results(i) = load(fload);
    catch ME
        warning([fload ' Not Loaded, ' getReport(ME)])
    end
    tE = toc(t0);
    if verbose
        disp(['Just loaded ' fload ', took ' num2str(tE) 's'])
    end
end
%load all of them!

