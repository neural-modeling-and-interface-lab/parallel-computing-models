function [Fs, size] = dir_matlab_files(p)

ms = rdir([p filesep '**' filesep '*.m']);
mats = rdir([p filesep '**' filesep '*.mat']);
Fs = {ms.name mats.name};


size=sum([ms.bytes]) + sum([mats.bytes]);

