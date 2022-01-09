function f = getMetaDataFile(fol,varargin)
saveSubFolder = process_options(varargin,'saveSubFolder', fullfile('exampleData','exampleResult',fol));
f = fullfile(saveSubFolder,[fol 'metaData.mat']);