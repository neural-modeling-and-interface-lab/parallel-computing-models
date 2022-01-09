function meta_exist = existMetaData(fol,stageLabel,varargin)
saveSubFolder = process_options(varargin,'saveSubFolder', fullfile('exampleData','exampleResult',fol));
f = getMetaDataFile(fol,'saveSubFolder',saveSubFolder);
listing = whos('-file',f);
var_names = {listing.name};
meta_exist  = sum(ismember(var_names,stageLabel));