%cellGroup is a class for creating and manipulating collections of cellOb class objects
%   cellGroup can be constructed with a .mat file, .csv table or an array of cellobjects
%
%   given different input arguments, one of several input parseing functions
%   is chosen, (e.g. DNMS_str2CellArray)
%
%   author: Brian Robinson, 2014-2016
classdef cellGroup

    properties
        cellArray %is an array of cellOb class objects
    end
    
    methods
        function obj = cellGroup(varargin)
            if nargin == 1   %valid constructor for one input is another cellArray or a mat file
                if iscellstr(varargin{1})     %this constructor allows a list of mat files to be sent and used to construct the cellGroup
                    cellArrayTemp = [];
                    Fs = varargin{1};
                    for i=1:length(Fs);
                        if isNex(Fs{i})
                            cellArrayTemp = [cellArrayTemp nex_str2CellArray(Fs{i})];
                        elseif  is_Converted_Nev(Fs{i})  %this condition has to be first!!!
                            cellArrayTemp = [cellArrayTemp rodentNEV_str2CellArray(Fs{i})];
                        elseif is_rodent_data(Fs{i})
                            cellArrayTemp = [cellArrayTemp rodent_str2CellArray(Fs{i})];
                        else
                            cellArrayTemp = [cellArrayTemp DNMS_str2CellArray(Fs{i})];
                        end
                    end
                    obj.cellArray=cellArrayTemp;
                elseif ischar(varargin{1})   %this constructor allows a single mat file to be sent and used to construct the cellGroup
                    if isNex(varargin{1})
                        cellArrayTemp = [nex_str2CellArray(varargin{1})];
                    elseif  is_Converted_Nev(varargin{1})  %this condition has to be first!!!
                        cellArrayTemp = rodentNEV_str2CellArray(varargin{1});
                    elseif is_rodent_data(varargin{1})
                        cellArrayTemp = rodent_str2CellArray(varargin{1});
                    else
                        cellArrayTemp = DNMS_str2CellArray(varargin{1});
                    end
                    
                    obj.cellArray=cellArrayTemp;
                else
                    obj.cellArray = varargin{1};
                end
            elseif nargin == 2 %valid constructor for two inputs is string for csv file with cell information and string for data Location, this used for Buzsaki data
                cellTableFile=varargin{1};
                dataDir=varargin{2};
                Cell_Data_All = importdata(cellTableFile,',',1);  %load metadata table that has all cell properties
                labels =  Cell_Data_All.textdata(1,:);            %load column heading from the tables
                for i=1:length(labels)                             %get rid of the first two chracters in label i.e. c.clu -> clu
                    labels{i} = labels{i}(3:end);
                end
              
                for data_row=1:size(Cell_Data_All.data,1)          %loop through all rows and create a new cellobject 
                    vars = Cell_Data_All.textdata(data_row+1,:);
                    vars{16} = Cell_Data_All.data(data_row,1);
                    vars{17} = Cell_Data_All.data(data_row,2);
                    vars{18} = Cell_Data_All.data(data_row,3);
                    vars{19} = Cell_Data_All.data(data_row,4);
                    cellArrayTemp(data_row)=cellOb(vars,labels, dataDir);
                end
                obj.cellArray=cellArrayTemp;
            end
        end
        function filteredCellGroup = filterGroup(obj, paramName, paramVals)
            if ~iscell(paramVals)   %filter to only one potential value!
                    filteredCellGroup = filterGroupInd(obj,paramName,paramVals);
            else   %filter with many potential values by concatenating each cellGroups with each individual value!
                filteredCellArrrayTemp = [];
                for i=1:length(paramVals)
                    tempCellGroup=filterGroupInd(obj,paramName,paramVals{i});
                    filteredCellArrrayTemp = [filteredCellArrrayTemp tempCellGroup.cellArray];
                end
                filteredCellGroup=cellGroup(filteredCellArrrayTemp);
            end

            function filteredIndGroup = filterGroupInd(obj,paramName,paramVal)   %filters group to include only those with matching parameter name or parameter value
                if ~ischar(paramVal)
%                     cellArrayTemp = obj.cellArray([obj.cellArray.(paramName)]==paramVal); % Original
                    cellArrayTemp = obj.cellArray(ismember([obj.cellArray.(paramName)], paramVal)); % Modified by Xiwei for tms mimo
                else
                    cellArrayTemp = obj.cellArray(strcmp(paramVal,{obj.cellArray.(paramName)}));
                end
                filteredIndGroup = cellGroup(cellArrayTemp);
            end

        end
        function filteredCellGroup = filterGroupMin(obj, paramName, minParamVal)
            cellArrayTemp = obj.cellArray([obj.cellArray.(paramName)]>=minParamVal);
            filteredCellGroup = cellGroup(cellArrayTemp);
        end
        function filteredCellGroup = filterGroupMax(obj, paramName, maxParamVal)
            cellArrayTemp = obj.cellArray([obj.cellArray.(paramName)]<=maxParamVal);
            filteredCellGroup = cellGroup(cellArrayTemp);
        end
        function filteredCellGroup = filterOutGroup(obj, paramName, paramVals)
            if ~iscell(paramVals)   %filter to only one potential value!
                    filteredCellGroup = filterOutGroupInd(obj,paramName,paramVals);
            else   %filter with many potential values by concatenating each cellGroups with each individual value!
                filteredCellArrrayTemp = [];
                for i=1:length(paramVals)
                    tempCellGroup=filterOutGroupInd(obj,paramName,paramVals{i});
                    filteredCellArrrayTemp = [filteredCellArrrayTemp tempCellGroup.cellArray];
                end
                filteredCellGroup=cellGroup(filteredCellArrrayTemp);
            end

            function filteredIndGroup = filterOutGroupInd(obj,paramName,paramVal)   %filters group to include only those with matching parameter name or parameter value
                if ~ischar(paramVal)
                    cellArrayTemp = obj.cellArray([obj.cellArray.(paramName)]~=paramVal);
                else
                    cellArrayTemp = obj.cellArray(~strcmp(paramVal,{obj.cellArray.(paramName)}));
                end
                filteredIndGroup = cellGroup(cellArrayTemp);
            end
        
        end
        function [yr, t] = getFilteredOutput(obj,varargin)
            [bin_size,downsample,sigma] = process_options(varargin,'bin_size',1,'downsample',50,'sigma',5000);
            maxT = obj.getMaxT(bin_size);
            for i=1:length(obj.cellArray)
                [y_temp, t] = obj.cellArray(i).getSpikeRate('bin_size',bin_size,'downsample',downsample,'sigma',sigma,'maxT',maxT);
                yr(:,i) = y_temp;
            end
        end
        function maxT = getMaxT(obj,bin_size)
            maxT = 0;
            for i=1:length(obj.cellArray)
                temp = floor(max(obj.cellArray(i).getSpTimes('bin_size',bin_size)));
                if temp>maxT, maxT=temp; end
            end
        end
        function cellIDs = getCellIDs(obj)
             for i=1:length(obj.cellArray)
                if ~isempty(obj.cellArray(i).id)
                    cellIDs(i) = obj.cellArray(i).id;
                else
                    cellIDs(i) = obj.cellArray(i).BRIndex;
                end
            end
        end
        function plotFilteredOutput(obj,varargin)
            [bin_size,downsample,sigma,cscl,sortFreq,clims] = process_options(varargin,'bin_size',1,'downsample',50,'sigma',5000,'cscl',[],'sortFreq',0,'clims',[]);
            %sigma is gaussian standard deviation of filter in bins
            %downsample (if memory used in plotting is too high, increase downsample)
            [yr, t] = obj.getFilteredOutput('bin_size',bin_size,'downsample',downsample,'sigma',sigma);
            if ~isempty(cscl)  %cscl overrides clims input argumnet!
                clims = [0 max(yr(:))*cscl];
            end
            cellIds = [obj.getCellIDs'];
            if sortFreq
                [~, sort_ind] = sort(sum(yr),'descend');
                yr = yr(:,sort_ind);
                cellIds = cellIds(sort_ind);
            end
            if ~isempty(clims)
                h=imagesc(t,1:length(obj.cellArray),yr',clims);
            else
                h=imagesc(t,1:length(obj.cellArray),yr');
            end
            set(gca,'ytick',1:length(obj.cellArray))
            set(gca,'yticklabel',num2str(cellIds))
            colorbar;
        end
        function freqs = getFreqs(obj)
            freqs=[];
            for i=1:length(obj.cellArray)
                freqs(i) = obj.cellArray(i).getFreqQuick;
            end
        end
        function filteredCellGroup = filterMinFreq(obj,minFreq)
            freqs = obj.getFreqs;
            filteredCellArrrayTemp = obj.cellArray(freqs>minFreq);
            filteredCellGroup=cellGroup(filteredCellArrrayTemp);

        end
        function filteredCellGroup = filterMaxFreq(obj,maxFreq)
            freqs = obj.getFreqs;
            filteredCellArrrayTemp = obj.cellArray(freqs<maxFreq);
            filteredCellGroup=cellGroup(filteredCellArrrayTemp);
        end
        function n = getN(obj)
            n=length(obj.cellArray);
        end
        function X = getTiming(obj,varargin)
            bin_size = process_options(varargin,'bin_size',1);
            X = zeros(obj.getMaxT(bin_size),obj.getN);
            for i=1:length(obj.getN)
                t_sp = round(obj.cellArray(i).getSpTimes('bin_size',bin_size));
                X(t_sp,i) = 1;
            end
        end
    end
    
end

