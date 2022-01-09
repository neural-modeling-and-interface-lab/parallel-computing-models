% cellOb is a class for tracking cell properties and loading spike times
% author: Brian Robinson, 2014-2016
classdef cellOb
    properties
        %these properties are directly loaded from metadata
        name
        topdir
        id
        clu
        ele
        excited
        exciting
        nexciting
        inhibited
        inhibiting
        ninhibiting
        session
        behavior
        region
        familiarity
        duration
        eDist
        RefracRatio
        RefracViol
        wire 
        wire_cell 
        FileBase
        dataDir
        nSpikes
        direction
        label
        animal
        BRChannel
        BRIndex
        BRUnit
        LR
        convertedNEVFormat
        tend
        comment
        %these properties get calculated and stored later...
        excitedBy
        %these properties are called dynamically to load data
        spTimes
        
        % rTMS MIMO properties
        channel
        unit
        part
        
    end
    methods
        function obj = cellOb(cellVars,cellVarsLabels, dataDir)
            for i=1:length(cellVarsLabels)
                %some values are fed in as strings, try to turn them into
                %numerical values
                numConvert = str2double(cellVars{i});
                if isnan(numConvert)
                    paramVal = cellVars{i};
                else
                    paramVal = numConvert;
                end
                obj.(cellVarsLabels{i}) = paramVal;
            end
            obj.dataDir=dataDir;
            if isempty(obj.FileBase)  %The filebase is specified when the generating file is a nex file
                if isempty(obj.ele)  % this comparison lets us know it is a DNMS task!  Human tasks are also processed in the same way
                    obj.FileBase = fullfile(obj.dataDir, [obj.session '.mat']);
                else
                    obj.FileBase = fullfile(obj.dataDir,obj.topdir,obj.session,obj.session);
                end
            end
            
        end
        function peTs = getPeriEventTimes(obj,eventTimes,window)
            %window should be specified in seconds
            %eventTimes is in seconds
            %peTs is 2 columns, the first colum is the relative event time,
            %the second column is which event it corresponds too
            sp = obj.getSpTimes('bin_size',1000); %this gets spike times in seconds
            peTs = [];
            for i=1:length(eventTimes)
                peT = sp-eventTimes(i);
                peT = peT(peT>-window & peT < window);
                peTs_thisEvent = [peT ones(length(peT),1)*i];
                peTs = [peTs; peTs_thisEvent];
            end
        end
        function spTimes = getSpTimes(obj,varargin)   %returns vector of spike times in bins (default is 1 ms)
            bin_size = process_options(varargin,'bin_size',1);  %default bin size is 1 ms
            if isNex(obj.FileBase)  %here, we read the spike times from a nex file
                nexFileData = readNexFile(obj.FileBase);
                thisNeuron = nexFileData.neurons{obj.id};
                spTimes = thisNeuron.timestamps*1000/bin_size;
            elseif isGenSTDPFile(obj.FileBase)
                S = load(obj.FileBase,obj.name);
                spTimes = find(S.(obj.name))/bin_size;  %here, we assume everything is in the format of a 1ms spike train in the data file
                if size(spTimes,1) < size(spTimes,2) %here, we are correcting because y trains in this format are saved with the transpose
                    spTimes = spTimes';
                end
            elseif isempty(obj.ele)  % this comparison lets us know it is a DNMS task! Human tasks are also processed in the same way, generated data can be set like this too
                S = load(obj.FileBase,obj.name);
                spTimes = S.(obj.name)*1000/bin_size;
            else
                [T,G,Map,Par]=LoadCluRes(obj.FileBase,obj.ele,obj.clu);
                spTimes=T/Par.SampleRate*1000/bin_size;
            end
        end
        function spTimes = getSpTimesRangeSeconds(obj,t_min,t_max,varargin)
            [bin_size, oldWRONGWAY] = process_options(varargin,'bin_size',1,'oldWRONGWAY',0);
            spTimes_all = obj.getSpTimes('bin_size',bin_size);
            bin_offset = 0;
            spTimes=[];
            for i=1:length(t_min)
                min_bin = round(t_min(i))*1000/bin_size;
                max_bin = round(t_max(i))*1000/bin_size;
                spTimes_temp = spTimes_all(spTimes_all>min_bin & spTimes_all<max_bin);
                spTimes_temp = spTimes_temp-min_bin;
                spTimes_temp=spTimes_temp+bin_offset;
                if oldWRONGWAY
                    bin_offset = max_bin;  %this is how we were running everything before!!!!! AND IT WAS WRONG!!!!!!
                else
                    bin_offset = (max_bin - min_bin) + bin_offset;
                end
                spTimes=[spTimes; spTimes_temp];
            end
            %             %%%old way of doing it
            %             bin_size = process_options(varargin,'bin_size',1);
            %             spTimes_all = obj.getSpTimes('bin_size',bin_size);
            %             min_bin = round(t_min)*1000/bin_size;
            %             max_bin = round(t_max)*1000/bin_size;
            %             spTimes = spTimes_all(spTimes_all>min_bin & spTimes_all<max_bin);
            %             spTimes = spTimes-min_bin;
            
        end
        function freq = getFreq(obj)
            %this function returns spike times, but is slow because it
            %loads all spike times first
            dur = obj.getDur;
            if ~isempty(dur)
                freq = obj.nSpikes/dur;
            else
                freq = 0;
            end
        end
        function freq = getFreqQuick(obj)
            %if tend is specified in properties, this function quickly
            %calculates the frequency
            freq = obj.nSpikes/obj.tend;
        end
        function dur = getDur(obj)
            %this returns the maximum duration in seconds. This is
            %relatively slow vecause it loads asll spike times
            dur = max(obj.getSpTimes)/1000;
        end
        function spSegs = getSpSegTimes(obj,tMaxs,varargin)
            %tMaxs is in bins (not seconds or ms)
            bin_size = process_options(varargin,'bin_size',1);  %default bin size is 1 ms
            spT = obj.getSpTimes('bin_size',bin_size);
            t_start=1;
            for i=1:length(tMaxs)
                t_end = tMaxs(i);
                spSegs{i}=spT(spT>=t_start && spT<=t_end);
                t_start=t_end+1;
            end
        end
        function plotSpikeTimes(obj,varargin)
            [bin_size, linc] = process_options(varargin,'bin_size',1,'linc','b');
            spTimes = obj.getSpTimes('bin_size',bin_size);
            stem(spTimes/1000*bin_size,ones(size(spTimes)),'marker','none','color',linc);
            xlabel('Time (s)')
        end
        function [r, t] = getSpikeRate(obj,varargin) %bin size is in ms, t is in s!
            [bin_size,downsample,sigma,maxT] = process_options(varargin,'bin_size',1,'downsample',50,'sigma',5000,'maxT',[]);
            %sigma is specified in bins
            spTimes = obj.getSpTimes('bin_size',bin_size);
            spTimes = floor(spTimes);
            if isempty(maxT)
                maxT = max(spTimes);
            end
            if spTimes(1)==0 %if there is a spike time smaller than the first bin, we place it in the next bin for plotting purposes
                spTimes(1) = 1;
            end
            train = zeros(maxT,1);
            train(spTimes) = 1;
            r = gaussFilt(train,sigma)./bin_size*1000; %converts to Hz
            r = r(1:downsample:end);
            num_bins = length(r);
            t = (1:num_bins)*bin_size*downsample/1000;
            
        end
        function h = plotSpikeRate(obj,varargin)
            [bin_size,downsample,sigma,color,linw] = process_options(varargin,'bin_size',1,'downsample',50,'sigma',5000,'color','b','linw',1);
            [r, t] = obj.getSpikeRate('bin_size',bin_size,'downsample',downsample,'sigma',sigma);
            line(t,r,'color',color,'linewidth',linw)
        end
        function spTrain = getSpikeTrain(obj,varargin)
            bin_size = process_options(varargin,'bin_size',1);  %default bin size is 1 ms
            spTimes = obj.getSpTimes('bin_size',bin_size);
            
        end
    end
end