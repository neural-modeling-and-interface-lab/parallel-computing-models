function data=nev2DataMat(NEV,outputName,CA3_indices,CA1_indices,R_channels,L_channels)
%this function prepares spiking data by merging NEV data with channel
%region and R vs. L

%first, use NEV data to make equivalent matrix
TS_All = [double(NEV.Data.Spikes.Electrode'), double(NEV.Data.Spikes.Unit'), (double(NEV.Data.Spikes.TimeStamp)./double(NEV.MetaTags.SampleRes))'];  %converts NEV into TS format (see get Timestamps section)
TS_InputIndices=TS_All(:,1)*10+TS_All(:,2);  %find the neuron indices that correspond to all time stamps

%add CA3 data
for i=1:length(CA3_indices)
region = 'CA3';
ind = CA3_indices(i);
un=mod(ind,10);
ch=floor(ind/10);
if ismember(ch,R_channels)
LR='R';
else
    LR='L';
end
name = [region '_n' sprintf('%04d',ind) '_ch' sprintf('%03d',ch) '_un' sprintf('%01d',un) '_' LR];
TS_indices = find(ismember(TS_InputIndices,ind)); %find the indices of the time stamps that correspond to the specified input neurons.
TS = TS_All(TS_indices,3);
data.(name)=TS;
end

%add CA1 data
%same exact as above except with CA1
for i=1:length(CA1_indices)
region = 'CA1';
ind = CA1_indices(i);
un=mod(ind,10);
ch=floor(ind/10);
if ismember(ch,R_channels)
    LR='R';
else
    LR='L';
end
name = [region '_n' sprintf('%04d',ind) '_ch' sprintf('%03d',ch) '_un' sprintf('%01d',un) '_' LR];
TS_indices = find(ismember(TS_InputIndices,ind)); %find the indices of the time stamps that correspond to the specified input neurons.
TS = TS_All(TS_indices,3);
data.(name)=TS;
end

%add additional NEV data to file
data.MetaTags = NEV.MetaTags;
data.IO_events=NEV.Data.SerialDigitalIO.UnparsedData;
data.IO_TS = NEV.Data.SerialDigitalIO.TimeStampSec;
data.source_indices=CA3_indices;
data.OutputIndices=CA1_indices;

%save file
save(outputName,'-struct','data');