function val = isGenSTDPFile(fname)
%checks to see whether the file is from a STDP fit by looking for an x, y,
%or gk variable
if exist(fname,'file')
    listing = whos('-file',fname);
    if ismember('x',{listing.name}) && ismember('y',{listing.name}) && ismember('gk',{listing.name})
        val = 1;
    else
        val = 0;
    end
else
    val=0;
end