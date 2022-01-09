function alignXAxes(ax)
%ax is an array of axes
%this makes it so that all of the x axes have the same width.
%this is important when we have timed data that should align between
%subplots
minXW = 1;
for i=1:length(ax)
    thisPos = get(ax(i),'position');
    xW = thisPos(3);
    if xW<minXW
        minXW = xW;
    end
end

for i=1:length(ax)
    thisPos = get(ax(i),'position');
    thisPos(3)=minXW;
    set(ax(i),'position',thisPos)
end


end