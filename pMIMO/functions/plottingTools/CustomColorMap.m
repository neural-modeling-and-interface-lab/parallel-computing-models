function CustomColorMap(MaxC, MinC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Example:
%%%% Z = peaks;
%%%% surf(Z);
%%%% CustomColorMap(max(max(Z)),min(min(Z)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% RESOLUTION is adjustable
Resolution = 100;

RangeBlue = linspace(255,0,abs(MinC)/2*Resolution)';
RangeRed = linspace(0,255,abs(MaxC)/2*Resolution)';

if MinC >= 0
    B = [RangeRed zeros(size(RangeRed)) zeros(size(RangeRed)) ; ...
        255*ones(size(RangeRed)) RangeRed  zeros(size(RangeRed))]/255;
elseif MaxC <= 0
    B = [zeros(size(RangeBlue)) zeros(size(RangeBlue)) RangeBlue]/255; ...
else
    B = [RangeBlue RangeBlue 255*ones(size(RangeBlue));...
        zeros(size(RangeBlue)) zeros(size(RangeBlue)) RangeBlue;...
        RangeRed zeros(size(RangeRed)) zeros(size(RangeRed)) ; ...
        255*ones(size(RangeRed)) RangeRed  zeros(size(RangeRed)); ...
        ]/255;
    
%         B = [zeros(size(RangeBlue)) zeros(size(RangeBlue)) RangeBlue;...
%         RangeRed zeros(size(RangeRed)) zeros(size(RangeRed)) ; ...
%         255*ones(size(RangeRed)) RangeRed  zeros(size(RangeRed)); ...
%         ]/255;

end
colormap(B);

% RangeBlue = linspace(255,0,abs(MinC)*Resolution)';
% RangeRed = linspace(0,255,abs(MaxC)/2*Resolution)';
%
% if MinC >= 0;
%     B = [RangeRed zeros(size(RangeRed)) zeros(size(RangeRed)) ; ...
%         255*ones(size(RangeRed)) RangeRed  zeros(size(RangeRed))]/255;
% elseif MaxC <= 0
%     B = [zeros(size(RangeBlue)) zeros(size(RangeBlue)) RangeBlue]/255; ...
% else
%     B = [zeros(size(RangeBlue)) zeros(size(RangeBlue)) RangeBlue;...
%         RangeRed zeros(size(RangeRed)) zeros(size(RangeRed)) ; ...
%         255*ones(size(RangeRed)) RangeRed  zeros(size(RangeRed)); ...
%         ]/255;
% end
% colormap(B);
