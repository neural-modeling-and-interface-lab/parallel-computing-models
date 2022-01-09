function xfilt = gaussFilt(x,sigmas)
%filters a spike train by convolving a Gaussian
%x can have dimensions T x 1 or T x N
%sigma is the standard deviation of the gaussian that is convolved with x
%multiple sigma values can be provided to perform multiple filtering
%operations
if length(sigmas)>1  %create matrix with many sigma values
    xfilt = nan(length(x),length(sigmas));
    for i=1:length(sigmas)
        sigma = sigmas(i);
        xfilt(:,i) = filtTrain(x,sigma);
    end
else   %create matrix that filters many concatenated input trains.
    sigma = sigmas;
    xfilt = nan(size(x));
    for i=1:size(x,2);
        xfilt(:,i) = filtTrain(x(:,i),sigma);
    end
end
end

function gaussFilt = filtTrain(x,sigma)
size = sigma*6;
xf = linspace(-size / 2, size / 2, size);
gaussFilter = exp(-xf .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); % normalize
gaussFilt = conv(x, gaussFilter, 'same');
end

