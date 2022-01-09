function w = lossFunc(y,X,c,link)
%calculates for logit or probit
%   
%   author: Brian Robinson, 2014-2016
%X = [ones(size(X,1),1) X];
switch link
    case 'logit',
        logit_theta = X*c;
        theta=exp(logit_theta)./(1+exp(logit_theta));
        %sqtR=sqrt(theta.*(1-theta));
        w = sqrt(theta.*(1-theta));  %this matches the notation in Breheny!
    case 'probit',
        probit_theta = X*c;
        theta = normcdf(probit_theta);
        part1 = -y.*normpdf(probit_theta).*probit_theta./theta./(1-theta);
        part2 = -y.*(normpdf(probit_theta).^2).*(1-2*theta)./theta./(1-theta)./theta./(1-theta);
        part3 = probit_theta.*normpdf(probit_theta)./(1-theta);
        part4 = -(normpdf(probit_theta).^2)./(1-theta)./(1-theta);
        %sqtR = sqrt(-part1-part2-part3-part4);
        w = sqrt(-part1-part2-part3-part4); %this matches the notation in Breheny!
    otherwise
        disp('no such link function');
end

w((theta<.001)|(theta>.999)) = .001;  %this condition is found in Breheny 2009 code