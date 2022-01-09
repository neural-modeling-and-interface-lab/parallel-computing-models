function loss = calcLoss(X,c,y,link)
%in the probit source code for loss it has something saying that y
%needs to be either 1 or -1, not 1 or 0.  Not sure if this is also true for the logit loss.
switch link
    case 'probit'
        loss = ProbitLoss(c,X,y);
    case 'logit'
        loss = LLoss(c,X,y);
end
end