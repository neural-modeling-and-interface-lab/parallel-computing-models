function c = mcc(m)
% MCC Matthews correlation coefficient
% m: confusion matrix
% c: Matthews correlation coefficient

if size(m,1) == 1 && size(m, 2) == 1
    c = 0;
else
    TP = m(1,1); % True Positive
    TN = m(2,2); % True Negative
    FP = m(1,2); % False Positive
    FN = m(2,1); % False Negative
    
    if ((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))== 0
        c = (TP*TN-FP*FN);
    else
        c = (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    end
end
