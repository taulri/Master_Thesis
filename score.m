function [X,Y,score] = score(inclusion_EEG, gt)
% calculate ROC, quantify how good these threshold variables select
outcome = strings(length(gt),1);
for snip = 1:length(gt)
    if gt(snip) == 0 && inclusion_EEG(snip) == 0
        outcome(snip) = 'TN';
    elseif gt(snip) == 1 && inclusion_EEG(snip) == 1
        outcome(snip) = 'TP';
    elseif gt(snip) == 1 && inclusion_EEG(snip) == 0
        outcome(snip) = 'FN';
    elseif gt(snip) == 0 && inclusion_EEG(snip) == 1
        outcome(snip) = 'FP';
    else
        disp('Error!')
        break
    end
end

TPR = sum(outcome(:) == 'TP')/(sum(outcome(:) == 'TP') + sum(outcome(:) == 'FN')); % Sensitivity
FPR = sum(outcome(:) == 'FP')/(sum(outcome(:) == 'FP') + sum(outcome(:) == 'TN')); % 1 - Specificity
X = FPR;
Y = TPR;
score =(0-X)*(0-X) + (1-Y)*(1-Y);
