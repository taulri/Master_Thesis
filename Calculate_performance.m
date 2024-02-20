%% Calculate the performance (sensitivity and specificity) of the different detectors using the classifications of the single patients
% For every detector version we use the classifications of every patient to
% calculate the performance measures. 



%%
clc
clear
close all

%% Path and parameters
s = 'Y:\NCH_FL_Forschungsprojekte\Epilepsy\_Master Students\Ulrich Tanja\Results\';
addpath('\\fs-home\ulrta$\Documents\Master_Biomed\Code\DataSet\Function_Needed_23_12_20_version2');

patient = 1:20;                                                
filterName = {'Reference', 'V1', 'V2'};
condName = {'Spike', 'Ripple', 'FastRipple', 'RippleAndFastRipple', 'SpikeAndRipple', 'SpikeAndFastRipple', 'SpikeAndRippleAndFastRipple'};
critName = {'All_EoI', 'Criteria1_AMP', 'Criteria2_REC'};
filter = {'Ref', 'V1', 'V2'};
cond = {'LocAndCHSpike','R','FR','RFR','SpR','SpFR','SpRFR'};
crit = {'all', 'AMP', 'REC'};
allsens = [];
allspec = [];
allAcc = [];
allPPV = [];
allNPV = [];
allnInc = [];
allmin = [];
allmax = [];
Filter = strings(1,0);
Overlay = strings(1,0);
Feature_Criteria = strings(1,0);
detversion = strings(1,0);

excludedpatietns = {};
classific_all = {};


for f = 1:size(filter,2) % loop over filters
    for co = 1:size(cond,2) % loop over condition 
        for cr = 1:size(crit,2) % loop over criteria
            classifications = strings(1,length(patient));
            trusts = strings(1,length(patient));
            HFOchannels = cell(1,length(patient));
            rates = cell(1,length(patient));
            for pat = 1:length(patient) % loop over patients
                path = strcat(s,'Pat',num2str(patient(pat)),'\',filter{f},'\',cond{co},'\',crit{cr},'\PlotsSummaryEQUAL');
                if (~exist(path))
                    break;
                end

                cd(path);
                load('Classification.mat');
                load('Rates.mat');
                load('Trust.mat');
                load('Temporal_consistency.mat');
                load('Selected_channel.mat');
                % Get the classification for the patient         
                classifications(pat) = classific;
                trusts(pat) = trust;
                HFOchannels{pat} = selected_channel;
                rates{pat} = MatrixCHcount; % EoI/5min

            end
            
            excludedpatietns{end+1} = trusts;
            classific_all{end+1} = classifications;
            %Calculate performance for this detector version
            nNo = 0;
            nYes = 0;
            nTN = 0;
            nTP = 0;
            nFP = 0;
            nFN = 0;
            HFOarearates = [];
     
            for pa = 1:length(patient)
                if strcmp(trusts(pa),'No')
                    nNo = nNo + 1;
                else 
                    HFOarearates = [HFOarearates, reshape(rates{1,pa}(HFOchannels{pa},:),1,[])]; % EoI/5min
                    nYes = nYes + 1;
                    if strcmp(classifications(pa), 'TN')
                        nTN = nTN + 1;
                    elseif strcmp(classifications(pa), 'TP')
                        nTP = nTP + 1;
                    elseif strcmp(classifications(pa), 'FP')
                        nFP = nFP + 1;
                    elseif strcmp(classifications(pa), 'FN')
                        nFN = nFN + 1;
                    end
                end
            end
            
            sensitivity = nTP/(nTP+nFN);    % TP/(TP+FN) How many of the positive do we detect?
            specificity = nTN/(nTN+nFP);     % TN/(TN+FP) How many of the negative do we detect?
            NPV = nTN/(nTN+nFN); % Negative predicitve value 
            PPV = nTP/(nTP+nFP); % Positive predictive value 
            Accuracy = (nTP+nTN)/nYes; % Accuracy 

            %Get rates of pathological channels
            minRate = min(HFOarearates)/5; %HFO/min
            maxRate = max(HFOarearates)/5; %HFO/min


            
            %Save the value 
            detectorpath = strcat(s,filter{f},'\',cond{co},'\',crit{cr});
            if (~exist(detectorpath))
                    mkdir(detectorpath);
            end
            cd(detectorpath);
            save('Sensitivity', "sensitivity");
            save('Specificity', "specificity");

            Filter = [Filter filterName{f}];
            Overlay = [Overlay condName{co}];
            Feature_Criteria = [Feature_Criteria critName{cr}];
            detversion = [detversion strcat(filter{f},'_',cond{co},'_',crit{cr})];       
            allsens = [allsens, sensitivity];
            allspec = [allspec, specificity];
            allNPV = [allNPV, NPV];
            allPPV = [allPPV, PPV];
            allAcc = [allAcc, Accuracy];
            allnInc = [allnInc, nYes];
            allmin = [allmin, minRate];
            allmax = [allmax, maxRate];
     
        end
    end
end

Performances_detectors = [Filter; Overlay; Feature_Criteria; allspec; allsens; allNPV; allPPV; allAcc; allnInc];

%% Getting the mean values and ranges 

cd(s);
load("Performance_all_detectors_EUQAL95.mat")

performance = str2double(Performances_detectors(5:9,:));
out = performance(:,~isnan(performance(1,:)));


%Specificity
meanspec = mean(out(1,:));
max(out(1,:))
min(out(1,:))

%Sensitivity
meansense = mean(out(2,:));
max(out(2,:))
min(out(2,:))

%NPV
meanNPV = mean(out(3,:));
max(out(3,:))
min(out(3,:))

%PPV
meanPPV = mean(out(4,all(~isnan(out))));
max(out(4,all(~isnan(out))))
min(out(4,all(~isnan(out))))

%Accuracy 
meanAcc = mean(out(5,:));
max(out(5,:))
min(out(5,:))

%% Standard deviation of the performance measures across all conditions 
sens = [];
spec = [];
NPV = [];
PPV = [];
ACC = [];

for i = 1:7
     spec = [spec dataoverlay{1,i}(:,1)];
     sens = [sens dataoverlay{1,i}(:,2)];
     NPV = [NPV dataoverlay{1,i}(:,3)];
     PPV = [PPV dataoverlay{1,i}(:,4)];
     ACC = [ACC dataoverlay{1,i}(:,5)];
end

nanstd(spec*100,0,'all')
nanstd(sens*100,0,'all')
nanstd(NPV*100,0,'all')
nanstd(PPV*100,0,'all')
nanstd(ACC*100,0,'all')



%% Figure

% all
boxplot(100*out','Labels', {'Specificity', 'Sensitivity', 'NPV', 'PPV', 'Accuracy'}); 
title('Performance of the detectos');
ylabel('Performance [%]');

% filter
figure(1);
% ref = 100*out(:,1:19)';
% v1 = 100*out(:,20:38)';
% v2 = 100*out(:,39:57)';
% filterdata = {ref, v1, v2};
load("Data_filter.mat");
c = {'Ref', 'V1', 'V2'};
s = {'Specificity', 'Sensitivity', 'NPV', 'PPV', 'Accuracy'};
boxplotGroup(filterdata, 'primaryLabels', c, 'secondaryLabels', s, 'groupLabelType', 'horizontal');
title('Performance of the detectos for the different filters');
ylabel('Performance [%]');

% overlay
figure(2);
% Sp = [out(:,[1,20,39])'; NaN(6,5)];
% R = out(:,[2:4,21:23,40:42])';
% FR = out(:,[5:7,24:26,43:45])';
% R_FR = out(:,[8:10,27:29,46:48])';
% Sp_R = out(:,[11:13,30:32,49:51])';
% Sp_FR = out(:,[14:16,33:35,52:54])';
% Sp_R_FR = out(:,[17:19,36:38,55:57])';
% dataoverlay = {Sp,R,FR,R_FR,Sp_R,Sp_FR,Sp_R_FR};
load("Data_overlay.mat");
c = {'Sp','R','FR','R_FR','Sp_R','Sp_FR','Sp_R_FR'};
s = {'Specificity', 'Sensitivity', 'NPV', 'PPV', 'Accuracy'};
boxplotGroup(dataoverlay, 'primaryLabels', c, 'secondaryLabels', s, 'groupLabelType', 'vertical');
title('Performance of the detectos for the overlay criteria');
ylabel('Performance [%]');

%feature_criteria 
figure(3);
% all = out(:,[2:3:19, 21:3:38, 40:3:57])';
% crit1AMP = out(:,[3:3:19, 22:3:38, 41:3:57])';
% crit2REC = out(:,[4:3:19, 23:3:38, 42:3:57])';
% dataCrit = {all, crit1AMP, crit2REC};
load("Data_featurecrit.mat");
c = {'All', 'AMP', 'REC'};
s = {'Specificity', 'Sensitivity', 'NPV', 'PPV', 'Accuracy'};
gray=[0.5 0.5 0.5];

% All
scatter(repelem(1,18)',dataCrit{1,1}(:,1),'filled', 'MarkerFaceColor', gray)
hold on;
scatter(repelem(5,18)',dataCrit{1,1}(:,2),'filled', 'MarkerFaceColor', gray)
scatter(repelem(9,18)',dataCrit{1,1}(:,3),'filled', 'MarkerFaceColor', gray)
scatter(repelem(13,18)',dataCrit{1,1}(:,4),'filled', 'MarkerFaceColor', gray)
scatter(repelem(17,18)',dataCrit{1,1}(:,5),'filled', 'MarkerFaceColor', gray)
% AMP
scatter(repelem(2,18)',dataCrit{1,2}(:,1),'filled', 'MarkerFaceColor', gray)
scatter(repelem(6,18)',dataCrit{1,2}(:,2),'filled', 'MarkerFaceColor', gray)
scatter(repelem(10,18)',dataCrit{1,2}(:,3),'filled', 'MarkerFaceColor', gray)
scatter(repelem(14,18)',dataCrit{1,2}(:,4),'filled', 'MarkerFaceColor', gray)
scatter(repelem(18,18)',dataCrit{1,2}(:,5),'filled', 'MarkerFaceColor', gray)
% REC
scatter(repelem(3,18)',dataCrit{1,3}(:,1),'filled', 'MarkerFaceColor', gray)
scatter(repelem(7,18)',dataCrit{1,3}(:,2),'filled', 'MarkerFaceColor', gray)
scatter(repelem(11,18)',dataCrit{1,3}(:,3),'filled', 'MarkerFaceColor', gray)
scatter(repelem(15,18)',dataCrit{1,3}(:,4),'filled', 'MarkerFaceColor', gray)
scatter(repelem(19,18)',dataCrit{1,3}(:,5),'filled', 'MarkerFaceColor', gray)
boxplotGroup(dataCrit, 'primaryLabels', c, 'secondaryLabels', s, 'groupLabelType', 'vertical');

hold off;
title('Performance of the detectos grouped by feature criteria');
ylabel('Performance [%]');

%save
% save("Data_filter", 'filterdata');
% save("Data_overlay", 'dataoverlay');
% save("Data_featurecrit", 'dataCrit');


% 
% classific_all_pat(patex == 'No') = 'NA'
%  index = setdiff(1:63, [2,3,23,24,44,45])

