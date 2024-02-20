%% Statistical testing
% With a 3-way ANOVA we test the influence of the three cathegorical
% variables: Filters (3 levels: Ref, V1, V2), overlap (7 levels: Sp, R, FR,
% SpR, SpFR, RFR, SPRFR) and criteria (3 levels: All, REF, AMP)

clc
clear
close all

%% Path and parameters
s = 'Y:\NCH_FL_Forschungsprojekte\Epilepsy\_Master Students\Ulrich Tanja\Results';
addpath('\\fs-home\ulrta$\Documents\Master_Biomed\Code\DataSet\Function_Needed_23_12_20_version2\');

cd(s);
load('Performance_all_detectors_EUQAL95.mat');
filters = repelem(["Ref" "V1" "V2"], [21 21 21]);
overlay = repmat(repelem(["LocAndCHSpike" "R" "FR" "RFR" "SpFR" "SpR" "SpRFR"], [3 3 3 3 3 3 3]), 1,3);
criteria = repmat(["all" "AMP" "REC"], 1, 21);

%% Test assumptions
performance = str2double(Performances_detectors(5:8,:));
indspikedet = [2,3,23,24,44,45];
variables = Performances_detectors(2:4,:);
variables(:,indspikedet) = [];
spec = Performances_detectors(5,:);
spec(indspikedet) = [];
sens = Performances_detectors(6,:);
sens(indspikedet) = [];


%% Linear mixed model 
% Sensitivity
tblsense = table(str2double(sens)', variables(3,:)', variables(2,:)', 'VariableNames', {'Sensitivity', 'FeatureCrit', 'Overlap'});
lmesense = fitlme(tblsense, 'Sensitivity~1+FeatureCrit+(1|Overlap)');

% Specificity 
tblspec = table(str2double(spec)', variables(3,:)', variables(2,:)', 'VariableNames', {'Specificity', 'FeatureCrit', 'Overlap'});
lmespec = fitlme(tblspec, 'Specificity~1+FeatureCrit+(1|Overlap)');



