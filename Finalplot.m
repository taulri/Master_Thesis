function [] = Finalplot(MatrixCHcount, labels, inc, pathfilt, selected_channel, trust)

ParameterFile = 'Y:\NCH_FL_Forschungsprojekte\Epilepsy\_Master Students\Ulrich Tanja\PersistencePlot\RMorphPara.mat';
load(ParameterFile);
Detector = 'Y:\NCH_FL_Forschungsprojekte\Epilepsy\_Master Students\Ulrich Tanja\PersistencePlot\hfod\';
addpath(Detector);

%% Necessary elements
SummaryMats.RBaseLMat = zeros(size(MatrixCHcount));
SummaryMats.RandFRHFOCountMat = MatrixCHcount;
SummaryMats.DurVec = ones(1,size(MatrixCHcount,2))*5;
SummaryMats.DurMat  = ones(size(MatrixCHcount,1),size(MatrixCHcount,2))*5;
SummaryMats.channelNames = labels';

HFOArea = zeros(size(inc,1),1);
HFOArea(selected_channel) = 1;
maskHFOArea = inc;
precenceVEc = round(sum(inc,2)/size(MatrixCHcount,2)*100);

%% Persistence plot
path = strcat(pathfilt , '\PlotsSummaryEQUAL\');
if ~isdir(path)
    mkdir(path)
end
Para = DetPara;
[~]                    = PostProc.ExportHFOResults(SummaryMats, path, 'AllChanEpoch', [], HFOArea, maskHFOArea, precenceVEc,Para,'Epoch', trust);

%% Histogram plot
figure;
bar(MatrixCHcount/5,'black');
hold on;
bar((MatrixCHcount.*inc)/5,'facecolor', 'red', 'edgecolor', 'none');



% Threshold across intervals 97.5
%prc9705 = prctile(MatrixCHcount,97.5,"all");
% th = prc9705/5;

% HFO area
%inc_975 = MatrixCHcount > prc9705;
%HFO_area = inc_975.* inc;
%HFOArea = selected_channel;

%h = yline(th);
% ylim([0 10]);
%set(h,'LineWidth',2,'Color','red');

xtikkis = 1:size(MatrixCHcount,1);
set(gca,'xtick',[],'FontSize',14,'Position',[0.150210084033613,0.16919959473151,0.726428571428571,0.757826747720369],'TickLength',[0.005,0.001]);
ax = copyobj(gca,gcf); % new axis on the figure
ax(1,1).XTickLabelRotation = 45;
set(ax(1,1),'xtick',xtikkis(find(HFOArea)), 'xticklabel', SummaryMats.channelNames(find(HFOArea)),'XColor', 'red','FontSize',11,'YTickLabel', [],'box','off','TickDir','out','TickLength',[0.005,0.001]);
set(ax(1,1),'Position',[0.150210084033613,0.16919959473151,0.726428571428571,0.757826747720369]);
ax2 = copyobj(gca,gcf);
ax2(1,1).XTickLabelRotation = 45;
set(ax2(1,1),'xtick',xtikkis(~HFOArea),'ytick',[], 'xticklabel', SummaryMats.channelNames(~HFOArea), 'Fontsize', 11,'YTickLabel',[],'box','off','TickDir','out','TickLength',[0.005,0.001]);
set(ax2(1,1),'Position',[0.150210084033613,0.16919959473151,0.726428571428571,0.757826747720369]);

xlabel(ax2(1,1),'Channels','FontSize',14);
ylabel('HFO/min','FontSize',14);
title('HFO rates per channel for every epoch')

hold off;



