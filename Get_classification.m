%% Get Classification for every detector versoin and patient
% We calculate the prediction for the channels for every patient taking into account the temporal consistency across 
% 5-minute intervals and compare the predicted channels to the resected
% channels and the seizure outome. This allows the classification of every
% patientes for every detector version into TP, TN, FP, FN. 

%%
clc
clear
close all

%% Path and parameters
s = 'Y:\NCH_FL_Forschungsprojekte\Epilepsy\_Master Students\Ulrich Tanja\Results\';
addpath('\\fs-home\ulrta$\Documents\Master_Biomed\Code\DataSet\Function_Needed_23_12_20_version2');
patient = 20;                                                              %********** choose patient ***********
nchannel = 14;                                                             %********** give number of channels ****
nintervals = 29;                                                           %********** give number of intervals ****
filter = {'Ref', 'V1', 'V2'};
cond = {'FR', 'LocAndCHSpike', 'R', 'RFR', 'SpFR', 'SpR', 'SpRFR'};
crit = {'all', 'AMP', 'REC'};
resectedChannels = [3,4,5,11,12];                                          %********** resected channels ******
seizure_outcome = 5; % 1 = freedom, 1< = recurrent seizures                %********** seizure outcome ********


% Excluded channels
exChannels = [1];                                                          %********** excluded channels **********
nrows = 1:nchannel;
Chkeep = setdiff(nrows,exChannels);


% Patient data label for plot
patdata = 'Y:\NCH_FL_Forschungsprojekte\Epilepsy\_Vasileios\Copied from Hard Drive\Vasileios\HFO Analysis\Fedele Data\Data';
cd(patdata);
if patient < 10
    load(strcat('pat_0',num2str(patient),'_night_01_interval_01.mat'));
else
    load(strcat('pat_',num2str(patient),'_night_01_interval_01.mat'));
                    end
labels = data.label_bipolar(1:nchannel);


labels = labels(Chkeep,:);


%% Get the right band
for f = 1:size(filter,2) % loop over filters
    for co = 1:size(cond,2) % loop over condition
        for cr = 1:size(crit,2) % loop over criteria
            pathfilt = strcat(s,'Pat',num2str(patient),'\',filter{f},'\',cond{co},'\',crit{cr});
            if (~exist(pathfilt))
                break;
            end
            cd(pathfilt);
            dirnames = dir('*_ChCount.mat');
            for n = 1:length(dirnames)
                Chcount{n} = load(dirnames(n).name);
            end
            MatrixCHcount = struct2array(cell2mat(Chcount));
            intcol = 2:2:2*length(dirnames);
            MatrixCHcount = MatrixCHcount(:,intcol);


            % Remove excluded channels 
            MatrixCHcount = MatrixCHcount(Chkeep,:);


            % Channel with counts above or equal 95th percentile for every interval 
            prc95 = prctile(MatrixCHcount, 95);
            for i = 1:size(MatrixCHcount,2)
                inc(:,i) = MatrixCHcount(:,i) >= prc95(i);
                if prc95(i) == 0
                    inc(:,i) = 0; % If it is 0, we do not want to include the channels 
                end
            end
            

            % Temporal consistancy: Which channels are in the 97.5
            % percentile for temporal consistancy?
            countabovethperchannel = sum(inc,2);  
            prc975 = prctile(countabovethperchannel, 97.5);
            selected_channel = find(countabovethperchannel >= prc975);

            trust = 'Yes';
            absolutetemcons = countabovethperchannel/size(MatrixCHcount,2)*100; % We do not trust if the median of the selected channels is below 50%
            if median(absolutetemcons(selected_channel)) <= 50
                trust = 'No';
            end

            %Plot 
            if ~(nnz(MatrixCHcount) == 0) && ~(nnz(inc) == 0)
                Finalplot(MatrixCHcount, labels, inc, pathfilt, selected_channel, trust);
            end


            % Classification: TP, TN, FP, FN
            if seizure_outcome == 1 
                    if any(~ismember(selected_channel,resectedChannels))
                        classific = 'FP';
                    else
                        classific = 'TN';
                    end 
            else
                    if any(~ismember(selected_channel,resectedChannels))
                        classific = 'TP';
                    else
                        classific = 'FN';
                    end              
            end 

            
            % Save results
            path = strcat(pathfilt , '\PlotsSummaryEQUAL\');
            if ~isdir(path)
                mkdir(path)
            end
            cd(strcat(pathfilt, '\PlotsSummaryEQUAL\'));
            % Channel
            save('Selected_channel',"selected_channel");
            % Temporal consistancy
            tempcons = absolutetemcons(selected_channel);
            save('Temporal_consistency',"tempcons");
            % Classification
            save('Classification', "classific");
            % Rates
            save('Rates', 'MatrixCHcount');
            % Trust
            save('Trust', 'trust');

        end %crit
    end %cond
end %filt




