%% Get the number of events for every detector version and the corresponding atoms for the reconstruction 
% We first detect events that overlay temporally in the different bands: R,
% FR, Spike and classifiy the events into different overlay-groups. Then we
% further classify the events whether they acheive the criteria of only
% amplitude, AMP anf/or REC. For all this different detector groups we save
% the event counts per channel and the information from the reconstruction
% atoms.

%%
clc
clear
close all

%% Path and parameters
s = 'Y:\NCH_FL_Forschungsprojekte\Epilepsy\_Master Students\Ulrich Tanja\Results\';
addpath('\\fs-home\ulrta$\Documents\Master_Biomed\Code\DataSet\Function_Needed_23_12_20_version2\');
patient = 'Pat20';                                                           %********** choose patient ***********
nchannel = 14;                                                              %********** give number of intervals ****
nintervals = 29;                                                            %********** give number of intervals ****
filter = {'Ref', 'V1', 'V2'};
snipsize = 512; % length of atoms and snippets

%% Get the right band
for f = 1:size(filter,2) % loop over filters
    pathfilt = strcat(s,patient,'\',filter{f});
    for int = 1:nintervals % loop over intervals
        cd(pathfilt);
        if int < 10
            ast = '*0';
        else
            ast = '*';
        end
        
        % R
        Rfile = dir(strcat('R_',filter{f},ast,num2str(int),'.mat'));
        load(Rfile.name);
        % FR
        FRfile = dir(strcat('FR_',filter{f},ast,num2str(int),'.mat'));
        load(FRfile.name);
        % Sp
        Spfile = dir(strcat('Spikes_',filter{f},ast,num2str(int),'.mat'));
        load(Spfile.name);
        
        % Overlaps
        % Find all ripples which overlay with a spike 
        Sp_Rtable = overlay(LocAndCHSpike(:,2), LocAndCHSpike(:,1), R{1,2}, R{1,3}, snipsize); % ChA, LocA, ChB, LocB, snipsize    
        % Find all fast ripples which overlay with a spike or a ripple
        Sp_FRtable = overlay(LocAndCHSpike(:,2), LocAndCHSpike(:,1), FR{1,2}, FR{1,3}, snipsize); 
        R_FRtable = overlay(R{1,2}, R{1,3}, FR{1,2}, FR{1,3}, snipsize); 
        % Find overlay of Spike, Ripple and Fast ripple
        index_R_Sp = find(Sp_Rtable);
        Sp_R_FRtable = overlay(R{1,2}(index_R_Sp), R{1,3}(index_R_Sp), FR{1,2}, FR{1,3}, snipsize); 

        % SpR
        index_Sp_R = find(Sp_Rtable);
        if isempty(index_Sp_R)
            SpR = {zeros(0,4), [{R{5}(index_Sp_R)}, {R{6}(index_Sp_R)}, {R{7}(index_Sp_R)}, {R{8}(index_Sp_R)}], zeros(0,4)};
        else 
            SpR = {[R{1}(index_Sp_R) R{2}(index_Sp_R) R{3}(index_Sp_R) R{4}(index_Sp_R)], [{R{5}(index_Sp_R)}, {R{6}(index_Sp_R)}, {R{7}(index_Sp_R)}, {R{8}(index_Sp_R)}], [R{9}(index_Sp_R) R{10}(index_Sp_R) R{11}(index_Sp_R) R{12}(index_Sp_R)]};
        end
        % SpFR
        index_Sp_FR = find(Sp_FRtable);
        if isempty(index_Sp_FR)
            SpFR = {zeros(0,4), [{FR{5}(index_Sp_FR)}, {FR{6}(index_Sp_FR)}, {FR{7}(index_Sp_FR)}, {FR{8}(index_Sp_FR)}], zeros(0,4)};
        else
            SpFR = {[FR{1}(index_Sp_FR) FR{2}(index_Sp_FR) FR{3}(index_Sp_FR) FR{4}(index_Sp_FR)],[{FR{5}(index_Sp_FR)}, {FR{6}(index_Sp_FR)}, {FR{7}(index_Sp_FR)}, {FR{8}(index_Sp_FR)}],[FR{9}(index_Sp_FR) FR{10}(index_Sp_FR) FR{11}(index_Sp_FR) FR{12}(index_Sp_FR)]};
        end
        % RFR
        index_R_FR = find(R_FRtable);
        if isempty(index_R_FR)
            RFR = {zeros(0,4), [{FR{5}(index_R_FR)}, {FR{6}(index_R_FR)}, {FR{7}(index_R_FR)}, {FR{8}(index_R_FR)}],zeros(0,4)};
        else
            RFR = {[FR{1}(index_R_FR) FR{2}(index_R_FR) FR{3}(index_R_FR) FR{4}(index_R_FR)],[{FR{5}(index_R_FR)}, {FR{6}(index_R_FR)}, {FR{7}(index_R_FR)}, {FR{8}(index_R_FR)}],[FR{9}(index_R_FR) FR{10}(index_R_FR) FR{11}(index_R_FR) FR{12}(index_R_FR)]};
        end
        % SpRFR
        index_Sp_R_FR = find(Sp_R_FRtable);
        if isempty(index_Sp_R_FR)
            SpRFR = {zeros(0,4), [{FR{5}(index_Sp_R_FR)}, {FR{6}(index_Sp_R_FR)}, {FR{7}(index_Sp_R_FR)}, {FR{8}(index_Sp_R_FR)}], zeros(0,4)};
        else
            SpRFR = {[FR{1}(index_Sp_R_FR) FR{2}(index_Sp_R_FR) FR{3}(index_Sp_R_FR) FR{4}(index_Sp_R_FR)],[{FR{5}(index_Sp_R_FR)}, {FR{6}(index_Sp_R_FR)}, {FR{7}(index_Sp_R_FR)}, {FR{8}(index_Sp_R_FR)}],[FR{9}(index_Sp_R_FR) FR{10}(index_Sp_R_FR) FR{11}(index_Sp_R_FR) FR{12}(index_Sp_R_FR)]};
        end

        % Info needed to run the loop over conditions and criteria
        locations = {LocAndCHSpike(:,1), R{1,3}, FR{1,3}, SpR{1,1}(:,3), SpFR{1,1}(:,3), RFR{1,1}(:,3), SpRFR{1,1}(:,3)};
        channels = {LocAndCHSpike(:,2), R{1,2}, FR{1,2}, SpR{1,1}(:,2), SpFR{1,1}(:,2), RFR{1,1}(:,2), SpRFR{1,1}(:,2)};
        IndexAtoms = {{}, R{1,5}, FR{1,5}, SpR{1,2}{1,1}, SpFR{1,2}{1,1}, RFR{1,2}{1,1}, SpRFR{1,2}{1,1}};
        CoefAtoms = {{}, R{1,6}, FR{1,6}, SpR{1,2}{1,2}, SpFR{1,2}{1,2}, RFR{1,2}{1,2}, SpRFR{1,2}{1,2}};
        INDEX_AMP = {{}, R{1,9}, FR{1,9}, SpR{1,3}(:,1), SpFR{1,3}(:,1), RFR{1,3}(:,1), SpRFR{1,3}(:,1)};
        INDEX_REC = {{}, R{1,11}, FR{1,11}, SpR{1,3}(:,3), SpFR{1,3}(:,3), RFR{1,3}(:,3), SpRFR{1,3}(:,3)};
        condname = {'LocAndCHSpike', 'R', 'FR', 'SpR', 'SpFR', 'RFR', 'SpRFR'};
        critname = {'all', 'AMP', 'REC'};
       
        for crit = 1:length(critname)
            switch crit
                case 1
                    for cond = 1:length(channels)
                        % change to folder for overlay conditiona and
                        % inclusion criteria 
                        path = strcat(pathfilt,'\',condname{cond},'\',critname{crit});
                        if (~exist(path))
                            mkdir(path);
                        end
                        cd(path);

                        ChCount = counts(channels{cond}, nchannel);
                        if cond > 1
                           if ~isempty(IndexAtoms{cond})
                                AtomCoef = atoms(IndexAtoms{cond}, CoefAtoms{cond}, channels{cond}, locations{cond});                         
                           else
                                AtomCoef = [];
                           end
                        save(strcat(num2str(int),'_AtomCoef'),"AtomCoef");   
                        end

                        save(strcat(num2str(int),'_ChCount'),"ChCount");
                                          
                    end
    
                case 2
                    for cond = 2:length(channels) % we start from 2 since case 2 and 3 do not have only spikes 
                        % change to folder for overlay conditiona and
                        % inclusion criteria 
                        path = strcat(pathfilt,'\',condname{cond},'\',critname{crit});
                        if (~exist(path))
                            mkdir(path);
                        end
                        cd(path);

                        index_AMP = find(INDEX_AMP{cond});

                        ChCount = counts(channels{cond}(index_AMP), nchannel);
                        
                        if ~isempty(IndexAtoms{cond}(index_AMP))
                            AtomCoef = atoms(IndexAtoms{cond}(index_AMP), CoefAtoms{cond}(index_AMP), channels{cond}(index_AMP), locations{cond}(index_AMP)); 
                        else
                            AtomCoef = [];
                        end

                        %save both
                        save(strcat(num2str(int),'_ChCount'),"ChCount");
                        save(strcat(num2str(int),'_AtomCoef'),"AtomCoef"); 
                        
                    end

                case 3
                    for cond = 2:length(channels)
                        % change to folder for overlay conditiona and
                        % inclusion criteria 
                        path = strcat(pathfilt,'\',condname{cond},'\',critname{crit});
                        if (~exist(path))
                            mkdir(path);
                        end
                        cd(path);

                        index_REC = find(INDEX_REC{cond});

                        ChCount = counts(channels{cond}(index_REC), nchannel);

                        if ~isempty(IndexAtoms{cond}(index_REC))
                            AtomCoef = atoms(IndexAtoms{cond}(index_REC), CoefAtoms{cond}(index_REC), channels{cond}(index_REC), locations{cond}(index_REC)); 
                        else
                            AtomCoef = [];
                        end

                        %save both
                        save(strcat(num2str(int),'_ChCount'),"ChCount");
                        save(strcat(num2str(int),'_AtomCoef'),"AtomCoef"); 
                        
                    end
             end
         end

    end % interval
end % filter


