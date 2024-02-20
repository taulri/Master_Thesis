clc
clear
close all
%% Path and parameters
s = 'Y:\NCH_FL_Forschungsprojekte\Epilepsy\_Master Students\Ulrich Tanja\Results\';
addpath('\\fs-home\ulrta$\Documents\Master_Biomed\Code\DataSet\Function_Needed_10_11_23');
patient = 'Pat1';                                                           %********** choose patient ***********
nintervals = 28;                                                            %********** give number of intervals ****
filter = {'Ref', 'V1', 'V2'};
%band = {'R', 'FR', 'Spike'}; 
snipsize = 512; % length of atoms and snippets

%% Get the right band
for f = 1:size(filter) % loop over filters
    pathfilt = strcat(s,patient,'\',filter{f});
    cd(pathfilt);
    for int = 1:nintervals % loop over intervals
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
        
        % Overlays 
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
        SpR = { [R{1}(index_Sp_R) R{2}(index_Sp_R) R{3}(index_Sp_R) R{4}(index_Sp_R)], [{R{5}(index_Sp_R)}, {R{6}(index_Sp_R)}, {R{7}(index_Sp_R)}, {R{8}(index_Sp_R)}], [R{9}(index_Sp_R) R{10}(index_Sp_R) R{11}(index_Sp_R) R{12}(index_Sp_R)]};
        % SpFR
        index_Sp_FR = find(Sp_FRtable);
        SpFR = {[FR{1}(index_Sp_FR) FR{2}(index_Sp_FR) FR{3}(index_Sp_FR) FR{4}(index_Sp_FR)],[{FR{5}(index_Sp_FR)}, {FR{6}(index_Sp_FR)}, {FR{7}(index_Sp_FR)}, {FR{8}(index_Sp_FR)}],[FR{9}(index_Sp_FR) FR{10}(index_Sp_FR) FR{11}(index_Sp_FR) FR{12}(index_Sp_FR)]};
        % RFR 
        index_R_FR = find(R_FRtable);
        RFR = {[FR{1}(index_R_FR) FR{2}(index_R_FR) FR{3}(index_R_FR) FR{4}(index_R_FR)],[{FR{5}(index_R_FR)}, {FR{6}(index_R_FR)}, {FR{7}(index_R_FR)}, {FR{8}(index_R_FR)}],[FR{9}(index_R_FR) FR{10}(index_R_FR) FR{11}(index_R_FR) FR{12}(index_R_FR)]};
        % SpRFR
        index_Sp_R_FR = find(Sp_R_FRtable);
        SpRFR = {[FR{1}(index_Sp_R_FR) FR{2}(index_Sp_R_FR) FR{3}(index_Sp_R_FR) FR{4}(index_Sp_R_FR)],[{FR{5}(index_Sp_R_FR)}, {FR{6}(index_Sp_R_FR)}, {FR{7}(index_Sp_R_FR)}, {FR{8}(index_Sp_R_FR)}],[FR{9}(index_Sp_R_FR) FR{10}(index_Sp_R_FR) FR{11}(index_Sp_R_FR) FR{12}(index_Sp_R_FR)]};


        conditions = {LocAndCHSpike R FR SpR SpFR RFR SpRFR};
        condname = {'LocAndCHSpike', 'R', 'FR', 'SpR', 'SpFR', 'RFR', 'SpRFR'};
        critname = {'all', 'AMP', 'REC'};
        
       
        for crit = 1:length(critname)
            switch crit
                case 1
                    for cond = 1:length(conditions)
                        counts(condname(cond), critname(crit), int, conditions{cond}(:,2));
                        if cond > 1
                            atoms(condname(cond), critname(crit),conditions{cond}(:,5), int, condition{cond}(:,6)); 
                        end
                    end
    
                case 2
                    for cond = 2:length(conditions)
                        index_AMP = conditions{cond}(9);
                        counts(condname(cond), critname(crit), int, conditions{cond}(index_AMP,2));
                        atoms(condname(cond), critname(crit),condition{cond}(index_AMP,5), int, condition{cond}(index_AMP,6)); 
                    end

                case 3
                    for cond = 2:length(conditions)
                        index_REC = conditions{cond}(11);
                        counts(condname(cond), critname(crit), int, conditions{cond}(index_AMP,2));
                        atoms(condname(cond), critname(crit),condition{cond}(index_AMP,5), int, condition{cond}(index_AMP,6)); 
                    end
             end
         end

    end % interval
end % filter


