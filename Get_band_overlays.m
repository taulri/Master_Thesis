clc
clear
close all
%% Path and parameters
s = 'Y:\NCH_FL_Forschungsprojekte\Epilepsy\_Master Students\Ulrich Tanja\Results\';
addpath('\\fs-home\ulrta$\Documents\Master_Biomed\Code\DataSet\Function_Needed_10_11_23');
patient = 'Pat1';                                                           %********** choose patient ***********
nintervals = 28;                                                            %********** give number of intervals ****
filter = {'Ref', 'V1', 'V2'};
snipsize = 512; % length of atoms and snippets

%% Load data
for f = 1:size(filter) % loop over filters
    pathfilt = strcat(s,patient,'\',filter{f});
    cd(pathfilt);
    for int = 1:nintervals % loop over intervals
        if int < 10
            ast = '*0';
        else
            ast = '*';
        end
        Rfile = dir(strcat('R_',filter{f},ast,num2str(int),'.mat'));
        FRfile = dir(strcat('FR_',filter{f},ast,num2str(int),'.mat'));
        Spfile = dir(strcat('Spikes_',filter{f},ast,num2str(int),'.mat'));
        load(Rfile.name);
        load(FRfile.name);
        load(Spfile.name);

        % Compare locations of R(ch(2),loc(3)), FR(ch(2),loc(3)) and LocAndCHSpike (loc, CH)

        % Find all ripples which overlay with a spike 
        Sp_Rtable = overlay(LocAndCHSpike(:,2), LocAndCHSpike(:,1), R{1,2}, R{1,3}, snipsize); % ChA, LocA, ChB, LocB, snipsize    

        % Find all fast ripples which overlay with a spike or a ripple
        Sp_FRtable = overlay(LocAndCHSpike(:,2), LocAndCHSpike(:,1), FR{1,2}, FR{1,3}, snipsize); 
        R_FRtable = overlay(R{1,2}, R{1,3}, FR{1,2}, FR{1,3}, snipsize); 

        % Find overlay of Spike, Ripple and Fast ripple
        index_R_Sp = find(Sp_Rtable);
        Sp_R_FRtable = overlay(R{1,2}(index_R_Sp), R{1,3}(index_R_Sp), FR{1,2}, FR{1,3}, snipsize); 
        Matrix_FR_Sp_R_SpFR = [Sp_FRtable, R_FRtable, Sp_R_FRtable];

        % Save results
        cd(pathfilt);
        savein_R = strcat('R_overlay_interval_',num2str(int),'_',filter{f},'_',patient);
        savein_FR = strcat('FR_overlay_interval_',num2str(int),'_',filter{f},'_',patient);
        save(savein_R, "Sp_Rtable");
        save(savein_FR, "Matrix_FR_Sp_R_SpFR");

    end
end






 