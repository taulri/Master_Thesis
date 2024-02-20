%% Initial event detection and OMP reconstruction 
% Detect spikes, ripples and fast ripples using an initial amplitude detector. 
% Then reconstruct the ripples and fast ripples using and OMP algorithm and
% save the reconstruction information. Select events based on two differents
% critera using (1) the reconstruction percentage and (2) the amplitue
% distribution of events before and after the first reconstruction. 

% Input data: 1x1 struct with 15 fields:
%        x (nchannelx600000 double)
%        label (1xnchannel cell)
%        patient_number
%        night
%        interval
%        nch 
%        sampling frequency
%        Rec_timestamp
%        BipChOrder
%        R
%        FR
%        FRandR
%        Lag_seconds
%        label_bipolar
%        Rec_datetime 



%% 
clc
clear
close all

%% Add path to functions
s = 'F:\Code_complete';
path_to_Data = strcat('Y:\NCH_FL_Forschungsprojekte\Epilepsy\_Vasileios\Copied from Hard Drive\Vasileios\HFO Analysis\Fedele Data\Data');
path_to_function = s;
addpath(path_to_Data);
addpath(path_to_function);
Base = path_to_function;
List = dir(fullfile(Base, '*.mat*'));

%% Load data
cd(path_to_Data);

%for pats = 1:20 
pats = 1;                   %*********** chose patient *****************
if pats < 10
    name = 'pat_0';
else
    name = 'pat_';
end
pat = strcat(name,num2str(pats),'_night_0');  

intnumb = 1;
nightnumb = 1;

while size(dir(strcat(pat,num2str(nightnumb),'*')),1)

for ints = 1:size(dir(strcat(pat,num2str(nightnumb),'*')),1)

    if intnumb < 10
        filenamesave = strcat(pat,num2str(nightnumb),'_interval_0',num2str(intnumb));          
    else
        filenamesave = strcat(pat,num2str(nightnumb),'_interval_',num2str(intnumb));          
    end
    
load(filenamesave)

bipol_data = DB_HFO_create_bipolar(data);

%% Required variables 
Nc = 8; % minimum number of crossings
%NcS = 2; % minimum number of crossings spike 
fs = 2000; % sampling frequency 
size_atom = 512; % length of atoms and snippets
% channels = [24,24,24,24,34,65,37,65,45,34,37,37,25,42,14];       % *********** to change for every patient ***************
% channeln = channels(pats-5);
channeln = 23;
min_freq = [80, 250]; % min frequency of R and FR
filters = {'Ref', 'V1', 'V2'};
sigmaREC = 2;
sigmaAMP = 3;

%% Filter loop
for type = 1:3
    filttype = filters{type};
                          
%% Patient path
path_to_pat = (strcat('Y:\NCH_FL_Forschungsprojekte\Epilepsy\_Master Students\Ulrich Tanja\Results\', 'Pat',num2str(pats),'\',filttype,'\')); % *********** to change for every patient ***************
if (~exist(path_to_pat))
    mkdir(path_to_pat);
end
addpath(path_to_pat);

%% Filtering HFO bands
cd(Base)
[Rfiltered, FRfiltered] = filter_HFOBands(bipol_data.x(1:channeln,:), fs); % RfilteredData_G RfilteredData_40_10 RfilteredData_60_20
switch filttype
    case 'Ref'
    filtdata = {Rfiltered{1,1} FRfiltered{1,1}}; 
    case 'V1'
    filtdata = {Rfiltered{1,2} FRfiltered{1,2}}; 
    case 'V2'
    filtdata = {Rfiltered{1,3} FRfiltered{1,3}};   
end

%% HFO snippets
[snipR, snipFR, LocAndCHR, LocAndCHFR, threshold] = snippets_HFO(filtdata, min_freq, fs, Nc, size_atom);
th_ripple = threshold(1:channeln,:);
th_fastripple = threshold(channeln+1:channeln*2,:);

% channelR = channel_HFO(1:size(snipR,1));
% channelFR = channel_HFO(size(snipR,1)+1:end);
% thresholdR = threshold_HFO(1:channeln,:);
% thresholdFR = threshold_HFO(channeln+1:end,:);
% locationR = location_HFO(1:size(snipR,1));
% locationFR = location_HFO(size(snipR,1)+1:end);

%% Spike snippets
[snipSpike_raw, snipSpike_R, snipSpike_FR, channel_Sp, location_Sp] = snippets_spike_detector(bipol_data.x(1:channeln,:),filtdata{1}, filtdata{2}, size_atom, fs);
LocAndCHSpike = [location_Sp, channel_Sp];

%% Plot snippets 
% % Plot spikes
% snippetlist_Sp = {channel_Sp, location_Sp};
% plot_snippets(snippetlist_Sp, bipol_data.x, [], 'r', ' Spike'); % r because we want 10 random snippets plotted, other option: all
% 
% % Plot R
% snippetlistR = {LocAndCHR(:,2), LocAndCHR(:,1)};
% plot_snippets(snippetlistR, filtdata{1,1} , [], 'r', ' Ripple');
% 
% Plot FR
% snippetlistFR = {LocAndCHFR(:,2), LocAndCHFR(:,1)};
% plot_snippets(snippetlistFR, filtdata{1,2} , [], 'r', ' Fast ripple', th_fastripple);
 
% % Plot RFR
% snippetlistRFR = {channel_HFO, location_HFO};
% plot_snippets(snippetlistRFR, filtdata{1,1}, filtdata{1,2}, 'r', ' Ripple and Fast ripple');
% 
% % Plot spikes and R
% plot_snippets(snippetlist_Sp, bipol_data.x, filtdata{1,1}, 'r', 'Spike and Ripple');

%% Create the dictionary
N = size_atom; % Length of dictionary
[Dictlocinfo,DL,DR,DF] = Create_Dictionary(N,2000); % create the Gabor dictionary at 2000 Hz sampling frequency
Dict.DL = DL.Atom; % High
Dict.DR = DR.Atom; % Medium 
Dict.DF = DF.Atom; % Low
Dict.frq.DL = DL.frq;
Dict.frq.DR = DR.frq;
Dict.frq.DF = DF.frq;
Locatoms = [Dictlocinfo.Aloc, Dictlocinfo.Bloc, Dictlocinfo.Cloc, Dictlocinfo.Dloc, Dictlocinfo.Eloc, Dictlocinfo.Floc, Dictlocinfo.Gloc, Dictlocinfo.Hloc,Dictlocinfo.Iloc,Dictlocinfo.Jloc,Dictlocinfo.Kloc,Dictlocinfo.Lloc];

% xlim([0,512]);
% %ylim([0,1]);
% title("Slow wave atom"); 
% plot(Dict.DR(:,8700))
% xlim([0,512]);
% %ylim([0,1]);
% title("Ripple atom");
% plot(Dict.DF(:,8000))
% xlim([0,512]);
% %ylim([0,1]);
% title("Fast ripple atom");


%% OMP for R and FR snippets
%Ripple
OMP_R_list = {};
V_R_list = [];
for i = 1:size(LocAndCHR,1) 
    [OMP_R, V_R] = OMP_reconst_draw(snipR(i,:)',Dict,fs,Nc,80, th_ripple(LocAndCHR(i,1),LocAndCHR(i,2)-(size_atom/2):LocAndCHR(i,2)+(size_atom/2)-1)); % OMP process with Gabor Dictionary for one snippet
    OMP_R_list = [OMP_R_list; OMP_R];
    V_R_list = [V_R_list; V_R]; 
end
percRec_R = []; %Percentage of reconstruction after the first atom 
for i = 1:size(LocAndCHR,1) 
    percRec_R(i) = 1 - OMP_R_list{i,1}.Error(1);
end

%R table results
Rcoef = {};
Rlocatom = {};
Rerror = {};
RV = {};
RnRec = [];
for i = 1:length(OMP_R_list)
    Rerror{i} = 1-OMP_R_list{i}.Error; % Percentage of reconstruction for each snippet, if more than one interation we have a vector 
    Rlocatom{i} = OMP_R_list{i}.loc; % Chosen atom 
    indexcoef = find(OMP_R_list{i}.coeff); % Index of the coefficient 
    Rcoef{i} = OMP_R_list{i}.coeff(indexcoef); % Coefficient 
    indexV = find(V_R_list(i,:));
    RV{i} = V_R_list(i, indexV); % V-factor 
    RnRec(i) = length(Rerror{i}); % Number of reconstructions 
end


%Fast ripple
OMP_FR_list = {};
V_FR_list = [];
for i = 1:size(LocAndCHFR,1) 
    [OMP_FR, V_FR] = OMP_reconst_draw(snipFR(i,:)',Dict,fs,Nc,250, th_fastripple(LocAndCHFR(i,1),LocAndCHFR(i,2)-(size_atom/2):LocAndCHFR(i,2)+(size_atom/2)-1)); % OMP process with Gabor Dictionary for one snippet
    OMP_FR_list = [OMP_FR_list; OMP_FR];
    V_FR_list = [V_FR_list; V_FR];  
end
percRec_FR = [];
for i = 1:size(LocAndCHFR,1) 
    percRec_FR(i) = 1 - OMP_FR_list{i,1}.Error(1);
end

%FR table results
FRcoef = {};
FRlocatom = {};
FRerror = {};
FRV = {};
FRnRec = [];
for i = 1:length(OMP_FR_list)
    FRerror{i} = 1-OMP_FR_list{i}.Error;
    FRlocatom{i} = OMP_FR_list{i}.loc;
    indexcoef = find(OMP_FR_list{i}.coeff);
    FRcoef{i} = OMP_FR_list{i}.coeff(indexcoef);
    indexV = find(V_FR_list(i,:));
    FRV{i} = V_FR_list(i, indexV);
    FRnRec(i) = length(FRerror{i});
end

%% Getting n_r random locations where there is no HFO 
n_r = 500;
% Ripple
[randomlocR] = randomnonHFOloc(n_r, LocAndCHR, size_atom, size(bipol_data.x,2), channeln);
%plot_snippets({randomlocR(:,2), randomlocR(:,1)}, filtdata{1,1} , [], 'all', ' Ripple');
randsnipR = cutsnippets(filtdata{1,1}, n_r, size_atom, randomlocR(:,1), randomlocR(:,2));
% Fast Ripple
[randomlocFR] = randomnonHFOloc(n_r, LocAndCHFR, size_atom, size(bipol_data.x,2), channeln);
%plot_snippets({randomlocFR(:,2), randomlocFR(:,1)}, filtdata{1,2} , [], 'all', ' Fast Ripple');
randsnipFR = cutsnippets(filtdata{1,2}, n_r, size_atom, randomlocFR(:,1), randomlocFR(:,2));

%% OMP for these random snippets 
%Ripple
Rand_OMP_R_list = {};
Rand_V_R_list = [];
for i = 1:size(randomlocR,1) 
    [rand_OMP_R, rand_V_R] = OMP_reconst_draw(randsnipR(i,:)',Dict,fs,Nc,80, th_ripple(randomlocR(i,2),randomlocR(i,1)-(size_atom/2):randomlocR(i,1)+(size_atom/2)-1)); % OMP process with Gabor Dictionary for one snippet
    Rand_OMP_R_list = [Rand_OMP_R_list; rand_OMP_R];
    Rand_V_R_list = [Rand_V_R_list; rand_V_R]; 
end
%Fast ripple
Rand_OMP_FR_list = {};  
Rand_V_FR_list = [];
for i = 1:size(randomlocFR,1) 
    [rand_OMP_FR, rand_V_FR] = OMP_reconst_draw(randsnipFR(i,:)',Dict,fs,Nc,250, th_fastripple(randomlocFR(i,2),randomlocFR(i,1)-(size_atom/2):randomlocFR(i,1)+(size_atom/2)-1)); % OMP process with Gabor Dictionary for one snippet
    Rand_OMP_FR_list = [Rand_OMP_FR_list; rand_OMP_FR];
    Rand_V_FR_list = [Rand_V_FR_list; rand_V_FR];  
end

%% Calculating the averege reconstrction error after one atom for random snippets
%Ripples
percRec_R_rand = [];
for i = 1:length(Rand_OMP_R_list)
    percRec_R_rand(i) = 1 - Rand_OMP_R_list{i}.Error(1);
end
%Fast ripples
percRec_FR_rand = [];
for i = 1:length(Rand_OMP_FR_list)
    percRec_FR_rand(i) = 1 - Rand_OMP_FR_list{i}.Error(1);
end

%% Selecting snippets which need only one rec. to go below threshold 

%Select the HFO snippets for which we have one and only one resonstruction 

%Ripple
if isempty(V_R_list)
    index_selected_R = [];
else
    index_selected_R = setdiff(find(V_R_list(:,2)),find(V_R_list(:,3)));
end
OMP_R_list_selected = {OMP_R_list{index_selected_R}};
snipR_selected = snipR(index_selected_R,:);
%Rec. error for snippets with only one reconstruction 
percRec_R_sel = [];
for i = 1:length(OMP_R_list_selected)
    percRec_R_sel(i) = 1 - OMP_R_list_selected{1,i}.Error(1);
end

%Fast Ripple 
if isempty(V_FR_list)
    index_selected_FR = [];
else
    index_selected_FR = setdiff(find(V_FR_list(:,2)),find(V_FR_list(:,3)));
end
OMP_FR_list_selected = {OMP_FR_list{index_selected_FR}};
snipFR_selected = snipFR(index_selected_FR,:);
%Rec. error for snippets with only one reconstruction 
percRec_FR_sel = [];
for i = 1:length(OMP_FR_list_selected)
    percRec_FR_sel(i) = 1 - OMP_FR_list_selected{1,i}.Error(1);
end

%% Variation criteria using the distribution of amplitudes for snippets with only one reconstruction

% Ripple
[range_std_R, inc_R] = sigmarange(snipR_selected, sigmaAMP, OMP_R_list_selected);
index_R = find(inc_R);

inc_AMP_R = zeros(size(LocAndCHR,1),1);
inc_AMP_R(index_selected_R(index_R)) = 1;
th_AMP_R = zeros(size(LocAndCHR,1),1);
th_AMP_R(find(inc_AMP_R == 1)) = range_std_R(index_R,2);

inc_snip_R = {OMP_R_list_selected{index_R}};

%Rec. error for snippets with only one reconstruction + amp. dist. criteria
percRec_R_sel_amp = [];
for i = 1:length(inc_snip_R)
    percRec_R_sel_amp(i) = 1 - inc_snip_R{1,i}.Error(1);
end


% Fast Ripple
[range_std_FR, inc_FR] = sigmarange(snipFR_selected, sigmaAMP, OMP_FR_list_selected);
index_FR = find(inc_FR);

inc_AMP_FR = zeros(size(LocAndCHFR,1),1);
inc_AMP_FR(index_selected_FR(index_FR)) = 1;
th_AMP_FR = zeros(size(LocAndCHFR,1),1);
th_AMP_FR(find(inc_AMP_FR == 1)) = range_std_FR(index_FR,2);

inc_snip_FR = {OMP_FR_list_selected{index_FR}};

%Rec. error for snippets with only one reconstruction + amp. dist. criteria
percRec_FR_sel_amp = [];
for i = 1:length(inc_snip_FR)
    percRec_FR_sel_amp(i) = 1 - inc_snip_FR{1,i}.Error(1);
end

%% Reconstrction error criteria for snippets with only one reconstruction

%Ripples
th_rec_error_R = mean(percRec_R_rand)+sigmaREC*std(percRec_R_rand);
inc_rec_r = zeros(length(OMP_R_list_selected),1);
for i = 1:length(OMP_R_list_selected)
    inc_rec_r(i) = (1-OMP_R_list_selected{1,i}.Error(1)) > th_rec_error_R;
end
index_R_rec = find(inc_rec_r);

inc_REC_R = zeros(size(LocAndCHR,1),1);
inc_REC_R(index_selected_R(index_R_rec)) = 1;
th_REC_R = zeros(size(LocAndCHR,1),1);
th_REC_R(find(inc_REC_R == 1)) = th_rec_error_R;


%Fast ripples
th_rec_error_FR = mean(percRec_FR_rand)+sigmaREC*std(percRec_FR_rand);
inc_rec_fr = zeros(length(OMP_FR_list_selected),1);
for i = 1:length(OMP_FR_list_selected)
    inc_rec_fr(i) = (1-OMP_FR_list_selected{1,i}.Error(1)) > th_rec_error_FR;
end
index_FR_rec = find(inc_rec_fr);

inc_REC_FR = zeros(size(LocAndCHFR,1),1);
inc_REC_FR(index_selected_FR(index_FR_rec)) = 1;
th_REC_FR = zeros(size(LocAndCHFR,1),1);
th_REC_FR(find(inc_REC_FR == 1)) = th_rec_error_FR;

%% Plot Hist of dist of rec error for different pools: Random snippets, all RandFR, all RandFR after one rec., all RandFR after one rec. + amp variation criteria

% Ripples
figure(1);
histogram(percRec_R_rand*100, 20);
hold on
histogram(percRec_R_sel*100, 20);
%histogram(percRec_R, 20);
%histogram(percRec_R_sel_amp, 20);
xline(mean(percRec_R_rand*100)+sigmaREC*std(percRec_R_rand*100), 'r', strcat(num2str(sigmaREC) ," \times \sigma"));
xline(mean(percRec_R_rand*100), 'r', '\mu random snippets');
hold off
xlabel("% of reconstruction (1 - reconstruction error)");
ylabel("Counts");
title("% of reconstruction after one atom: Random vs. EoI (ripples)");
legend("Random", "EoI reconstructed with one atom");

%Fast ripples
figure(2);
histogram(percRec_FR_rand, 20);
hold on
histogram(percRec_FR_sel, 20);
% histogram(percRec_FR, 20);
% histogram(percRec_FR_sel_amp, 20);
xline(mean(percRec_FR_rand)+sigmaREC*std(percRec_FR_rand), 'r', strcat(num2str(sigmaREC) ," \times \sigma"));
hold off
xlabel("% of reconstruction");
ylabel("Counts");
title("% of reconstruction after one atom: Random vs. EoI (fast ripples)");
legend("Random", "EoI reconstructed with one atom");

%% Save results 
cd(path_to_pat)

% Ripples
R = {[1:length(LocAndCHR)]',LocAndCHR(:,1),LocAndCHR(:,2),RnRec',Rlocatom,Rcoef,Rerror,RV,inc_AMP_R,th_AMP_R,inc_REC_R,th_REC_R};
savein_R = strcat("R_",filttype,"_",filenamesave);
save(savein_R, "R");

% Fast ripples
FR = {[1:length(LocAndCHFR)]',LocAndCHFR(:,1),LocAndCHFR(:,2),FRnRec',FRlocatom,FRcoef,FRerror,FRV,inc_AMP_FR,th_AMP_FR,inc_REC_FR,th_REC_FR};
savein_FR = strcat("FR_",filttype,"_",filenamesave);
save(savein_FR, "FR");

% Spikes
savein_Spikes = strcat("Spikes_",filttype,"_",filenamesave);
save(savein_Spikes, "LocAndCHSpike");

cd(s);

end % filter 

intnumb = intnumb+1;
end % interval 

nightnumb = nightnumb+1;
cd(path_to_Data);
end % night

%end

