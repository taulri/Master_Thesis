% view_data
%
% This script plots raw and filtered data of the events.
% Click on the figure to proceed to the next event.
%
% It calls
% DB_HFO_create_bipolar.m
% to create the bipolar montage and
% FIR_2kHz.mat
% to filter in the R and FR frequency bands.
%
% It was written in Matlab 2016b and uses datetime (optional)
%
% When using the data, please cite the publication:
%
% Fedele T, Burnos S, Boran E, Krayenbühl N, Hilfiker P, Grunwald T, Sarnthein J.
% Resection of high frequency oscillations predicts seizure outcome
% in the individual patient. Scientific Reports. 2017;7(1):13836.
%
% and the data set:
%
% Fedele T, Burnos S, Boran E, Krayenbühl N, Hilfiker P, Grunwald T, Sarnthein J.
% (2017); “High frequency oscillations detected in the intracranial EEG of 
% epilepsy patients during interictal sleep, patients’ electrode location and 
% outcome of epilepsy surgery.” CRCNS.org
% http://dx.doi.org/10.6080/K06Q1VD5
%
% 171025
% tommaso.fedele@usz.ch
% johannes.sarnthein@usz.ch
%


%% 1. load a data set
if 1
    clear all
    close all
    clc
    
    load FIR_2kHz
    fname = 'pat_01_night_01_interval_01.mat'; % some arbitrary interval
    % fname = 'iEEG_sleepHFO\pat_07_night_01_interval_01.mat'; % some arbitrary interval
    % fname = 'iEEG_sleepHFO\pat_08_night_01_interval_03.mat'; % some arbitrary interval
    % fname = 'iEEG_sleepHFO\pat_20_night_04_interval_29.mat'; % some arbitrary interval   
    % fname = sprintf('pat_%02d_night_%02d_interval_%02d',data.patient_number,data.night,data.interval)
    fprintf('now loading %s \n',fname)
    load(fname)
end

%% 2. re-reference to bipolar montage
data_bip = DB_HFO_create_bipolar(data);


%% 3. select and plot event

for iev = 1:100 % some arbitrary indices
    
    switch 3
        case 1 % ripple (R)
            ch     = data_bip.R(iev,1);
            estart = data_bip.R(iev,2);
            estop  = data_bip.R(iev,3);
            etype  = 'R';
        case 2 % fast ripple (FR)
            ch     = data_bip.FR(iev,1);
            estart = round(data_bip.FR(iev,2));
            estop  = round(data_bip.FR(iev,3));
            etype  = 'FR';
        case 3 % simultaneous occurrence of ripple and fast ripple (FRandR)
            ch     = data_bip.FRandR(iev,1);
            estart = round(data_bip.FRandR(iev,2));
            estop  = round(data_bip.FRandR(iev,3));
            etype  = 'FRandR';
    end
    
    dt = 1/data_bip.fs;
    t = dt:dt:dt*length(data_bip.x);
    
    figure('units','normalized','outerposition',[0 0  1 1])
    ax(1) = subplot(311);
    plot(t, data_bip.x(ch,:))
    ylabel('Raw Amplitude [\muV]')
    hold on;
    plot(t(estart:estop),data_bip.x(ch,round(estart):round(estop)),'Color','r');
    if exist('datetime') %#ok<EXIST>
        event_onset = data.Rec_datetime + seconds(data.Lag_seconds) + seconds(estart/data.fs);
        title(sprintf('Patient %02d, Night %02d, Interval %02d, %s #%i,  Channel %s,  Event onset time %s',...
            data.patient_number, data.night, data.interval,etype,iev,data_bip.label{ch},...
            datestr(event_onset,'yyyy-mm-dd hh:MM:ss'))) % .SSS gives error
    else
        title(sprintf('Patient %02d, Night %02d, Interval %02d, %s #%i,  Channel %s,  Recording onset time %s',...
            data.patient_number, data.night, data.interval,etype,iev,data_bip.label{ch},...
            datestr(data.Rec_timestamp)))
    end
    
    ax(2) = subplot(312);
    RfilteredData =  filtfilt(filter.Rb, filter.Ra, data_bip.x(ch,:)')';
    plot(t, RfilteredData);
    ylabel('Ripple Amplitude [\muV]')
    hold on;
    plot(t(estart:estop),RfilteredData(round(estart):round(estop)),'Color','r');
      
    ax(3) = subplot(313);
    FRfilteredData =  filtfilt(filter.FRb, filter.FRa, data_bip.x(ch,:)')';
    hilbFilt = abs(hilbert(FRfilteredData));
    plot(t, FRfilteredData);
    ylabel('Fast Ripple Amplitude [\muV]')
    xlabel('Time in interval [s]')
    hold on;
    plot(t(estart:estop),FRfilteredData(round(estart):round(estop)),'Color','r');
    linkaxes(ax,'x')
   
    xlim([t(estart)-.5  t(estop)+.5])
    ylim([-20 +20])
    
    line([t(estart)   t(estop) ],[-10 -10],'color','r')
    
    disp('Please click on figure to proceed to next event')
    try
        waitforbuttonpress
    catch
        return;
    end
    close all
    
end