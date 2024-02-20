function [all_snipraw, all_snipEEG, all_snipHFO, location, channel, threshold] = initial_detection(rawdata, min_freq, fs, Nc, snipsize)
% Takes as input the raw data (5 minutes) in bipolar montage (23 channels)

% Filter in the two bands
load FIR_2kHz;
load PBfilt_HFO.mat;
load PBfilt_EEG.mat;
RfilteredData =  filtfilt(filter.Rb, filter.Ra, rawdata(:,:)')'; 
FRfilteredData =  filtfilt(filter.FRb, filter.FRa, rawdata(:,:)')'; 
fulldatafilt = {RfilteredData FRfilteredData}; 
data_filt_EEG = filtfilt(PBfilt_EEG, 1, rawdata(:,:)')';
data_filt_HFO = filtfilt(PBfilt_HFO, 1, rawdata(:,:)')'; 

% Compute SD for threshold for each channel separatly and the time index of
% threshold crossing 
channel = [];
location = [];
threshold = [];
all_snipraw = [];
all_snipEEG = [];
all_snipHFO = [];

for i = 1:length(min_freq) % two different filter bands
    data = fulldatafilt{1,i};
    for ch = 1:size(data,1) % for one channel
        snippets = [];
        SDs = buffered_stats(data(ch,:),200,100,'std');
        th = 3*median(SDs);
        N = round(1/min_freq(i)*fs);
        h = ones(Nc-1,1);
        index_HFO = [];
        [~,~,ix] = zerocrossrate((data(ch,:)-th)); % logical of positions where we have a crossing
        index_cross = find(ix); % find indicies of crossing 
        if ~isempty(index_cross)
            dist_cross = diff(index_cross); % distance of indicies of crossing 
            if ~isempty(dist_cross)
                x_HFO = conv(dist_cross<N,h); % makes sure distance is smaller than 25 within HFO
                n_HFO = sum(x_HFO>=(Nc-1)); % number of detected HFO
                index_HFO = x_HFO(3:end);
                index_HFO = find(index_HFO == 3);
                index_HFO_final = index_cross(index_HFO); % time index of HFO in filtered signal       
            end

        end
        location = [location, index_HFO_final]; % save index of HFO
        channel = [channel, repelem(ch,n_HFO)]; % save channel of HFO 
        threshold = [threshold, repelem(th, n_HFO)]; % save threshold for channel

 % Put snippets of this channel into list --> switch to choose in which
 % filter we want the snippets: raw, slow waves (4-80), HFO (80 -250) 
        for type = 1:3
            switch type
                 case 1
                    targetdata = rawdata;
                 case 2
                    targetdata = data_filt_EEG; 
%                 targetdata = zeros(size(rawdata,1), size(rawdata,2));
%                 for j = 1:size(rawdata,1)
%                     targetdata(j,:) = bandpass(rawdata(j,:),[4 80],fs);
                 case 3
                    targetdata = data_filt_HFO; 
%                 targetdata = zeros(size(rawdata,1), size(rawdata,2));
%                 for k = 1:size(rawdata,1)
%                     targetdata(k,:) = bandpass(rawdata(k,:),[4 80],fs);        
            end 
            % Create snippets
            snippets = zeros(n_HFO,snipsize);
            for event = 1:length(index_HFO_final) 
                if index_HFO_final(event)<=(snipsize/2) % if it's at the beginning of the recording, take the beginning 
                    snippets(event,:) = targetdata(ch,(1:snipsize));
                elseif index_HFO_final(event)>(size(data,2)-(snipsize/2)) % if it's at the end of the recording, take the end
                    snippets(event,:) = targetdata(ch,(end-snipsize+1):end);
                else % else take the index minus half the snippet size until index plus half the snippets size 
                    snippets(event,:) = targetdata(ch,(index_HFO_final(event)-(snipsize/2)):(index_HFO_final(event)+(snipsize/2)-1));
                end
            end 
            switch type
                case 1
                    all_snipraw = [all_snipraw; snippets];
                case 2
                    all_snipEEG = [all_snipEEG; snippets];
                case 3
                    all_snipHFO = [all_snipHFO; snippets];  
            end
        end
    end
end