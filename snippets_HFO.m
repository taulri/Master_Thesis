function [snipR, snipFR, LocAndCHR, LocAndCHFR, threshold] = snippets_HFO(data_R_FR, min_freq, fs, Nc, snipsize)

channel = [];
location = [];
threshold = [];
snipR = [];
snipFR = [];

for i = 1:length(min_freq) % two different filter bands
    data = data_R_FR{1,i};
    for ch = 1:size(data,1) % for one channel in one filter band
        SDs = buffered_stats(data(ch,:),2000,1000,'std'); % frame = 1s, overlap = 0.5 s
        th = repelem(3*SDs,1000); % th is adaptive in 0.5 s 
        N = round(1/min_freq(i)*fs); % Minimal distance so that it's not the same event 
        h = ones(Nc-1,1); % Array with minimal number of crossings as ones 
        index_HFO = [];
        [~,~,ix] = zerocrossrate((data(ch,:)-th)); % logical of positions where we have a crossing
        index_cross = find(ix); % find indicies of crossing 
        if ~isempty(index_cross) % if there are crossings
            dist_cross = diff(index_cross); % distance of indicies of crossing 
            if ~isempty(dist_cross)
                x_HFO = conv(dist_cross<N,h); % makes sure distance is smaller than 25 within HFO
                index_HFO = find(x_HFO(5:end) == Nc-1); % index of HFO in crossing array 
                index_HFO_final = index_cross(index_HFO); % time index of HFO in filtered signal       
                % Check that the indexes are not the same event
                if length(index_HFO_final)>= 2
                    index_HFO_final = index_HFO_final([1,find((diff(index_HFO_final) > N))+1]);
                end
                n_HFO = length(index_HFO_final); % number of detected HFO

            end
        end

        % Create snippets
        snippets = zeros(n_HFO,snipsize);
        for event = 1:length(index_HFO_final) 
            if index_HFO_final(event)<=(snipsize/2) % if it's at the beginning of the recording, take the beginning 
                snippets(event,:) = data(ch,(1:snipsize));
                index_HFO_final(event) = (snipsize/2)+1;
            elseif index_HFO_final(event)>(size(data,2)-(snipsize/2)) % if it's at the end of the recording, take the end
                snippets(event,:) = data(ch,(end-snipsize+1):end);
                index_HFO_final(event) = size(data,2)-snipsize/2;
            else % else take the index minus half the snippet size until index plus half the snippets size 
                snippets(event,:) = data(ch,(index_HFO_final(event)-(snipsize/2)):(index_HFO_final(event)+(snipsize/2)-1));
            end
        end 
        if i == 1
            snipR = [snipR; snippets];
        else
            snipFR = [snipFR; snippets];
        end
        location = [location, index_HFO_final]; % save index of HFO
        channel = [channel, repelem(ch,n_HFO)]; % save channel of HFO 
        threshold = [threshold; th]; % save threshold for channel
    end
end

% Snippets R and FR
LocAndCh = zeros(2,length(channel));
LocAndCh(1,:) = channel;
LocAndCh(2,:) = location;
LocAndCHR = LocAndCh(:,1:size(snipR,1))';
LocAndCHFR = LocAndCh(:,size(snipR,1)+1:end)';

% % Snippets RFR
% locRFR = [];
% chRFR = [];
% for q = 1:size(LocAndCHR,1)
%   spot = ((LocAndCHFR(:,1) >= LocAndCHR(q,1)-(snipsize/2)) + (LocAndCHFR(:,1) <= LocAndCHR(q,1)+(snipsize/2)));
%   index = find(spot == 2);
%   ch = LocAndCHR(q,2);
%   indexmatch = find(LocAndCHFR(index,2) == ch);
%   locRFR = [locRFR, LocAndCHFR(index(indexmatch),1)'];
%   chRFR = [chRFR, LocAndCHFR(index(indexmatch),2)'];
% end
% LocAndCHRFR = [locRFR;chRFR]';
% 
% 
% snipRandFR_R = zeros(size(LocAndCHRFR,1),snipsize);
% for event = 1:size(LocAndCHRFR,1)
%      if LocAndCHRFR(event,1)<=(snipsize/2) % if it's at the beginning of the recording, take the beginning 
%           snipRandFR_R(event,:) = data_R_FR{1,1}(LocAndCHRFR(event,2),1:snipsize);
%      elseif LocAndCHRFR(event,1)>(size(data,2)-(snipsize/2)) % if it's at the end of the recording, take the end
%           snipRandFR_R(event,:) = data_R_FR{1,1}(LocAndCHRFR(event,2),(end-snipsize+1):end);
%      else % else take the index minus half the snippet size until index plus half the snippets size 
%           snipRandFR_R(event,:) = data_R_FR{1,1}(LocAndCHRFR(event,2),LocAndCHRFR(event,1)-(snipsize/2):LocAndCHRFR(event,1)+(snipsize/2)-1);
%      end
%  end 
% 
% snipRandFR_FR = zeros(size(LocAndCHRFR,1),snipsize);
% for event = 1:size(LocAndCHRFR,1)
%      if LocAndCHRFR(event,1)<=(snipsize/2) % if it's at the beginning of the recording, take the beginning 
%           snipRandFR_FR(event,:) = data_R_FR{1,2}(LocAndCHRFR(event,2),1:snipsize);
%      elseif LocAndCHRFR(event,1)>(size(data,2)-(snipsize/2)) % if it's at the end of the recording, take the end
%           snipRandFR_FR(event,:) = data_R_FR{1,2}(LocAndCHRFR(event,2),(end-snipsize+1):end);
%      else % else take the index minus half the snippet size until index plus half the snippets size 
%           snipRandFR_FR(event,:) = data_R_FR{1,2}(LocAndCHRFR(event,2),LocAndCHRFR(event,1)-(snipsize/2):LocAndCHRFR(event,1)+(snipsize/2)-1);
%      end
% end 
end
