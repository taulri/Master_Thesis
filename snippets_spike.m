function [snipSpike, channel, threshold, location] = snippets_spike(data, snipsize, fs)

channel = [];
location = [];
threshold = [];
snipSpike = [];
overlap = 499;
window = 500;
%th = 1600; 

for ch = 1:size(data,1) % for one channel in one filter band

    % Calculate length in one channel 
    segments = buffer(data(ch,:),window, overlap);
    segments = segments(:,window:end);
    len = zeros(size(segments,2),1);
    for i = 1:length(len)
        seg = segments(:,i);
        dy = abs(diff(seg));
        len(i) = sum(dy);
    end
    %plot(len);
    th = 30*std(len);
    ix = len > th;
    n_Sp = nnz(ix);
    index_Sp = find(ix)';

    location = [location, index_Sp]; % save index of HFO
    channel = [channel, repelem(ch,n_Sp)]; % save channel of HFO 
    threshold = [threshold; th]; % save threshold for channel

    % Create snippets 
    snippets = zeros(n_Sp,snipsize);
    for event = 1:length(index_Sp) 
        if index_Sp(event)<=(snipsize/2) % if it's at the beginning of the recording, take the beginning 
            snippets(event,:) = data(ch,(1:snipsize));
        elseif index_Sp(event)>(size(data,2)-(snipsize/2)) % if it's at the end of the recording, take the end
            snippets(event,:) = data(ch,(end-snipsize+1):end);
        else % else take the index minus half the snippet size until index plus half the snippets size 
            snippets(event,:) = data(ch,(index_Sp(event)-(snipsize/2)):(index_Sp(event)+(snipsize/2)-1));
        end
    end 
    snipSpike = [snipSpike; snippets];

end

end
