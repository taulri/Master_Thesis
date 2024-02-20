function [snippets_raw, snippets_R, snippets_FR, channel, location] = snippets_spike_detector(data_raw, data_R, data_FR, snipsize, fs)

%% Running spike detector
settings = '-bl 10 -bh 60 -h 60 -jl 3.65 -dec 200'; % default settings of Janca detector
out = spike_detector_hilbert_v25(data_raw',fs,settings);

%% Snippets
channel = out.chan;
location = round(out.pos*fs);

snippets_raw = zeros(length(location),snipsize);
snippets_R = zeros(length(location),snipsize);
snippets_FR = zeros(length(location),snipsize);

    for event = 1:length(location) 
        if location(event)<=(snipsize/2) % if it's at the beginning of the recording, take the beginning 
            snippets_raw(event,:) = data_raw(channel(event),(1:snipsize));
            snippets_R(event,:) = data_R(channel(event),(1:snipsize));
            snippets_FR(event,:) = data_FR(channel(event),(1:snipsize));
            location(event) = (snipsize/2)+1;

        elseif location(event)>(size(data_raw,2)-(snipsize/2)) % if it's at the end of the recording, take the end
            snippets_raw(event,:) = data_raw(channel(event),(end-snipsize+1):end);
            snippets_R(event,:) = data_R(channel(event),(end-snipsize+1):end);
            snippets_FR(event,:) = data_FR(channel(event),(end-snipsize+1):end);
            location(event) = size(data_raw,2)-snipsize/2;

        else % else take the index minus half the snippet size until index plus half the snippets size 
            snippets_raw(event,:) = data_raw(channel(event),(location(event)-(snipsize/2)):(location(event)+(snipsize/2)-1));
            snippets_R(event,:) = data_R(channel(event),(location(event)-(snipsize/2)):(location(event)+(snipsize/2)-1));
            snippets_FR(event,:) = data_FR(channel(event),(location(event)-(snipsize/2)):(location(event)+(snipsize/2)-1));
        end
    end 
end


