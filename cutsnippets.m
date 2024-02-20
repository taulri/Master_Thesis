function [snippets] = cutsnippets(data, n_snip, snipsize, index, channel)

snippets = zeros(n_snip,snipsize);
for event = 1:length(index) 
    if index(event)<=(snipsize/2) % if it's at the beginning of the recording, take the beginning 
        snippets(event,:) = data(channel(event),(1:snipsize));
        index(event) = (snipsize/2)+1;
    elseif index(event)>(size(data,2)-(snipsize/2)) % if it's at the end of the recording, take the end
        snippets(event,:) = data(channel(event),(end-snipsize+1):end);
        index(event) = size(data,2)-snipsize/2;
    else % else take the index minus half the snippet size until index plus half the snippets size 
        snippets(event,:) = data(channel(event),(index(event)-(snipsize/2)):(index(event)+(snipsize/2)-1));
    end
end 
