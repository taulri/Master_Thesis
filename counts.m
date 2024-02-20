function [finalcount] = counts(channels, channelnumber)

channelcount = [unique(channels), histc(channels, unique(channels))];
finalcount = [[1:channelnumber]' zeros(channelnumber,1)];
for i = 1:size(channelcount,1)
    finalcount(channelcount(i),2) = channelcount(i,2);
end