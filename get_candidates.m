function [finalcount, valsnip] = get_candidates(all_snip_list, channelnumber)
% get list of channels with choosen snippets and count of events per
% channel
valsnip = all_snip_list(all_snip_list(:,end) == 1,:);
channelcount = [unique(valsnip(:,1)), histc(valsnip(:,1), unique(valsnip(:,1)))];
channeln = [1:channelnumber]';
finalcount = [channeln zeros(size(channeln,1),1)];
for i = 1:size(channelcount,1)
    finalcount(channelcount(i),2) = channelcount(i,2);
end