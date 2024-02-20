function plot_signal(rawdata, snippetlist, fs)
% plots trace of some channels with the indicies of the snippets
n = 1;
t = 1/fs:1/fs:10000/fs;
for channel = 1:3 % for each channel 
    plot(t,rawdata(channel,1:10000),'black');
    xlabel('Time [s]')
    ylabel('Raw Amplitude [uV]')
%     title(['Channel ',num2str(channel)]);
%     while snippetlist(n,1) == channel
%         xline(snippetlist(n,2)/fs, 'r');
%         n = n+1;
%     end
end
