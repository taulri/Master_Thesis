function plot_snippets(snippetlist, Data, Data2, whattoplot, name, th)
% plots snippets in filtered trace 

switch whattoplot
    case 'r'
        z = randi(size(snippetlist{1},1),10,1);
    case 'all'
        z = 1:size(snippetlist{1},1);
end

for i = 1:length(z) %how many snippets
    figure(i);
    r = z(i);
    channel = snippetlist{2}(r,1);
    loc = snippetlist{1}(r,1);
    %threshold = snippetlist{3}(:,channel);
    if loc < 512
        loc = 513;
    end
    plot(loc-512:loc+512, Data(channel,loc-512:loc+512), 'black');
    hold on;

    plot(loc-512:loc+512, th(channel,loc-512:loc+512), 'red');
    if strcmp(name, 'Ripple and Fast ripple') || strcmp(name, 'Spike and Ripple')
        plot(loc-512:loc+512, Data2(channel,loc-512:loc+512), 'red');
    end
    %plot(loc-512:loc+512, threshold(loc-512:loc+512,1), 'red');
    xlim([loc-512 loc+512])
    set(gca,'XTick',loc-512:200:loc+512)
    set(gca,'XTickLabel',round(((loc-512)/2000),2):0.01:round(((loc+512)/2000),2)) % change xlab to s        
    title(strcat('Snippet ', num2str(r), name));
    ylabel('Amplitude [\muV]');
    xlabel('Time [s]');
    hold off;
end



