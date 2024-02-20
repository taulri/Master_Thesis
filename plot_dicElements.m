function [] = plot_dicElements(ele1)

plot(1:512, ele1);
hold on;
xlim([1 514])
set(gca,'XTick',1:64:514)
set(gca,'XTickLabel',round(((1)/2000),2):round(((64)/2000),2):round(((512)/2000),2)) % change xlab to s        
title('Dictionary element');
ylabel('Amplitude [\muV]');
xlabel('Time [s]');
hold off;

end