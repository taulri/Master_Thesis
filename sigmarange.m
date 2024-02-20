function [sigmarange, inc] = sigmarange(snippets, alpha, OMP_rec)

% Calculating the range within apha* std() of the signal 

sigma = std(snippets,0,2); 
meanamp = mean(snippets,2);
sigmarange = [meanamp-(alpha*sigma),meanamp+(alpha*sigma)]; 
inc = zeros(length(OMP_rec),1);
for i = 1:length(OMP_rec)
    inc(i) = all(sigmarange(i,1) < OMP_rec{i}.Residual(:,2)) && all(OMP_rec{i}.Residual(:,2) < sigmarange(i,2));
end

figure(1);
n = 30; 
h1 = histogram(snippets(n,:));
hold on
h2 = histogram(OMP_rec{n}.Residual(:,2));
xline(sigmarange(n,1), 'r', strcat("-", num2str(alpha)," \times \sigma"));
xline(sigmarange(n,2), 'r', strcat(num2str(alpha)," \times \sigma"));
xline(meanamp(n), 'r', '\mu');

hold off
xlabel("Amplitude [\muV]");
ylabel("Counts");
title("Distribution of amplitudes within snippet");
legend("Snippet before first reconstruction", "Residual after first reconstruction");

end