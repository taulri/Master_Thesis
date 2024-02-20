function [OMP, V] = OMP_reconst_draw(dat,Dict,fs,Nc,minHz,th)

% Input: dat: raw event
%        Dict: Dictionary of atoms (Create_Dictionary function for
%        details)
%        fs: Sampling frequency
%        Nc: Number of crossing the threshold
%        minHz: Lower boundary of the HFO interval 
%        th: threshold for that snippet 

% Initial parameters
stop.Res = 0.05; % Stopping criteria-1
stop.iter = 10; % Stopping criteria-2

t = 1/fs:1/fs:length(dat)/fs; % time index
rng = 1.05*max(abs(dat)); % range of raw event
DLSize = size(Dict.DL,2);% low band dictionary size
frq = [Dict.frq.DL Dict.frq.DR Dict.frq.DF];% Dictionary frequency
Dictionary = [Dict.DL Dict.DR Dict.DF];%Dictionary atoms

% OMP Process for one snippet
[OMP.reconstructed,OMP.coeff,OMP.loc,OMP.Residual,OMP.Error]...
    = OMP_process(Dictionary,dat,stop.iter,stop.Res, fs, minHz, Nc, th);
% Compute the V-Factor
V = zeros(1, stop.iter);
for k=1:min(stop.iter, size(OMP.Residual,2))
    [~,~,V(k)] = VFactor(OMP.Residual(:,k));
end
 
% figure;
% % tcl = tiledlayout(2,1);
% % title(tcl, 'OMP decomposition', 'Fontsize', 26, 'FontName', 'Calibri');
% k = 1;
% subplot(2,1,1);%plot signal and approximation
% plot(t,dat,'b','LineWidth',1.2);
% xlim([1/fs length(dat)/fs]);
% %yticks([-floor(rng) -floor(rng/2) 0 ceil(rng/2) ceil(rng)]);
% xticks([]);
% ylim([-rng rng]);
% hold on
% plot(t,OMP.reconstructed(:,k),'-r','LineWidth',1.2); % atom
% plot(t,th,'m');
% set(gca,'FontSize',22, 'FontName','Calibri');
% ylabel('Amplitude [\muV]', FontSize=22);
% hold off;
% 
% subplot(2,1,2);%plot the Residual waveform
% plot(t,OMP.Residual(:,k+1),'b','LineWidth',1.2);
% hold on;
% plot(t,th,'m');
% hold off;
% xlim([1/fs length(dat)/fs]);
% set(gca,'FontSize',22, 'FontName','Calibri');
% xlabel('Time [s]', FontSize=22);
% ylabel('Amplitude [\muV]', FontSize=22);


end


    

