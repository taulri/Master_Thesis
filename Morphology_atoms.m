%% Shapes of chosen atoms 
clc
clear
close all

%% Path and parameters
s = 'Y:\NCH_FL_Forschungsprojekte\Epilepsy\_Master Students\Ulrich Tanja\Results';
addpath('\\fs-home\ulrta$\Documents\Master_Biomed\Code\DataSet\Function_Needed_23_12_20_version2\');
cd(s);

%% Get the atom IDs of detectors

% We chose detectors using R or FR
band = {'R', 'FR'}; 
filter = {'Ref'};
criteria = {'all','REC','AMP'}; % 'all', 'AMP', 'REC'

conditions = cell(1,6);
conditionsname = cell(1,6);

i = 0;

for HFOband =  1:length(band)
    for filt = 1:length(filter)
        for crit = 1:length(criteria)
            i = i+1;
            conditionsname{1,i} = strcat(filter{filt},"\",band{HFOband},"\",criteria{crit});
            for pat = 18
                cd(strcat(s,"\Pat",num2str(pat),"\",filter{filt},"\",band{HFOband},"\",criteria{crit}));
                files = dir("*_AtomCoef.mat");
                for fi = 1:length(files)
                    data = load(files(fi).name); 
                    conditions{1,i} = [conditions{1,i} data.AtomCoef];                    
                end 
            end  
        end
    end
end




%% The dictionary

N = 512; % Length of dictionary
[Dictlocinfo,DL,DR,DF] = Create_Dictionary(N,2000); % create the Gabor dictionary at 2000 Hz sampling frequency
Dict.DL = DL.Atom; % High
Dict.DR = DR.Atom; % Medium 
Dict.DF = DF.Atom; % Low
Dict.frq.DL = DL.frq;
Dict.frq.DR = DR.frq;
Dict.frq.DF = DF.frq;


Locatoms = [Dictlocinfo.Iloc, Dictlocinfo.Jloc, Dictlocinfo.Kloc, Dictlocinfo.Lloc, Dictlocinfo.Eloc, Dictlocinfo.Floc, Dictlocinfo.Gloc, Dictlocinfo.Hloc,Dictlocinfo.Aloc,Dictlocinfo.Bloc,Dictlocinfo.Cloc,Dictlocinfo.Dloc];
frequencies = [DL.frq DR.frq DF.frq];
windowing = [Dictlocinfo.Iwnd Dictlocinfo.Jwnd Dictlocinfo.Kwnd Dictlocinfo.Lwnd Dictlocinfo.Ewnd Dictlocinfo.Fwnd Dictlocinfo.Gwnd Dictlocinfo.Hwnd Dictlocinfo.Awnd Dictlocinfo.Bwnd Dictlocinfo.Cwnd Dictlocinfo.Dwnd];

%% Groups

% Start Dict with high ones, DL (I, J, K, L)
sizeI = size(Dictlocinfo.I,2);
part1I = [1:sizeI/4 (2*sizeI/4)+1:3*sizeI/4];
part2I = [(sizeI/4)+1:2*sizeI/4 (3*sizeI/4)+1:4*sizeI/4];

sizeJ = size(Dictlocinfo.J,2);
part1J = [1:sizeJ/4 (2*sizeJ/4)+1:3*sizeJ/4] + sizeI;
part2J = [(sizeJ/4)+1:2*sizeJ/4 (3*sizeJ/4)+1:4*sizeJ/4] + sizeI;

sizeK = size(Dictlocinfo.K,2); % 2
partK = [1:sizeK] + sizeI + sizeJ;

sizeL = size(Dictlocinfo.L,2); % 2
partL = [1:sizeL] + sizeI + sizeJ + sizeK;


% Continue with medium ones, DR (E, F, G, H)

sizeE = size(Dictlocinfo.E,2);
part1E = [1:sizeE/4 (2*sizeE/4)+1:3*sizeE/4] + sizeI + sizeJ + sizeK +sizeL;
part2E = [(sizeE/4)+1:2*sizeE/4 (3*sizeE/4)+1:4*sizeE/4] + sizeI + sizeJ + sizeK +sizeL;

sizeF = size(Dictlocinfo.F,2);
part1F = [1:sizeF/4 (2*sizeF/4)+1:3*sizeF/4] + sizeI + sizeJ + sizeK +sizeL +sizeE;
part2F = [(sizeF/4)+1:2*sizeF/4 (3*sizeF/4)+1:4*sizeF/4] + sizeI + sizeJ + sizeK +sizeL + sizeE;

sizeG = size(Dictlocinfo.G,2); % 2
partG = [1:sizeG] + sizeI + sizeJ + sizeK +sizeL + sizeE + sizeF;

sizeH = size(Dictlocinfo.H,2); % 2
partH = [1:sizeH] + sizeI + sizeJ + sizeK +sizeL + sizeE + sizeF + sizeG;



% Continue with low ones, DF (A, B, C, D)

sizeA = size(Dictlocinfo.A,2);
part1A = [1:sizeA/4 (2*sizeA/4)+1:3*sizeA/4] + sizeI + sizeJ + sizeK +sizeL + sizeE + sizeF + sizeG +sizeH;
part2A = [(sizeA/4)+1:2*sizeA/4 (3*sizeA/4)+1:4*sizeA/4] + sizeI + sizeJ + sizeK +sizeL + sizeE + sizeF + sizeG +sizeH;

sizeB = size(Dictlocinfo.B,2);
part1B = [1:sizeB/4 (2*sizeB/4)+1:3*sizeB/4] + sizeI + sizeJ + sizeK +sizeL + sizeE + sizeF + sizeG +sizeH+sizeA;
part2B = [(sizeB/4)+1:2*sizeB/4 (3*sizeB/4)+1:4*sizeB/4] + sizeI + sizeJ + sizeK +sizeL + sizeE + sizeF + sizeG +sizeH+sizeA;

sizeC = size(Dictlocinfo.C,2);
part1C = [1:sizeC/4 (2*sizeC/4)+1:3*sizeC/4]+ sizeI + sizeJ + sizeK +sizeL + sizeE + sizeF + sizeG +sizeH+sizeA+sizeB;
part2C = [(sizeC/4)+1:2*sizeC/4 (3*sizeC/4)+1:4*sizeC/4]+ sizeI + sizeJ + sizeK +sizeL + sizeE + sizeF + sizeG +sizeH+sizeA+sizeB;

sizeD = size(Dictlocinfo.D,2);
part1D = [1:sizeD/4 (2*sizeD/4)+1:3*sizeD/4]+ sizeI + sizeJ + sizeK +sizeL + sizeE + sizeF + sizeG +sizeH+sizeA+sizeB+sizeC;
part2D = [(sizeD/4)+1:2*sizeD/4 (3*sizeD/4)+1:4*sizeD/4]+ sizeI + sizeJ + sizeK +sizeL + sizeE + sizeF + sizeG +sizeH+sizeA+sizeB+sizeC;


% % Low frequency 
% 
% sizeA = size(Dictlocinfo.A,2);
% part1A = [1:sizeA/4 (2*sizeA/4)+1:3*sizeA/4];
% part2A = [(sizeA/4)+1:2*sizeA/4 (3*sizeA/4)+1:4*sizeA/4];
% 
% sizeB = size(Dictlocinfo.B,2);
% part1B = [1:sizeB/4 (2*sizeB/4)+1:3*sizeB/4] + sizeA;
% part2B = [(sizeB/4)+1:2*sizeB/4 (3*sizeB/4)+1:4*sizeB/4] + sizeA;
% 
% sizeC = size(Dictlocinfo.C,2);
% part1C = [1:sizeC/4 (2*sizeC/4)+1:3*sizeC/4]+ sizeA + sizeB;
% part2C = [(sizeC/4)+1:2*sizeC/4 (3*sizeC/4)+1:4*sizeC/4]+ sizeA + sizeB;
% 
% sizeD = size(Dictlocinfo.D,2);
% part1D = [1:sizeD/4 (2*sizeD/4)+1:3*sizeD/4]+ sizeA + sizeB + sizeC;
% part2D = [(sizeD/4)+1:2*sizeD/4 (3*sizeD/4)+1:4*sizeD/4]+ sizeA + sizeB + sizeC;
% 
% % Medium frequency 
% 
% sizeE = size(Dictlocinfo.E,2);
% part1E = [1:sizeE/4 (2*sizeE/4)+1:3*sizeE/4]+ sizeA + sizeB + sizeC + sizeD;
% part2E = [(sizeE/4)+1:2*sizeE/4 (3*sizeE/4)+1:4*sizeE/4]+ sizeA + sizeB + sizeC + sizeD;
% 
% sizeF = size(Dictlocinfo.F,2);
% part1F = [1:sizeF/4 (2*sizeF/4)+1:3*sizeF/4]+ sizeA + sizeB + sizeC + sizeD + sizeE;
% part2F = [(sizeF/4)+1:2*sizeF/4 (3*sizeF/4)+1:4*sizeF/4]+ sizeA + sizeB + sizeC + sizeD + sizeE;
% 
% sizeG = size(Dictlocinfo.G,2); % 2
% partG = [1:sizeG]+ sizeA + sizeB + sizeC + sizeD + sizeE + sizeF;
% 
% sizeH = size(Dictlocinfo.H,2); % 2
% partH = [1:sizeH]+ sizeA + sizeB + sizeC + sizeD + sizeE + sizeF + sizeG;
% 
% % High frequency 
% 
% sizeI = size(Dictlocinfo.I,2);
% part1I = [1:sizeI/4 (2*sizeI/4)+1:3*sizeI/4]+ sizeA + sizeB + sizeC + sizeD + sizeE + sizeF + sizeG + sizeH;
% part2I = [(sizeI/4)+1:2*sizeI/4 (3*sizeI/4)+1:4*sizeI/4]+ sizeA + sizeB + sizeC + sizeD + sizeE + sizeF + sizeG + sizeH;
% 
% sizeJ = size(Dictlocinfo.J,2);
% part1J = [1:sizeJ/4 (2*sizeJ/4)+1:3*sizeJ/4]+ sizeA + sizeB + sizeC + sizeD + sizeE + sizeF + sizeG + sizeH + sizeI;
% part2J = [(sizeJ/4)+1:2*sizeJ/4 (3*sizeJ/4)+1:4*sizeJ/4]+ sizeA + sizeB + sizeC + sizeD + sizeE + sizeF + sizeG + sizeH + sizeI;
% 
% sizeK = size(Dictlocinfo.K,2); % 2
% partK = [1:sizeK]+ sizeA + sizeB + sizeC + sizeD + sizeE + sizeF + sizeG + sizeH + sizeI + sizeJ;
% 
% sizeL = size(Dictlocinfo.L,2); % 2
% partL = [1:sizeL]+ sizeA + sizeB + sizeC + sizeD + sizeE + sizeF + sizeG + sizeH + sizeI + sizeJ + sizeK;



%% Counts per group 

for cond = 1:length(conditions)
    
    p1A = 0; p2A = 0; p1B = 0; p2B = 0; p1C = 0; p2C = 0; p1D = 0; p2D = 0; p1E = 0; p2E = 0; p1F = 0; p2F = 0; pG = 0; pH = 0; 
    p1I = 0; p2I = 0; p1J = 0; p2J = 0; pK = 0; pL = 0;
    
    freq = zeros(1,size(conditions{1,cond},2));
    wnd = zeros(1,size(conditions{1,cond},2));
    r = zeros(1,size(conditions{1,cond},2));
    amp = conditions{1,cond}(2,:);
    
    for i = 1:size(conditions{1,cond},2)
        freq(i) = frequencies(conditions{1,cond}(1,i));
        wnd(i) = windowing(conditions{1,cond}(1,i));
       
    
        if ismember(conditions{1,cond}(1,i), part1A)
            p1A = p1A + 1;
            r(i) = 0.2;
        elseif ismember(conditions{1,cond}(1,i), part2A)
            p2A = p2A + 1;
            r(i) = 0.2;
    
        elseif ismember(conditions{1,cond}(1,i), part1B)
            p1B = p1B + 1;
            r(i) = 0.8;
        elseif ismember(conditions{1,cond}(1,i), part2B)
            p2B = p2B + 1;
            r(i) = 0.8;
    
        elseif ismember(conditions{1,cond}(1,i), part1C)
            p1C = p1C + 1;
            r(i) = 0.2;
        elseif ismember(conditions{1,cond}(1,i), part2C)
            p2C = p2C + 1;
            r(i) = 0.2;
    
        elseif ismember(conditions{1,cond}(1,i), part1D)
            p1D = p1D + 1;
            r(i) = 0.8;
        elseif ismember(conditions{1,cond}(1,i), part2D)
            p2D = p2D + 1;
            r(i) = 0.8;
    
        elseif ismember(conditions{1,cond}(1,i), part1E)
            p1E = p1E + 1;
            r(i) = 0.2;
        elseif ismember(conditions{1,cond}(1,i), part2E)
            p2E = p2E + 1;
            r(i) = 0.2;
    
        elseif ismember(conditions{1,cond}(1,i), part1F)
            p1F = p1F + 1;
            r(i) = 0.8;
        elseif ismember(conditions{1,cond}(1,i), part2F)
            p2F = p2F + 1;
            r(i) = 0.8;
    
        elseif ismember(conditions{1,cond}(1,i), partG)
            pG = pG + 1;
            r(i) = 0.2;
    
        elseif ismember(conditions{1,cond}(1,i), partH)
            pH = pH + 1;   
            r(i) = 0.8;
    
        elseif ismember(conditions{1,cond}(1,i), part1I)
            p1I = p1I + 1;
            r(i) = 0.2;
        elseif ismember(conditions{1,cond}(1,i), part2I)
            p2I = p2I + 1;
            r(i) = 0.2;
    
        elseif ismember(conditions{1,cond}(1,i), part1J)
            p1J = p1J + 1;
            r(i) = 0.8;
        elseif ismember(conditions{1,cond}(1,i), part2J)
            p2J = p2J + 1;
            r(i) = 0.8;
    
        elseif ismember(conditions{1,cond}(1,i), partK)
            pK = pK + 1;
            r(i) = 0.2;
    
        elseif ismember(conditions{1,cond}(1,i), partL)
            pL = pL + 1;
            r(i) = 0.8;
    
        else
            "Error"
        end
    end

    conditions{1, cond}(5,:) = freq;
    conditions{1, cond}(6,:) = wnd;
    conditions{1, cond}(7,:) = r;
    conditions{1, cond}(8,:) = amp;
    
end

%% Save atom features
% 
 cd(s);
% descr = 'atomindex, coef, channel, locationindex, freq, wnd, r, amp';
 save('Atom_features', "conditions");
% save('Atom_features_conditions', "conditionsname");
% save('Table_description_atom_features', "descr");



%% Plots



% Histograms
faceColor1 = [0.500000000000000,0.500000000000000,0.500000000000000]; % grey
faceColor2 = [0.200000000000000,0.800000000000000,0.200000000000000]; % green
faceColor3 = [0.800000000000000,0.200000000000000,0.200000000000000]; % red

% Frequencies R 
figure(1);
%histogram(conditions{1,1}(5,:),'FaceColor', faceColor1, 'EdgeColor', faceColor1, 'BinWidth', 10,'FaceAlpha', 1);
bar(unique(frequencies),histc(conditions{1,1}(5,:), unique(frequencies)), 1, 'FaceColor', faceColor1, 'EdgeColor', faceColor1,'FaceAlpha', 1);
hold on;
%histogram(conditions{1,3}(5,:),'FaceColor', faceColor2, 'EdgeColor', faceColor2, 'BinWidth', 10,'FaceAlpha', 0.9);
bar(unique(frequencies),histc(conditions{1,3}(5,:), unique(frequencies)),1,'FaceColor', faceColor2, 'EdgeColor', faceColor2,'FaceAlpha', 0.9);
hold on;
%histogram(conditions{1,2}(5,:),'FaceColor', faceColor3, 'EdgeColor', faceColor3, 'BinWidth', 10,'FaceAlpha', 0.8);
bar(unique(frequencies),histc(conditions{1,2}(5,:), unique(frequencies)),1,'FaceColor', faceColor3, 'EdgeColor', faceColor3,'FaceAlpha', 0.8);

hold off;
set(gca,'FontName','Calibri','FontSize',22);
ylabel('Bin Count','FontSize',22,'FontName','Calibri');
xlabel('Frequency [Hz]','FontSize',22,'FontName','Calibri');
title('Frequency of atoms decomposing ripples','FontSize',26, 'FontName','Calibri');
legend('All EoI', 'AMP', 'REC', 'FontSize', 22, 'FontName','Calibri');
 
% Frequencies FR
figure(2);
bar(unique(frequencies),histc(conditions{1,4}(5,:), unique(frequencies)), 1, 'FaceColor', faceColor1, 'EdgeColor', faceColor1,'FaceAlpha', 1);
%histogram(conditions{1,4}(5,:),'FaceColor', faceColor1, 'EdgeColor', faceColor1, 'BinWidth', 10,'FaceAlpha', 1);
hold on;
bar(unique(frequencies),histc(conditions{1,5}(5,:), unique(frequencies)), 1, 'FaceColor', faceColor3, 'EdgeColor', faceColor3,'FaceAlpha', 0.9);
bar(unique(frequencies),histc(conditions{1,6}(5,:), unique(frequencies)), 1, 'FaceColor', faceColor2, 'EdgeColor', faceColor2,'FaceAlpha', 0.8);
%histogram(conditions{1,5}(5,:), 'FaceColor', faceColor3, 'EdgeColor', faceColor3, 'BinWidth', 10,'FaceAlpha', 0.9);
%histogram(conditions{1,6}(5,:), 'FaceColor', faceColor2, 'EdgeColor', faceColor2, 'BinWidth', 10,'FaceAlpha', 0.8);

hold off;
set(gca,'FontName','Calibri','FontSize',22);
ylabel('Bin Count','FontSize',22,'FontName','Calibri');
xlabel('Frequency [Hz]','FontSize',22,'FontName','Calibri');
title('Frequency of atoms decomposing fast ripples','FontSize',26, 'FontName','Calibri');
legend('All EoI', 'REC', 'AMP', 'FontSize', 22, 'FontName','Calibri');

% % Width R 
% % wnd:  n = 4 (0.0312    0.0625    0.1250    0.2500)
% figure(2);
% histogram(conditions{1,1}(6,:));
% hold on;
% histogram(conditions{1,2}(6,:));
% hold on;
% histogram(conditions{1,3}(6,:));
% hold off;
% 
% % Cosine fraction R 
% % r: n = 2 (0.2000    0.8000)
% figure(3);
% histogram(conditions{1,1}(7,:));
% hold on;
% histogram(conditions{1,2}(7,:));
% hold on;
% histogram(conditions{1,3}(7,:));
% hold off;

% Amplitude Ripples
% amp: n = 519'389
figure(3);
histogram(abs(conditions{1,1}(8,:)),'FaceColor', faceColor1, 'EdgeColor', faceColor1,'FaceAlpha', 1);
hold on;
histogram(abs(conditions{1,3}(8,:)), 'FaceColor', faceColor2, 'EdgeColor', faceColor2,'FaceAlpha', 0.9);
hold on;
histogram(abs(conditions{1,2}(8,:)), 'FaceColor', faceColor3, 'EdgeColor', faceColor3,'FaceAlpha', 0.8);
hold off;
set(gca,'FontName','Calibri','FontSize',22);
ylabel('Bin Count','FontSize',22,'FontName','Calibri');
xlabel('Amplitude [\muV]','FontSize',22,'FontName','Calibri');
title('Amplitude of atoms decomposing ripples','FontSize',26, 'FontName','Calibri');
legend('All EoI', 'AMP', 'REC', 'FontSize', 22, 'FontName','Calibri');

% Amplitude Fast ripple 
figure(4);
histogram(abs(conditions{1,4}(8,:)),'FaceColor', faceColor1, 'EdgeColor', faceColor1,'FaceAlpha', 1);
hold on;
hold on;
histogram(abs(conditions{1,5}(8,:)),'FaceColor', faceColor3, 'EdgeColor', faceColor3,'FaceAlpha', 0.9);
histogram(abs(conditions{1,6}(8,:)),'FaceColor', faceColor2, 'EdgeColor', faceColor2,'FaceAlpha', 0.8);

hold off;
set(gca,'FontName','Calibri','FontSize',22);
ylabel('Bin Count','FontSize',22,'FontName','Calibri');
xlabel('Amplitude [\muV]','FontSize',22,'FontName','Calibri');
title('Amplitude of atoms decomposing fast ripples','FontSize',26, 'FontName','Calibri');
legend('All EoI', 'REC', 'AMP', 'FontSize', 22, 'FontName','Calibri');

% Scatter plots 

% Frequency vs Amp
% Ripples
figure(5);
scatter(conditions{1,1}(5,:)-1, abs(conditions{1,1}(8,:)),50, faceColor1, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
scatter(conditions{1,3}(5,:), abs(conditions{1,3}(8,:)),50, faceColor2, 'filled', 'MarkerFaceAlpha', 0.5);
scatter(conditions{1,2}(5,:)+1, abs(conditions{1,2}(8,:)),50, faceColor3, 'filled', 'MarkerFaceAlpha', 0.5);
hold off;
ylabel('Amplitude [\muV]','FontSize',22,'FontName','Calibri');
xlabel('Frequency [Hz]','FontSize',22,'FontName','Calibri');
title('Frequency vs amplitude of atoms decomposing ripples','FontSize',26, 'FontName','Calibri');
legend('All EoI', 'AMP', 'REC', 'FontSize', 22, 'FontName','Calibri');
set(gca,'FontName','Calibri','FontSize',22);


% Fast ripples 
figure(6);
scatter(conditions{1,4}(5,:)-2, abs(conditions{1,4}(8,:)),30, faceColor1, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
scatter(conditions{1,6}(5,:), abs(conditions{1,6}(8,:)),30, faceColor2, 'filled', 'MarkerFaceAlpha', 0.5);
scatter(conditions{1,5}(5,:)+2, abs(conditions{1,5}(8,:)),30, faceColor3, 'filled', 'MarkerFaceAlpha', 0.5);
hold off;
set(gca,'FontName','Calibri','FontSize',22);
ylabel('Amplitude [\muV]','FontSize',22,'FontName','Calibri');
xlabel('Frequency [Hz]','FontSize',22,'FontName','Calibri');
title('Frequency vs amplitude of atoms decomposing fast ripples','FontSize',26, 'FontName','Calibri');
legend('All EoI', 'AMP', 'REC', 'FontSize', 22, 'FontName','Calibri');



% % 3D plots
% 
% figure(7);
% %plot3(abs(conditions{1,4}(6,:))-0.005, abs(conditions{1,4}(5,:)), abs(conditions{1,4}(8,:)),'o','LineWidth',1, 'color', faceColor1, 'MarkerFaceColor', faceColor1);
% scatter3((abs(conditions{1,1}(6,:))*100)-0.5, abs(conditions{1,1}(5,:)), abs(conditions{1,1}(8,:)),30, faceColor1, 'filled', 'MarkerFaceAlpha', 0.5);
% hold on;
% scatter3((abs(conditions{1,3}(6,:))*100)+0.5, abs(conditions{1,3}(5,:)), abs(conditions{1,3}(8,:)),30, faceColor2, 'filled', 'MarkerFaceAlpha', 0.5);
% scatter3((abs(conditions{1,2}(6,:))*100), abs(conditions{1,2}(5,:)), abs(conditions{1,2}(8,:)),30, faceColor3, 'filled', 'MarkerFaceAlpha', 0.5);
% 
% %plot3(abs(conditions{1,6}(6,:)), abs(conditions{1,6}(5,:)), abs(conditions{1,6}(8,:)),'o','LineWidth',1, 'color', faceColor2, 'MarkerFaceColor', faceColor2);
% %plot3(abs(conditions{1,5}(6,:))+0.005, abs(conditions{1,5}(5,:)), abs(conditions{1,5}(8,:)),'o','LineWidth',1, 'color', faceColor3, 'MarkerFaceColor', faceColor3);
% hold off;
% legend('All EoI', 'AMP', 'REC', 'FontSize', 22, 'FontName','Calibri');
% set(gca,'FontName','Calibri','FontSize',22);
% ylabel('Frequency [Hz]','FontSize',22,'FontName','Calibri');
% xlabel('Oscillatory fraction of the atom [%]','FontSize',22,'FontName','Calibri');
% zlabel('Amplitude [\muV]','FontSize',22,'FontName','Calibri');
% hold(gca,'off');
% title('Amplitude vs frequency vs oscillation width of atoms decomposing ripples','FontSize',26, 'FontName','Calibri');
% %view(gca,[7.4825840958656 29.478171859023]);
% %view(gca,[22.8748882338293 29.4900391472166]);
% %view(gca,[4.84058646938357 9.27052661539923]);
% view(gca,[12.4307209985316 13.6554697972174]);
% 
% 
% hold(gca,'off');
% set(gca,'FontName','Calibri','FontSize',22);
% legend1 = legend(gca,'show');
% set(legend1,'FontSize',22);
% 
% 
% figure(8);
% scatter3((abs(conditions{1,4}(6,:))*100)-0.5, abs(conditions{1,4}(5,:)), abs(conditions{1,4}(8,:)),30, faceColor1, 'filled', 'MarkerFaceAlpha', 0.5);
% hold on;
% scatter3((abs(conditions{1,6}(6,:))*100)+0.5, abs(conditions{1,6}(5,:)), abs(conditions{1,6}(8,:)),30, faceColor2, 'filled', 'MarkerFaceAlpha', 0.5);
% scatter3((abs(conditions{1,5}(6,:))*100), abs(conditions{1,5}(5,:)), abs(conditions{1,5}(8,:)),30, faceColor3, 'filled', 'MarkerFaceAlpha', 0.5);
% % plot3(abs(conditions{1,1}(6,:)), abs(conditions{1,1}(5,:)), abs(conditions{1,1}(8,:)), 'o','LineWidth',1);
% % hold on;
% % plot3(abs(conditions{1,3}(6,:)), abs(conditions{1,3}(5,:)), abs(conditions{1,3}(8,:)), 'o','LineWidth',1);
% % plot3(abs(conditions{1,2}(6,:)), abs(conditions{1,2}(5,:)), abs(conditions{1,2}(8,:)), 'o','LineWidth',1);
% hold off;
% legend('All EoI', 'AMP', 'REC', 'FontSize', 22, 'FontName','Calibri');
% set(gca,'FontName','Calibri','FontSize',22);
% ylabel('Frequency [Hz]','FontSize',22,'FontName','Calibri');
% xlabel('Oscillatory fraction of the atom [%]','FontSize',22,'FontName','Calibri');
% zlabel('Amplitude [\muV]','FontSize',22,'FontName','Calibri');
% hold(gca,'off');
% title('Amplitude vs frequency vs width of atoms decomposing fast ripples','FontSize',26, 'FontName','Calibri');
% %view(gca,[7.4825840958656 29.478171859023]);
% %view(gca,[22.8748882338293 29.4900391472166]);
% %view(gca,[4.84058646938357 9.27052661539923]);
% view(gca,[12.4307209985316 13.6554697972174]);
% 
% 
% hold(gca,'off');
% set(gca,'FontName','Calibri','FontSize',22);
% legend1 = legend(gca,'show');
% set(legend1,'FontSize',22);




