function [Rfiltered, FRfiltered] = filter_HFOBands(rawdata, fs)
% Takes as input the raw data (5 minutes) in bipolar montage (23 channels)

%cd('\\fs-home\ulrta$\Documents\Master_Biomed\Code\DataSet\Filter');
load('Gotman_FR_coef.mat');
load('Gotman_R_coef.mat');
load('At40_ATB10_FR_coef.mat');
load('At40_ATB10_R_coef.mat');
load('At60_ATB20_FR_coef.mat');
load('At60_ATB20_R_coef.mat');

% Filter Reference
RfilteredData_G = filtfilt(Gotman_R, 1, rawdata(:,:)')'; 
FRfilteredData_G = filtfilt(Gotman_FR, 1, rawdata(:,:)')'; 

% Filter At40 ATB10, Version1
RfilteredData_40_10 = filtfilt(At40_ATB10_R, 1, rawdata(:,:)')'; 
FRfilteredData_40_10 = filtfilt(At40_ATB10_FR, 1, rawdata(:,:)')'; 

% Filter At60 ATB20, Version2
RfilteredData_60_20 = filtfilt(At60_ATB20_R, 1, rawdata(:,:)')'; 
FRfilteredData_60_20 = filtfilt(At60_ATB20_FR, 1, rawdata(:,:)')'; 

Rfiltered = {RfilteredData_G RfilteredData_40_10 RfilteredData_60_20};
FRfiltered = {FRfilteredData_G FRfilteredData_40_10 FRfilteredData_60_20};

end