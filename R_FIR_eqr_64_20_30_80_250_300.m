function Hd = R_FIR_eqr_64_20_30_80_250_300
%R_FIR_EQR_64_20_30_80_250_300 Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.13 and Signal Processing Toolbox 9.1.
% Generated on: 09-Oct-2023 16:15:34

% Equiripple Bandpass filter designed using the FIRPM function.

% All frequency values are in Hz.
Fs = 2000;  % Sampling Frequency

N      = 64;   % Order
Fstop1 = 30;   % First Stopband Frequency
Fpass1 = 80;   % First Passband Frequency
Fpass2 = 250;  % Second Passband Frequency
Fstop2 = 300;  % Second Stopband Frequency
Wstop1 = 1;    % First Stopband Weight
Wpass  = 1;    % Passband Weight
Wstop2 = 1;    % Second Stopband Weight
dens   = 20;   % Density Factor

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, [0 Fstop1 Fpass1 Fpass2 Fstop2 Fs/2]/(Fs/2), [0 0 1 1 0 ...
           0], [Wstop1 Wpass Wstop2], {dens});
Hd = dfilt.dffir(b);

% [EOF]
