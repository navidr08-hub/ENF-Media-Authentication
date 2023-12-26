%% Step 2: ENF Signal Harmonics

% @60Hz

% Reading file data and setting inputs
[recording, Fs] = audioread('audio_files\recording.wav');
BlockSize = Fs*16;
ZeroPad = 0;
Overlap = 0.5;
Window = hanning(BlockSize, 'periodic');
Frequency = 60;

% Calling enf function and generating plot
[y1, ~, ~] = enf(recording, Fs, BlockSize, ZeroPad, Overlap, Window, Frequency);
[r,c] = size(y1);
k1 = Frequency - 1; k2 = Frequency + 1;
increment = (k2 - k1)/c;
x = k1:increment:(k2-increment);
y = 1:r;

surf(x,y,y1); title("ENF Magnitude Spectrum @60 Hz of 'recording.wav'"); ylabel('Block Number'); xlabel('Frequency (Hz)');

% @120Hz

% Reading file data and setting inputs
[recording, Fs] = audioread('audio_files\recording.wav');
BlockSize = Fs*16;
ZeroPad = 0;
Overlap = 0.5;
Window = hanning(BlockSize, 'periodic');
Frequency = 120;

% Calling enf function and generating plot
[y1, ~, ~] = enf(recording, Fs, BlockSize, ZeroPad, Overlap, Window, Frequency);
[r,c] = size(y1);
k1 = Frequency - 1; k2 = Frequency + 1;
increment = (k2 - k1)/c;
x = k1:increment:(k2-increment);
y = 1:r;

figure;
surf(x,y,y1); title("ENF Magnitude Spectrum @120 Hz of 'recording.wav'"); ylabel('Block Number'); xlabel('Frequency (Hz)');

% @180Hz

% Reading file data and setting inputs
[recording, Fs] = audioread('audio_files\recording.wav');
BlockSize = Fs*16;
ZeroPad = 0;
Overlap = 0.5;
Window = hanning(BlockSize, 'periodic');
Frequency = 180;

% Calling enf function and generating plot
[y1, ~, ~] = enf(recording, Fs, BlockSize, ZeroPad, Overlap, Window, Frequency);
[r,c] = size(y1);
k1 = Frequency - 1; k2 = Frequency + 1;
increment = (k2 - k1)/c;
x = k1:increment:(k2-increment);
y = 1:r;

figure;
surf(x,y,y1); title("ENF Magnitude Spectrum @180 Hz of 'recording.wav'"); ylabel('Block Number'); xlabel('Frequency (Hz)');

% @240Hz

% Reading file data and setting inputs
[recording, Fs] = audioread('audio_files\recording.wav');
BlockSize = Fs*16;
ZeroPad = 0;
Overlap = 0.5;
Window = hanning(BlockSize, 'periodic');
Frequency = 240;

% Calling enf function and generating plot
[y1, ~, ~] = enf(recording, Fs, BlockSize, ZeroPad, Overlap, Window, Frequency);
[r,c] = size(y1);
k1 = Frequency - 1; k2 = Frequency + 1;
increment = (k2 - k1)/c;
x = k1:increment:(k2-increment);
y = 1:r;

figure;
surf(x,y,y1); title("ENF Magnitude Spectrum @240 Hz of 'recording.wav'"); ylabel('Block Number'); xlabel('Frequency (Hz)');

%% Step 3: Comparison of 2 ENF Signals without preprocessing

% Reading file data and setting inputs
[recording, Fs] = audioread('audio_files\recording.wav');
[ground, ~] = audioread('audio_files\ground truth.wav');
BlockSize = Fs*16;
ZeroPad = 0;
Overlap = 0.5;
Window = hanning(BlockSize, 'periodic');
Frequency = 120;

k1 = Frequency - 1; k2 = Frequency + 1;

% Calling enf function and generating plot for reference signal
[g1, g2, g3] = enf(ground, Fs, BlockSize, ZeroPad, Overlap, Window, Frequency);
[g_row, g_col] = size(g1);

g_incr = (k2 - k1)/g_col;
x_g = k1:g_incr:(k2-g_incr);
y_g = 1:g_row;

surf(x_g,y_g,g1); title("Reference Signal ENF 'ground truth.wav'"); ylabel('Block Number'); xlabel('Frequency (Hz)');

% Calling enf function and generating plot for sample signal
[r1, r2, r3] = enf(recording, Fs, BlockSize, ZeroPad, Overlap, Window, Frequency);
[r_row, r_col] = size(r1);

r_incr = (k2 - k1)/r_col;
x_r = k1:r_incr:(k2-r_incr);
y_r = 1:r_row;

figure;
surf(x_r,y_r,r1); title("Recording Signal ENF 'recording.wav'"); ylabel('Block Number'); xlabel('Frequency (Hz)');

% Subtracting mean from weighted energies
g_mean = mean(g3);
r_mean = mean(r3);
ENF_rec = r3 - r_mean;
ENF_gt = g3 - g_mean;
[c, lags] = xcorr(ENF_gt, ENF_rec);

% Plot Cross-Correlation Plot
figure;
plot(lags, c);
title('Time lag between Recording and Reference ENF singals using Cross-Correlation Plot'); xlabel('Lags'); ylabel('Cross-Correlation');

% Zero-padding to align the sample signal by lag of 2 windows
% aligned_r2 = padarray(r2',2,r_mean);

% % Plotting Max Magnitudes of ENF pre-alignment
% figure;
% plot(g2);
% hold on
% plot(r2);
% title('Max Magnitude of Recording vs Reference ENF signals (pre-alignment)'); xlabel('Block Number'); ylabel('Magnitude');
% 
% % Plotting Weighted Magnitudes post-alignment
% figure;
% plot(g2);
% hold on;
% plot(aligned_r2);
% title('Max Magnitude of Recording vs Reference ENF signals (post-alignment)'); xlabel('Block Number'); ylabel('Magnitude');

% Plotting Weighted Magnitudes pre-alignment
figure;
plot(g3); 
hold on; 
plot(r3); 
title('Weighted Magnitude of Recording vs Reference ENF signals (pre-alignment)'); xlabel('Block Number'); ylabel('Magnitude');

% Zero-padding to align the sample signal by lag of 2 windows
aligned_r3 = padarray(r3',2,r_mean);

% Plotting Weighted Magnitudes post-alignment
figure;
plot(g3);
hold on;
plot(aligned_r3);
title('Weighted Magnitude of Recording vs Reference ENF signals (post-alignment)'); xlabel('Block Number'); ylabel('Magnitude');

%% Step 4: Comparison of 2 ENF Signals with pre-processing

[recording, ~] = audioread('audio_files\recording.wav');
[ground, ~] = audioread('audio_files\ground truth.wav');

load('SOS.mat');
load('G.mat');
LPF_rec = downsample(filtfilt(SOS, G, recording), 100);
LPF_gt = downsample(filtfilt(SOS, G, ground), 100);

Fs = 441;
BlockSize = Fs*16;
ZeroPad = 16384 - (Fs*16);
Overlap = 0.5;
Window = hanning(BlockSize, 'periodic');
Frequency = 120; 

k1 = Frequency - 1; k2 = Frequency + 1;

% Calling enf function and generating plot for reference signal
[g1, ~, g3] = enf(LPF_gt, Fs, BlockSize, ZeroPad, Overlap, Window, Frequency);
[g_row, g_col] = size(g1);

g_incr = (k2 - k1)/g_col;
x_g = k1:g_incr:(k2-g_incr);
y_g = 1:g_row;

surf(x_g,y_g,g1); title("Reference Signal ENF 'ground truth.wav' after pre-processing"); ylabel('Block Number'); xlabel('Frequency (Hz)');

% Calling enf function and generating plot for sample signal
[r1, ~, r3] = enf(LPF_rec, Fs, BlockSize, ZeroPad, Overlap, Window, Frequency);
[r_row, r_col] = size(r1);

r_incr = (k2 - k1)/r_col;
x_r = k1:r_incr:(k2-r_incr);
y_r = 1:r_row;

figure;
surf(x_r,y_r,r1); title("Recording Signal ENF 'recording.wav' after pre-processing"); ylabel('Block Number'); xlabel('Frequency (Hz)');

% Subtracting mean from weighted energies
g_mean = mean(g3);
r_mean = mean(r3);
ENF_rec = r3 - r_mean;
ENF_gt = g3 - g_mean;
[c, lags] = xcorr(ENF_gt, ENF_rec);

% Plot Cross-Correlation Plot
figure;
plot(lags, c);
title('Time lag between Recording and Reference ENF singals using Cross-Correlation Plot'); xlabel('Lags'); ylabel('Cross-Correlation');

% Plotting Weighted Magnitudes pre-alignment
figure;
plot(g3); 
hold on; 
plot(r3); 
title('Weighted Magnitude of Recording vs Reference ENF signals (pre-alignment)'); xlabel('Block Number'); ylabel('Magnitude');

% Zero-padding to align the sample signal by lag of 2 windows
aligned_r3 = padarray(r3',2,r_mean);

% Plotting Weighted Magnitudes post-alignment
figure;
plot(g3);
hold on;
plot(aligned_r3);
title('Weighted Magnitude of Recording vs Reference ENF signals (post-alignment)'); xlabel('Block Number'); ylabel('Magnitude');

%% Step 5: Comparison 2nd set of signals with pre-processing

[recording, ~] = audioread('audio_files\recording 2.wav');
[ground, ~] = audioread('audio_files\ground truth 2.wav');

load('SOS.mat');
load('G.mat');
LPF_rec = downsample(filtfilt(SOS, G, recording), 100);
LPF_gt = downsample(filtfilt(SOS, G, ground), 100);

Fs = 441;
BlockSize = Fs*16;
ZeroPad = 16384 - (Fs*16);
Overlap = 0.5;
Window = hanning(BlockSize, 'periodic');
Frequency = 120; 

k1 = Frequency - 1; k2 = Frequency + 1;

% Calling enf function and generating plot for reference signal
[g1, ~, g3] = enf(LPF_gt, Fs, BlockSize, ZeroPad, Overlap, Window, Frequency);
[g_row, g_col] = size(g1);

g_incr = (k2 - k1)/g_col;
x_g = k1:g_incr:(k2-g_incr);
y_g = 1:g_row;

surf(x_g,y_g,g1); title("Reference Signal ENF 'ground truth 2.wav' after pre-processing"); ylabel('Block Number'); xlabel('Frequency (Hz)');

% Calling enf function and generating plot for sample signal
[r1, ~, r3] = enf(LPF_rec, Fs, BlockSize, ZeroPad, Overlap, Window, Frequency);
[r_row, r_col] = size(r1);

r_incr = (k2 - k1)/r_col;
x_r = k1:r_incr:(k2-r_incr);
y_r = 1:r_row;

figure;
surf(x_r,y_r,r1); title("Recording Signal ENF 'recording 2.wav' after pre-processing"); ylabel('Block Number'); xlabel('Frequency (Hz)');

% Subtracting mean from weighted energies
g_mean = mean(g3);
r_mean = mean(r3);
ENF_rec = r3 - r_mean;
ENF_gt = g3 - g_mean;
[c, lags] = xcorr(ENF_gt, ENF_rec);

% Plot Cross-Correlation Plot
figure;
plot(lags, c);
title('Time lag between Recording and Reference ENF singals using Cross-Correlation Plot'); xlabel('Lags'); ylabel('Cross-Correlation');

% Plotting Weighted Magnitudes pre-alignment
figure;
plot(g3); 
hold on; 
plot(r3); 
title('Weighted Magnitude of Recording vs Reference ENF signals (pre-alignment)'); xlabel('Block Number'); ylabel('Magnitude');

% Zero-padding to align the sample signal by lag of 2 windows
aligned_r3 = padarray(r3',2,r_mean);

% Plotting Weighted Magnitudes post-alignment
figure;
plot(g3);
hold on;
plot(aligned_r3);
title('Weighted Magnitude of Recording vs Reference ENF signals (post-alignment)'); xlabel('Block Number'); ylabel('Magnitude');
