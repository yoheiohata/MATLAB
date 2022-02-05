clear;close all;clc

set(0,'defaultTextFontName','Times New Roman');
set(0,'defaultLegendFontName','Times New Roman');
set(0,'defaultAxesFontName','times new roman');
set(0,'defaultAxesGridLineStyle',':');
set(0,'defaultAxesBox','on');
set(0,'defaultAxesYGrid','on');
set(0,'defaultAxesXGrid','on');
set(0,'defaultAxesFontSize',16);
set(0,'DefaultAxesLineWidth', 1, 'DefaultLineLineWidth', 2); 

%% params
% filter params
N = 30;
f_left  = 0.05;
f_right = 0.95;
Weight = [1 1];
fs = 4000;
fs2 = fs/2;
L=15000;
ts=1/fs;
t=(0:L-1)*ts;
F=0:0.1:fs2;
f = [f_left f_right];
a = [1 1];

y_low = -150;
y_high = 20;

%% design filter

% conventional
[h_HT,err,res]=firpm(N,f,a,{10000}, "Hilbert");
[H_HT,w]=freqz(h_HT,1,F,fs);
%% draw
figure
hold on
plot(w/fs2, 20*log(abs(H_HT)),"DisplayName","N="+N);
hold off
xlim([0 1]);
ylim([y_low y_high])
xlabel("Normalized Angular Frequency(\times\pi)   [rad/sample]");
ylabel("Magnitude   [dB] ");
legend("location", "northeast")
exportgraphics(gca,strcat('.\figure\amp_ht.pdf'),'ContentType','vector');

%% SNR

% signal params
A = 1;
A_noize_rate = 0.7;
sig_f = 0.46*fs2;
sig_f_n = 0.1*fs2;

% define signal
sig = A*cos(2*pi*sig_f*t);
sig_n = sig + A_noize_rate*A*cos(2*pi*sig_f_n*t);

% filtering
sig_n_zeros = [sig_n zeros(1,N/2)];
sig_n_Re = sig_n_zeros(1:end-N/2);
sig_n_Im = filter(h_HT,1,sig_n_zeros);

% cal SNrate
SN = snr(sig_n_Im, fs)