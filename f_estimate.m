clear;close all;clc

color1=[0 0.4470 0.7410];
color2=[0.8500 0.3250 0.0980];
color3=[0.9290 0.6940 0.1250];
color4=[0.4940 0.1840 0.5560];
color5=[0.4660 0.6740 0.1880];
color6=[0.3010 0.7450 0.9330];

set(0,'defaultTextFontName','Times New Roman');
set(0,'defaultLegendFontName','Times New Roman');
set(0,'defaultAxesFontName','times new roman');
set(0,'defaultAxesGridLineStyle',':');
set(0,'defaultAxesBox','on');
set(0,'defaultAxesYGrid','on');
set(0,'defaultAxesXGrid','on');
set(0,'defaultAxesFontSize',16);
set(0, 'DefaultAxesLineWidth', 1, 'DefaultLineLineWidth', 2); 

%% %% params
% filter params
N = 28;
f_left_stop  = 0.2;
f_left_pass  = 0.3;
f_right_pass = 0.7;
f_right_stop = 0.8;
notch_freq_1 = 0.07;
notch_freq_2 = 0.15;
Weight = [1 1 1];
fs = 4000;
fs2 = fs/2;
L=15000;
ts=1/fs;
t=(0:L-1)*ts;
F=0:0.1:fs2;
f = [0 f_left_stop f_left_pass f_right_pass f_right_stop 1];
a = [0 0 1 1 0 0];

y_low = -200;
y_high = 20;

%% design filter

% notch filter
nf1 = 3600*180/fs2;    r1=1;    D2=1;
nf2 = 7200*180/fs2;    r2=1;    D4=1;
nf3 = 3600*180/fs2;    r3=1;

notch_weight1=[1 1 1];
notch_amplitude1=[0 0 1 1 0 0];

[notch1,notch2,notch3]=notch(nf1,nf2,nf3,D2,D4,fs2,r1,r2,r3);
notch1 = [1 -2*cos(notch_freq_1*pi) 1*1];
notch2 = [1 -2*cos(notch_freq_2*pi) 1*1];
notch3 = conv(notch1, notch2);

notch1 = notch1/sum(notch1);
notch2 = notch2/sum(notch2);
notch3 = notch3/sum(notch3);
save notch1.cof notch1 -ascii;
save notch2.cof notch2 -ascii;
save notch3.cof notch3 -ascii;
fresp1 = {'fresp_h1',a};
fresp2 = {'fresp_h2',a};
fresp3 = {'fresp_h3',a};
h_notch = load('notch3.cof');
[H_notch,w]=freqz(h_notch,1,F,fs);   

% hilbert filter 
h_HT=firpm(N,f,fresp3,Weight,'hilbert');
[H_HT,w]=freqz(h_HT,1,F,fs);

% convolution
h_CHT=conv(h_notch,h_HT);
[H_CHT,w]=freqz(h_CHT,1,F,fs);

%% signal process
% In this section, a signal with params below input into the filter designed
% above, and fetch output a signals.
sig_f = 1000;
sig_f_n = 140;
A = 1;
A_noize_rate = 0.7;

sig = A*cos(2*pi*sig_f*t);
sig_n = sig + A_noize_rate*A*cos(2*pi*sig_f_n*t);

sig_n_zeros = [sig_n zeros(1,N/2)];
sig_n_Im = filter(h_CHT,1,sig_n_zeros);
sig_n_Re = sig_n_zeros(1:end-N/2);
sig_analysis = sig_n_Re + 1i*sig_n_Im((N/2+1):end);


figure;
plot(t, sig_n);
hold on
plot(t, sig_n_Im((N/2+1):end));
hold off
xlim([0.01 0.02]);
ylim([-1.5 1.5])
legend()

%% estimate frequency
phase=angle(sig_analysis); %解析信号の位相角を計算
phase=unwrap(phase); %滑らかになるよう位相角の修正
omega_tilde_low=diff(phase)/ts; %位相角の微分
IF=omega_tilde_low/(2*pi); %瞬時周波数
plotfft(IF,fs); %瞬時周波数のスペクトルを表示
