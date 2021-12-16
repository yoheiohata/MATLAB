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
N1 = 60;
f_left_stop  = 0.19;
f_left_pass  = 0.2;
f_right_pass = 0.8;
f_right_stop = 0.81;
Weight = [1 1 1];
fs = 4000;
fs2 = fs/2;
L=15000;
ts=1/fs;
t=(0:L-1)*ts;
F=0:0.1:fs2;
f = [0 f_left_stop f_left_pass f_right_pass f_right_stop 1];
a = [0 0 1 1 0 0];

y_low = -100;
y_high = 20;

%% design filter

% conventional
[h_HT,err,res]=firpm(N,f,a,{10000},Weight);
[H_HT,w]=freqz(h_HT,1,F,fs);

[h_HT1,err,res]=firpm(N1,f,a,{10000},Weight);
[H_HT1,w]=freqz(h_HT1,1,F,fs);

%% draw
figure
hold on
plot(w/fs2, 20*log(abs(H_HT)),"DisplayName","N="+N);
plot(w/fs2, 20*log(abs(H_HT1)),"DisplayName","N="+N1);
hold off
xlim([0 1]);
ylim([y_low y_high])
legend("location", "northeast")
exportgraphics(gca,strcat('.\figure\amp_bpht.pdf'),'ContentType','vector');

