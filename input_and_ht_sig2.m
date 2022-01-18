clear;close all;clc; % clear:ワークスペースの変数をクリア clc:コマンドウィンドウのクリア
%% setting default params
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
set(0,'defaultAxesYGrid','on');
set(0,'defaultAxesXGrid','on');
set(0,'defaultAxesFontSize',16);
set(0, 'DefaultAxesLineWidth', 1, 'DefaultLineLineWidth', 2); 

%% setting params
fs = 50000;
fs2 = fs/2;
N=30;
L=15000;
ts=1/fs;
t=(0:L-1)*ts;
F=0:0.1:fs2;
p_area1=0.05;
p_area2=1-p_area1;
sig=sin(2*pi*1800*t);

%% design filter
[h,err,res]=firpm(N,[p_area1 p_area2],[1 1],{10000},'Hilbert');    
[H,W]=freqz(h,1,F,fs); 

%% draw filter
figure
plot(W/fs2,20*log(abs(H)));
ylim([-100 20])
xlabel 'Normalized Frequency (\times\pi rad/sample)'
ylabel 'Magnitude   [dB]'
exportgraphics(gca,strcat('.\figure\amp_ht_N=',num2str(N),'.pdf'),'ContentType','vector');

%% filtering
sig_1=[sig zeros(1,N/2)];
sig_Im=filter(h,1,sig_1);
sig_Re=sig_1(1:end-N/2);
sig_analysis=sig_Re+1i*sig_Im((N/2+1):end);

%% draw in/out signals
figure;
plot(t,sig);
hold on;
plot(t,sig_Im((N/2+1):end));
xlim([0 0.003]);
ylim([-1.2 2]);
legend('Input','Hilbert(N=42)');
xlabel("Time[s]");
ylabel("Amplitude");
