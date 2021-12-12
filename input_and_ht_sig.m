clear;close all;clc;

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


fs = 50000;          % Sampling frequency;
fs2 = fs/2;
N=42; %ヒルベルト変換器の次数
N_high=60; %ヒルベルト変換器の次数(低いやつ)

N_Low=26;
N_High=36;
L=15000*1; %信号長 15000
ts=1/fs; %サンプリング周期
t=(0:L-1)*ts; %サンプリング時間
F=0:0.1:fs2;     %　周波数設定
p_area1=0.05; %遷移領域(f≒0)
p_area2=1-p_area1; %遷移領域(f≒1)
sig_f=1800; %正弦波信号の周波数
sig=sin(2*pi*sig_f*t); %正弦波信号


[h,err,res]=firpm(N,[p_area1 p_area2],[1 1],{10000},'Hilbert'); %フィルタ係数の算出
[H,W]=freqz(h,1,F,fs);
[h60,err60,res60]=firpm(N_high,[p_area1 p_area2],[1 1],{10000},'Hilbert'); %フィルタ係数の算出
[H60,W60]=freqz(h60,1,F,fs);


%%{
% figure を作成
figure1=figure;

% axes を作成
axes1 = axes('Parent',figure1);
hold(axes1,'on');
plot(W/fs2,abs(H),'LineWidth',2,'Parent',axes1,"color",color2);
plot(W60/fs2,abs(H60),'LineWidth',2,'Parent',axes1,"color",color5);
ylabel('Magnitude');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylim([0 1.05]);
legend('N=42','N=60')
box(axes1,'on');
set(axes1,'FontSize',16,'FontName','Times New Roman');

% axes を作成
axes2 = axes('Parent',figure1,...
    'Position',[0.3 0.3 0.45 0.4]);
hold(axes2,'on');
% plot(W/fs2,abs(H),W60/fs2,abs(H60),'LineWidth',2,'Parent',axes2,"color",color2,"color",color5);
plot(W/fs2,abs(H),'LineWidth',2,'Parent',axes2,"color",color2);
plot(W60/fs2,abs(H60),'LineWidth',2,'Parent',axes2,"color",color5);
plot(W/fs2,20*log(abs(H)),'LineWidth',2,'Parent',axes2,"color",color2);
plot(W60/fs2,20*log(abs(H60)),'LineWidth',2,'Parent',axes2,"color",color5);


axis(axes2, [0 0.3 0.5 1]);
box(axes2,'on');
ylim([0.98 1.02]);
xlim([0.04 0.1]);
set(axes2,'FontSize',14);
%}


figure

hold on
plot(W/fs2,20*log(abs(H)),"DisplayName","N="+N_Low);
plot(W60/fs2,20*log(abs(H60)),"DisplayName","N="+N_High);
hold off

xlabel("Normalized Angular Frequency(\times\pi)　　　(rad/sample)");
ylabel("Magnitude       [dB] ");
legend()
xlim([0 1]);
ylim([-100 20]);


%%遅延補正
sig_1=[sig zeros(1,N/2)]; 
%%%%
sig_Im=filter(h,1,sig_1); %ヒルベルト変換
sig_Re=sig_1(1:end-N/2); %実部
sig_analysis=sig_Re+1i*sig_Im((N/2+1):end); %解析信号

%%遅延補正
sig_high=[sig zeros(1,N_high/2)]; 
%%%%
sig_Im_high=filter(h60,1,sig_high); %ヒルベルト変換
sig_Re_high=sig_high(1:end-N_high/2); %実部
sig_analysis_high=sig_Re_high+1i*sig_Im_high((N_high/2+1):end); %解析信号


figure;
plot(t,sig);
hold on;
plot(t,sig_Im((N/2+1):end));
plot(t,sig_Im_high((N_high/2+1):end));
xlim([0 0.003]);
ylim([-1.2 2]);
legend('Input','Hilbert(N=42)','Hilbert(N=60)');
xlabel("Time[s]");
ylabel("Amplitude");

figure;
plot(t,sig);
hold on;
plot(t,sig_Im((N/2+1):end));
plot(t,sig_Im_high((N_high/2+1):end));
xlim([0.0014 0.003]);
ylim([0.98 1.04]);
legend('Input','Hilbert(N=42)','Hilbert(N=60)');
xlabel("Time[s]");
ylabel("Amplitude");

rate=10;

sig_up=interp(sig,rate);
sig_Im_up=interp(sig_Im,rate);
sig_Im_high_up=interp(sig_Im_high,rate);
t_up=interp(t,rate);

figure;
plot(t_up,sig_up);
xlim([0 0.003]);
ylim([-1.2 2]);
legend('Input');
xlabel("Time[s]");
ylabel("Amplitude");

figure;
plot(t_up,sig_up);
xlim([0.0016 0.003]);
ylim([0.98 1.025]);
legend('Input');
xlabel("Time[s]");
ylabel("Amplitude");


figure;
plot(t_up,sig_up);
hold on;
plot(t_up,sig_Im_up((rate*N/2+1):end));
plot(t_up,sig_Im_high_up((rate*N_high/2+1):end),"-.",'color',color5);
xlim([0 0.003]);
ylim([-1.2 2]);
legend('Input','Hilbert(N=42)','Hilbert(N=60)');
xlabel("Time[s]");
ylabel("Amplitude");

figure;
plot(t_up,sig_up);
hold on;
plot(t_up,sig_Im_up((rate*N/2+1):end));
plot(t_up,sig_Im_high_up((rate*N_high/2+1):end),'color',color5);
xlim([0.0016 0.003]);
ylim([0.98 1.025]);
legend('Input','Hilbert(N=42)','Hilbert(N=60)');
xlabel("Time[s]");
ylabel("Amplitude");

figure;
plot(t_up,sig_up);
hold on;
plot(t_up,sig_Im_up((rate*N/2+1):end));
plot(t_up,sig_Im_high_up((rate*N_high/2+1):end),'color',color3);
xlim([0.0016 0.003]);
ylim([0.98 1.025]);
legend('Input','Hilbert(N=42)','Hilbert(N=60)');
xlabel("Time[s]");
ylabel("Amplitude")