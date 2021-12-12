clear;close all;clc; % clear:ワークスペースの変数をクリア clc:コマンドウィンドウのクリア


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
fs2 = fs/2;          % ナイキスト周波数
N=40;                %ヒルベルト変換器の次数
L=15000;             %信号長 15000
ts=1/fs;             %サンプリング周期
t=(0:L-1)*ts;        %サンプリング時間
F=0:0.1:fs2;         %周波数設定
p_area1=0.05;        %遷移領域(f≒0)
p_area2=1-p_area1;   %遷移領域(f≒1)
sig=sin(2*pi*1800*t); %正弦波信号

%フィルタ係数の算出
[h,err,res]=firpm(N,[p_area1 p_area2],[1 1],{10000},'Hilbert');    
[H,W]=freqz(h,1,F,fs); 

sig_1=[sig zeros(1,N/2)]; %%遅延補正
sig_Im=filter(h,1,sig_1); %ヒルベルト変換
sig_Re=sig_1(1:end-N/2); %実部
sig_analysis=sig_Re+1i*sig_Im((N/2+1):end); %解析信号


f = [0 0.3 0.4 0.6 0.7 1];
a = [0 0 1 1 0 0];
bandpass_coefficient = firpm(N,f,a,{10000});
[Bandpass,w_bandpass] = freqz(bandpass_coefficient,1,F,fs);

%% 描画(バンドパス)
figure
plot(w_bandpass/fs2,20*log(abs(Bandpass)));
ylim([-100 20])
xlabel 'Normalized Frequency (\times\pi rad/sample)'
ylabel 'Magnitude   [dB]'
exportgraphics(gca,strcat('..\CAS_final\fig\normal_bandpass','.pdf'),'ContentType','vector');

%% 描画(ヒルベルト変換)

figure
plot(W/fs2,20*log(abs(H)));
ylim([-100 20])
xlabel 'Normalized Frequency (\times\pi rad/sample)'
ylabel 'Magnitude   [dB]'
exportgraphics(gca,strcat('..\CAS_final\fig\hilbert_amplitude','.pdf'),'ContentType','vector');
%%{
% figure を作成
figure1=figure;
% axes を作成
axes1 = axes('Parent',figure1);
hold(axes1,'on');
plot(W/fs2,20*log(abs(H)),'LineWidth',2,'Parent',axes1,"color",color2);
% plot(W/fs2,abs(H),'LineWidth',2,'Parent',axes1,"color",color2);
ylabel('Magnitude');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylim([-100 20]);
legend('N=42','N=60')
box(axes1,'on');
set(axes1,'FontSize',16,'FontName','Times New Roman');
% axes を作成
axes2 = axes('Parent',figure1,...
    'Position',[0.3 0.3 0.45 0.4]);
hold(axes2,'on');
% plot(W/fs2,abs(H),'LineWidth',2,'Parent',axes2,"color",color2);
plot(W/fs2,20*log(abs(H)),'LineWidth',2,'Parent',axes2,"color",color2);
axis(axes2, [0 0.3 0.5 1]);
box(axes2,'on');
ylim([-1 1]);
xlim([0.2 0.8]);
set(axes2,'FontSize',14);
%}}



figure;
plot(t,sig);
hold on;
plot(t,sig_Im((N/2+1):end));
xlim([0 0.003]);
ylim([-1.2 2]);
legend('Input','Hilbert(N=42)');
xlabel("Time[s]");
ylabel("Amplitude");
