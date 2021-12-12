clear;close all;clc
%**********************************************************************

%---------   SN比計算   -----------

%**********************************************************************
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
set(0,'DefaultAxesLineWidth', 1, 'DefaultLineLineWidth', 2); 

disp('------------------------------------------------------')
disp('       f                     SNR')
disp('------------------------------|-----------------------')
disp('                   HT         |       nCHT')
disp('------------------------------|-----------------------')
disp('            normal     noize  | normal    noize')
disp('==============================|=======================')
%% フィルタ設計部

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%---------------------フィルタのパラメータ変更部分------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%各フィルタの次数部
N_rearHT=26;

%正規化周波数帯域設定
f_left_stop  = 0.2;
f_right_stop = 0.8;
f_left_pass  = 0.3;
f_right_pass = 0.7;



%ノッチの場所
notch_freq_1 = 0.1;       
notch_freq_2 = 0.15;



%重み[左側阻止域　通過域　右側阻止域]
Weight = [1 1 1];         


%%--------------------------------------------------------------------

%ヒルベルトフィルタ設計部
fs = 50000;          %Sampling frequency;
fs2 = fs/2;          %ナイキスト周波数
L=15000;             %信号長 15000
ts=1/fs;             %サンプリング周期
t=(0:L-1)*ts;        %サンプリング時間
F=0:0.1:fs2;         %周波数設定
f = [0 f_left_stop f_left_pass f_right_pass f_right_stop 1]; %振幅特性の通過域、遷移域、阻止域の周波数指定
a = [0 0 1 1 0 0];         %理想周波数特性[阻止域　通過域　阻止域]
%振幅特性0、振幅特性1、振幅特性0、
%繋がってないところは遷移域


%% 信号定義  

%振幅定義
A=1;
A_noize_rate_1=0.7;
A_noize_rate_2=0.5;

% インクリメント
p=0;
q=0;
r=0;
% figureの設定
y_low = -100;
y_high = 20;
%% ノッチヒルベルト用ヒルベルト変換器(手法1と3に用いるやつ)
%ノッチヒルベルト用ヒルベルト変換器
[h_HT,err,res]=firpm(N_rearHT, f, a, {10000}, Weight,'hilbert'); %フィルタ係数の算出
[H_HT,w]=freqz(h_HT,1,F,fs);
%% ノッチフィルタ

%　ノッチフィルタ
nf1 = 3600*180/fs2;    r1=1;        %notch1 & on,off
nf2 = 7200*180/fs2;    r2=1;        %notch2 & on,off
nf3 = 3600*180/fs2;    r3=1;        %notch3 & on,off

D2=1;
D4=1;

notch_weight1=[1 1 1];     %伝達関数後段の重み[第一帯域 第二帯域　第三帯域]
notch_amplitude1=[0 0 1 1 0 0];

[notch1,notch2,notch3]=notch(nf1,nf2,nf3,D2,D4,fs2,r1,r2,r3);
notch1 = [1 -2*cos(notch_freq_1*pi) 1*1];
notch2 = [1 -2*cos(notch_freq_2*pi) 1*1];
notch3 = conv(notch1,notch2);

notch1 = notch1/sum(notch1);
notch2 = notch2/sum(notch2);
notch3 = notch3/sum(notch3);

save notch1.cof notch1 -ascii;
save notch2.cof notch2 -ascii;
save notch3.cof notch3 -ascii;


h_coefficient_notch = load('notch3.cof');
[h_amplitude_notch,w]=freqz(h_coefficient_notch,1,F,fs);   
fresp3 = {'fresp_h3',a};
%% フィルタ合成(ノッチフィルタ＋補正ヒルベルト変換)


%補正ヒルベルトフィルタの作成
h_coefficient_compensationHT=firpm(N_rearHT,f,fresp3,Weight,'hilbert');
[h_amplitude_compensationCHT,w]=freqz(h_coefficient_compensationHT,1,F,fs);               %補正のフィルタの振幅特性


%フィルタ合成(ノッチフィルタ＋ヒルベルト変換器)
h_coefficient_notchCHT=conv(h_coefficient_notch,h_coefficient_compensationHT);
[h_amplitude_notchCHT,w]=freqz(h_coefficient_notchCHT,1,F,fs);
h_amplitude_nCHT_evaluation(:,1) = w/fs2;
h_amplitude_nCHT_evaluation(:,2) = h_amplitude_notchCHT;

%% 設計したフィルタでヒルベルト変換
for j=(f_left_pass*fs2):2000:(f_right_pass*fs2)
    p=p+1;
    
%% 信号定義
%周波数定義
sig_f=j;                        %周波数
normalized_sig_f = sig_f/fs2;   %正規化周波数
sig_f_noize1=notch_freq_1*fs2;
sig_f_noize2=notch_freq_2*fs2;

%正弦波信号
sig=A*cos(2*pi*sig_f*t);
%ノイズあり信号

sig_n=A*cos(2*pi*sig_f*t)+A_noize_rate_1*A*cos(2*pi*sig_f_noize1*t); 
% sig_n=A*cos(2*pi*sig_f*t)+A_noize_rate_2*A*cos(2*pi*sig_f_noize2*t);
% sig_n=A*cos(2*pi*sig_f*t)+A_noize_rate_1*A*cos(2*pi*sig_f_noize1*t)+...
%       A_noize_rate_2*A*cos(2*pi*sig_f_noize2*t); 

%ゼロ値付加
sig_1=[sig zeros(1,N_rearHT/2)];
sig_n_1=[sig_n zeros(1,N_rearHT/2)];

%% ヒルベルト変換とSN比計算

%阻止域を有するのFIRヒルベルト変換器とSN比
sig_htout=filter(h_HT,1,sig_1);
sig_noize_htout=filter(h_HT,1,sig_n_1);
SNR_out=snr(sig_htout);
SNR_noize_out=snr(sig_noize_htout);

%nCHTとそのSN比
sig_Im_nCHT=filter(h_coefficient_notchCHT,1,sig_1);
sig_noize_Im_nCHT=filter(h_coefficient_notchCHT,1,sig_n_1);
SNR_nCHT_out=snr(sig_Im_nCHT);
SNR_noize_nCHT_out=snr(sig_noize_Im_nCHT);

snr_outcome(p,1) = normalized_sig_f;
snr_outcome(p,2) = SNR_out;
snr_outcome(p,3) = SNR_noize_out;
snr_outcome(p,4) = SNR_nCHT_out;
snr_outcome(p,5) = SNR_noize_nCHT_out;
%% 入出力描画
% figure;
% plot(t,sig,"DisplayName","in");
% hold on;
% plot(t,sig_Im_nCHT((N_rearHT/2+1):end),"DisplayName","nCHT");
% plot(t,sig_noize_Im_nCHT((N_rearHT/2+1):end),"DisplayName","nCHT_w_n");
% plot(t,sig_htout((N_rearHT/2+1):end),"DisplayName","HT");
% plot(t,sig_noize_htout((N_rearHT/2+1):end),"DisplayName","HT_w_n");
% hold off;
% xlim([0 0.0005]);
% ylim([-1 1]);
% legend();
% xlabel("time   [s]");
% ylabel("Amplitude");
% % exportgraphics(gca,strcat('..\CAS_final\fig\SNR_normalized\f=',num2str(normalized_sig_f),'.pdf'),'ContentType','vector');
% % exportgraphics(gca,strcat('..\CAS_final\fig\SNR\f=',num2str(sig_f),'.pdf'),'ContentType','vector');
% % 
end
disp(snr_outcome)
disp('------------------------------------------------------')
SNR_outcome(:,1)=mean(snr_outcome(:,2),1);
SNR_outcome(:,2)=mean(snr_outcome(:,4),1);
SNR_outcome_noize(:,1)=mean(snr_outcome(:,3),1);
SNR_outcome_noize(:,2)=mean(snr_outcome(:,5),1);
disp(SNR_outcome)
disp(SNR_outcome_noize)


