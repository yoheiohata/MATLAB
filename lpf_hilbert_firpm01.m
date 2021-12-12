clear;close all;clc
%%

%　ノッチフィルタが一個の場合

%ヒルベルトフィルタ設計部
fs = 50000;          %Sampling frequency;
fs2 = fs/2;          %ナイキスト周波数
N_HT=40;             %ヒルベルト変換器の次数
L=15000;             %信号長 15000
ts=1/fs;             %サンプリング周期
t=(0:L-1)*ts;        %サンプリング時間
F=0:0.1:fs2;         %周波数設定
f = [0 0.1 0.2 0.8 0.9 1]; %振幅特性の通過域、遷移域、阻止域の周波数指定
a = [0 0 1 1 0 0];         %理想周波数特性[阻止域　通過域　阻止域]
%振幅特性0、振幅特性1、振幅特性1、
%繋がってないところは遷移域
Weight = [3 30 3];         %重み[左側阻止域　通過域　右側阻止域]

%%
%ヒルベルト変換器

[h_HT,err,res]=firpm(N_HT, f,a,{10000},Weight,'Hilbert'); %フィルタ係数の算出
[H_HT,W]=freqz(h_HT,1,F,fs);

%%
%　ノッチフィルタ

N_rear=30;                 %フィルタ次数
notch_freq_1 = 0.1;        %ノッチの場所
notch_weight1=[1 1 1];     %伝達関数後段の重み[第一帯域 第二帯域　第三帯域]

notch1 = [1 -2*cos(notch_freq_1*pi) 1*1];
notch1 = notch1/sum(notch1);
save notch1.cof notch1 -ascii;

%% フィルタ合成

%ノッチフィルタ
h_coeff_notch = load('notch1.cof');
[h_notch,w]=freqz(h_coeff_notch,1,F,fs);

fresp1 = {'fresp_h1',a};      %ノッチ1を考慮した設計

%前段を考慮したフィルタ
h_coeff_rear=firpm(N_rear,f,fresp1,notch_weight1);
[h_rear,w]=freqz(h_coeff_rear,1,F,fs);

%フィルタ合成
h_coeff_all=conv(h_coeff_notch,h_coeff_rear);
[h_all,w]=freqz(h_coeff_all,1,F,fs);


%%
% 描画
figure
hold on
plot(w/fs2,20*log10(abs(H_HT)));
plot(w/fs2,20*log10(abs(h_notch)));
plot(w/fs2,20*log10(abs(h_rear)));
plot(w/fs2,20*log10(abs(h_all)));
xlim([0 1]);
ylim([-60 20]);
legend('ヒルベルトのみ','ノッチ','後段','ノッチヒルベルト');
xlabel("Amplitude Normalized Frequency (times/pi rad/sample)");
ylabel("Magnitude");
