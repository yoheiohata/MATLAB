clear;close all;clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------   ノッチを複数入れたフィルタ   -----------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%

%　ノッチフィルタが複数個の場合

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
%振幅特性0、振幅特性1、振幅特性0、
%繋がってないところは遷移域
Weight = [3 30 3];         %重み[左側阻止域　通過域　右側阻止域]

%%
%ヒルベルト変換器

[h_HT,err,res]=firpm(N_HT, f,a,{10000},Weight,'Hilbert'); %フィルタ係数の算出
[H_HT,W]=freqz(h_HT,1,F,fs);

%%
%　ノッチフィルタ
N1=12; %フィルタ次数
f1 = [0 0.1 0.2 0.8 0.9 1 ];
% f1=[   0.0        0.1       ... % 第 1 帯域の下限周波数 上限周波数（通過域）
%        3000/fs2     1             ... % 第 2 帯域の下限周波数 上限周波数（阻止域）
%    ]; 

%  2nd フィルタ ------------------
N2=28; %フィルタ次数
% f2=[   0.0             0.0001        ... % 第 1 帯域の下限周波数 上限周波数（通過域）
%        3000/fs2     1              ... % 第 2 帯域の下限周波数 上限周波数（阻止域）
%        ];
f2 = [0 0.1 0.2 0.8 0.9 1 ];

% 3rd フィルタ ----------------------
N3=8; %フィルタ次数
% f3=[   0.0            0.0001         ... % 第 1 帯域の下限周波数 上限周波数（通過域）
%        3000/fs2    1               ... % 第 2 帯域の下限周波数 上限周波数（阻止域）
%        ];
f3 = [0 0.1 0.2 0.8 0.9 1 ];

nf1 = 3600*180/fs2;    r1=1;        %notch1 & on,off
nf2 = 7200*180/fs2;    r2=1;        %notch2 & on,off
nf3 = 3600*180/fs2;    r3=1;        %notch3 & on,off

D2=1;
D4=1;

N_rear=30;                 %フィルタ次数
notch_freq_1 = 0.1;       %ノッチの場所
notch_freq_2 = 0.2;
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

%%
%フィルタ合成
%h_coeff_notch = load('notch1.cof');
%fresp1 = {'fresp_h1',a};      %ノッチ1を考慮した設計
%h_coeff_rear=firpm(N_rear,f,fresp1,notch_weight1);
h_coeff_rear1=firpm(N1,f1,notch_amplitude1,notch_weight1);
h_coeff_rear1=h_coeff_rear1./sum(h_coeff_rear1);
H_coeff_rear1=freqz(h_coeff_rear1,1,F,fs);
save h_coeff_rear1.cof h_coeff_rear1 -ascii;

h_coeff_rear2 = firpm(N2,f2,notch_amplitude1,notch_weight1);
h_coeff_rear2 = h_coeff_rear2./sum(h_coeff_rear2);

hh2=conv(h_coeff_rear2,notch3);

hh2=hh2./sum(hh2);
H_notch3=freqz(notch3,1,F,fs/D2);
H_coeff_rear2=freqz(h_coeff_rear2,1,F,fs/D2);
HH2=freqz(hh2,1,F,fs/D2);
save hh2.cof hh2 -ascii;


h_coeff_rear3=firpm(N3,f3,notch_amplitude1,notch_weight1);
h_coeff_rear3=h_coeff_rear3./sum(h_coeff_rear3);

hh3=conv(h_coeff_rear3,notch3);
hh3=hh3./sum(hh3);

HH3 = freqz(hh3,1,F,fs/D4);
H3_after=freqz(h_coeff_rear3,1,F,fs/D4);
H3_before=freqz(notch3,1,F,fs/D4);
save hh3.cof hh3 -ascii;

h_coeff_notch = load('hh3.cof');
fresp3 = {'fresp_h3',a};
h_coeff_rear = H3_after;

h_coeff_all=conv(h_coeff_notch,h_coeff_rear);      %
[h_notch,w]=freqz(h_coeff_notch,1,F,fs);           %freqz:frequency response
[h_rear,w]=freqz(h_coeff_rear,1,F,fs);
[h_all,w]=freqz(h_coeff_all,1,F,fs);

%%
% 描画
figure
hold on
plot(w/fs2,20*log10(abs(h_notch)));
plot(w/fs2,20*log10(abs(h_rear)));
plot(w/fs2,20*log10(abs(h_all)));
xlim([0 1]);
ylim([-60 20]);
legend('Notch','Rear','All');
xlabel("Amplitude Normalized Frequency (times/pi rad/sample)");
ylabel("Magnitude");
