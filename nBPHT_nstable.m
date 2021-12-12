clear;close all;clc
%**********************************************************************

%---------   ノッチを複数入れたフィルタ   -----------

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
set(0, 'DefaultAxesLineWidth', 1, 'DefaultLineLineWidth', 2); 

%% フィルタ設計部


%各フィルタの次数部
N_HT=20;             %ヒルベルト変換器の次数
N_rear=16;           %後段のフィルタ次数
N_rearHT=N_HT+N_rear;

%各フィルタのfor文
for i=2:24
    Ni_rear = i+2;
    Ni_HT = 30 - Ni_rear;
end
%周波数帯域設定
f_left_stopband  = 0.2;
f_right_stopband = 0.8;
f_left_passband  = 0.3;
f_right_passband = 0.7;



%重み[左側阻止域　通過域　右側阻止域]
Weight = [1 1 1];         


%ヒルベルトフィルタ設計部
fs = 50000;          %Sampling frequency;
fs2 = fs/2;          %ナイキスト周波数
L=15000;             %信号長 15000
ts=1/fs;             %サンプリング周期
t=(0:L-1)*ts;        %サンプリング時間
F=0:0.1:fs2;         %周波数設定
f = [0 f_left_stopband f_left_passband f_right_passband f_right_stopband 1]; %振幅特性の通過域、遷移域、阻止域の周波数指定
a = [0 0 1 1 0 0];         %理想周波数特性[阻止域　通過域　阻止域]
%振幅特性0、振幅特性1、振幅特性0、
%繋がってないところは遷移域


%%
%ヒルベルト変換器

[h_HT,err,res]=firpm(N_HT, f,a,{10000},Weight,'hilbert'); %フィルタ係数の算出
[H_HT,w]=freqz(h_HT,1,F,fs);

%%
%　ノッチフィルタ
nf1 = 3600*180/fs2;    r1=1;        %notch1 & on,off
nf2 = 7200*180/fs2;    r2=1;        %notch2 & on,off
nf3 = 3600*180/fs2;    r3=1;        %notch3 & on,off

D2=1;
D4=1;


notch_freq_1 = 0.1;       %ノッチの場所
notch_freq_2 = 0.15;
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

%% フィルタ合成(ノッチバンドパスフィルタ+ヒルベルト変換)

%ノッチフィルタ
h_coefficient_notch = load('notch3.cof');
[h_amplitude_notch,w]=freqz(h_coefficient_notch,1,F,fs);             %ノッチフィルタの振幅特性

fresp3 = {'fresp_h3',a};      %ノッチ1を考慮した設計

%補正フィルタ
h_coefficient_compensation=firpm(N_rear,f,fresp3,Weight);
[h_amplitude_compensation,w]=freqz(h_coefficient_compensation,1,F,fs);               %補正のフィルタの振幅特性

%フィルタ合成(ノッチフィルタ＋補正バンドパス)
h_coefficient_all=conv(h_coefficient_notch,h_coefficient_compensation);
[h_amplitude_all,w]=freqz(h_coefficient_all,1,F,fs);                 %ノッチバンドパスの振幅特性

%フィルタ合成(ノッチバンドパス＋ヒルベルト変換器)
h_coeff_notchBPHT=conv(h_coefficient_all,h_HT);
[h_amplitude_notchBPHT,w]=freqz(h_coeff_notchBPHT,1,F,fs);
h_amplitude_nBPHT_evaluation(:,1) = w/fs2;
h_amplitude_nBPHT_evaluation(:,2) = h_amplitude_notchBPHT;

%% フィルタ合成(ノッチフィルタ＋ヒルベルト)

%補正ヒルベルトフィルタの作成
h_coefficient_compensationHT=firpm(N_rearHT,f,fresp3,Weight,'hilbert');
[h_amplitude_compensation,w]=freqz(h_coefficient_compensationHT,1,F,fs);               %補正のフィルタの振幅特性

%フィルタ合成(ノッチフィルタ＋ヒルベルト変換器)
h_coefficient_notchHT=conv(h_coefficient_notch,h_coefficient_compensationHT);
[h_amplitude_notchHT,w]=freqz(h_coefficient_notchHT,1,F,fs);
h_amplitude_nHT_evaluation(:,1) = w/fs2;
h_amplitude_nHT_evaluation(:,2) = h_amplitude_notchHT;


%% ノッチバンドパス＋ヒルベルト変換器の最大誤差

%全評価点を3バンドのインデックスに分ける．(左側阻止域　通過域　右側阻止域)
index_nBPHT_pass = find((f_left_passband <= h_amplitude_nBPHT_evaluation(:,1)) & (h_amplitude_nBPHT_evaluation(:,1) <= f_right_passband));
index_nBPHT_leftstop = find((0 <= h_amplitude_nBPHT_evaluation(:,1)) & (h_amplitude_nBPHT_evaluation(:,1) <= f_left_stopband));
index_nBPHT_rightstop = find((f_right_stopband <= h_amplitude_nBPHT_evaluation(:,1)) & (h_amplitude_nBPHT_evaluation(:,1) <= 1));

%3バンドの振幅特性を代入
nBPHT_pass_evaluation(:,2) = h_amplitude_nBPHT_evaluation(index_nBPHT_pass,2);
nBPHT_leftstop_evaluation(:,2) = h_amplitude_nBPHT_evaluation(index_nBPHT_leftstop,2);
nBPHT_rightstop_evaluation(:,2) = h_amplitude_nBPHT_evaluation(index_nBPHT_rightstop,2);

%理想振幅特性と実際の振幅特性の差を求めた後，誤差の最大値を抽出
max_rem_nBPHT_pass=max(abs(nBPHT_pass_evaluation(:,2))-1);
max_rem_nBPHT_leftstop=max(abs(nBPHT_leftstop_evaluation(:,2))-1);
max_rem_nBPHT_rightstop=max(abs(nBPHT_rightstop_evaluation(:,2))-1);
max_rem_nBPHT_all=[abs(max_rem_nBPHT_pass) abs(max_rem_nBPHT_leftstop) abs(max_rem_nBPHT_rightstop)];
max_rem_nBPHT = max(max_rem_nBPHT_all)




%% ノッチバンドパス＋ヒルベルト変換器の最大誤差

%全評価点を3バンドのインデックスに分ける．(左側阻止域　通過域　右側阻止域)
index_nHT_pass = find((f_left_passband <= h_amplitude_nHT_evaluation(:,1)) & (h_amplitude_nHT_evaluation(:,1) <= f_right_passband));
index_nHT_leftstop = find((0 <= h_amplitude_nHT_evaluation(:,1)) & (h_amplitude_nHT_evaluation(:,1) <= f_left_stopband));
index_nHT_rightstop = find((f_right_stopband <= h_amplitude_nHT_evaluation(:,1)) & (h_amplitude_nHT_evaluation(:,1) <= 1));


%3バンドの振幅特性を代入
nHT_pass_evaluation(:,2) = h_amplitude_nHT_evaluation(index_nHT_pass,2);
nHT_leftstop_evaluation(:,2) = h_amplitude_nHT_evaluation(index_nHT_leftstop,2);
nHT_rightstop_evaluation(:,2) = h_amplitude_nHT_evaluation(index_nHT_rightstop,2);


%理想振幅特性と実際の振幅特性の差を求めた後，誤差の最大値を抽出
max_rem_nHT_pass=max(abs(nHT_pass_evaluation(:,2))-1);
max_rem_nHT_leftstop=max(abs(nHT_leftstop_evaluation(:,2))-1);
max_rem_nHT_rightstop=max(abs(nHT_rightstop_evaluation(:,2))-1);
max_rem_nHT_all=[abs(max_rem_nHT_pass) abs(max_rem_nHT_leftstop) abs(max_rem_nHT_rightstop)];
max_rem_nHT = max(max_rem_nHT_all)

%%
% 描画

number_nBPHT = numel(h_coeff_notchBPHT)
number_nHT = numel(h_coefficient_notchHT)


%ヒルベルト変換器
figure
plot(w/fs2,abs(H_HT),"DisplayName","N="+N_HT);
ylim([-0.02 1.1]);
legend()
xlabel("Amplitude Normalized Frequency (\times\pi rad/sample)");
ylabel("Magnitude");

%　ノッチフィルタ
figure
plot(w/fs2,20*log(abs(h_amplitude_notch)));
xlim([0 1]);
ylim([-100 50]);
xlabel("Amplitude Normalized Frequency (\times\pi rad/sample)");
ylabel("Magnitude");


%ノッチバンドパスフィルタ
figure

y_low = -100;
y_high = 50;
x = [notch_freq_1 notch_freq_1];
y = [y_low y_high];
line(x,y,'Color','green','linestyle','--',"DisplayName","notch="+notch_freq_1);
x = [notch_freq_2 notch_freq_2];
y = [y_low y_high];
line(x,y,'Color','green','linestyle','--',"DisplayName","notch="+notch_freq_2);
legend()
xlim([0 1]);
ylim([y_low y_high]);


hold on
plot(w/fs2,20*log10(abs(h_amplitude_notch)),"DisplayName","notch");
plot(w/fs2,20*log10(abs(h_amplitude_compensation)),"DisplayName","bandpass");
plot(w/fs2,20*log10(abs(h_amplitude_all)),"DisplayName","all");
hold off

legend();
% legend('Location','southeast');
xlabel("Normalized Angular Frequency(\times\pi)   (rad/sample)");
ylabel("Magnitude   [dB] ");
exportgraphics(gca,strcat('nbp','.pdf'),'ContentType','vector');



%ノッチバンドパスフィルタ＋ヒルベルト変換　vs ノッチヒルベルト変換
figure1 = figure(4);
close 4;
figure1 = figure(4);

axes1 = axes('Parent',figure1);
hold(axes1,'on');

%ノッチバンドパスフィルタt＋ヒルベルト変換器
plot(w/fs2, 20*log10(abs(h_amplitude_notchBPHT)),"DisplayName","nBPHT");
hold on;
%ノッチヒルベルト変換器
plot(w/fs2, 20*log10(abs(h_amplitude_notchHT)),"DisplayName","nCHT");
x = [notch_freq_1 notch_freq_1];
y = [-200 20];
line(x,y,'Color','green','linestyle','--',"DisplayName","notch="+notch_freq_1);
x = [notch_freq_2 notch_freq_2];
y = [-200 20];
line(x,y,'Color','green','linestyle','--',"DisplayName","notch="+notch_freq_2);
xlim([0 1]);
ylim([-200 20]);
xlabel("Amplitude Normalized Frequency (\times\pi rad/sample)");
ylabel("Magnitude");
box(axes1,'on');
set(axes1,'FontSize',16,'FontName','Times New Roman');
legend('Location','southeast');


axes2 = axes('Parent',figure1,'Position',[0.4 0.68 0.38 0.2]);
hold(axes2,'on');
%ノッチバンドパスフィルタt＋ヒルベルト変換器
plot(w/fs2, 20*log10(abs(h_amplitude_notchBPHT)));
hold on;
%ノッチヒルベルト変換器
plot(w/fs2, 20*log10(abs(h_amplitude_notchHT)));
set(axes2,'FontSize',14,'FontName','Times New Roman');
set(axes2,'Color','w');
box(axes2,'on');
xlim([0.45 0.6]);
ylim([-4 6]);





% %理想特性の図
% figure
% hold on
% grid off
% 
% xticks([0 0.1 0.2 0.25 0.3 0.7 0.8 1]);
% xticklabels({'0','\omega_{b1}','\omega_{b2}','\omega_{s1}','\omega_{p1}','\omega_{p2}','\omega_{s2}','\pi'});
% yticks(1);
% yticklabels({'0'});
% xlim([0 1]);
% ylim([0 1.03]);
% x = [0 0.01  0.09 0.1 0.11 0.19 0.2 0.21 0.25 0.25];
% y = [0 0.1   0.1  0   0.1  0.1  0   0.1  0.1  0];
% line(x,y,'Color','black');
% x = [0.3 0.3 0.7 0.7];
% y = [0 1 1 0];
% line(x,y,'Color','black');
% x = [0.8 0.8 0.9 0.99 1];
% y = [0   0.1 0.1 0.1  0];
% line(x,y,'Color','black');
% 
% xlabel("Normalized Angular Frequency   (rad/sample)");
% ylabel("Magnitude   [dB] ");
% set(gca,'TickLength',[0 0])                 %内向きの目盛りを消す
% exportgraphics(gca,strcat('ideal_resp','.pdf'),'ContentType','vector');

