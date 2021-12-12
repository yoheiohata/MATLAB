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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%---------------------フィルタのパラメータ変更部分------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p=0;
%各フィルタの次数部
for i=0:1:11
    p=p+1;
N_rearHT=26 - i*2;
N_rearHT2=36;
N=N_rearHT+4;

%正規化周波数帯域設定
f_left_stopband  = 0.2;
f_right_stopband = 0.8;
f_left_passband  = 0.3;
f_right_passband = 0.7;


%ノッチの場所
notch_freq_1 = 0.1;       
notch_freq_2 = 0.9;



%重み[左側阻止域　通過域　右側阻止域]
weight = [0.05 0.9 0.05];         

%%--------------------------------------------------------------------

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


%% 信号定義

%振幅
A = 1;
A_noize_rate_1 = 0.7;
A_noize_rate_2 = 0.5;

% インクリメント
q=0;
r=0;


%%

% figureの設定
y_low = -100;
y_high = 20;
limy_high = 140;
limy_low = -150;
%% ノッチヒルベルト用ヒルベルト変換器(手法1と3に用いるやつ)

%ノッチヒルベルト用ヒルベルト変換器
[h_HT,err,res]=firpm(N_rearHT, f,a,{10000},weight,'hilbert'); %フィルタ係数の算出
[H_HT,w]=freqz(h_HT,1,F,fs);

%阻止域を入れた時の振幅特性を比較する用
[h_HT2,err,res]=firpm(N_rearHT2, f,a,{10000},weight,'hilbert'); %フィルタ係数の算出
[H_HT2,w]=freqz(h_HT2,1,F,fs);

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


%% フィルタ合成(ノッチ＋ヒルベルト変換）


h_coefficient_notchHT = conv(h_coefficient_notch,h_HT);
[h_amplitude_notchHT,w] = freqz(h_coefficient_notchHT,1,F,fs);
h_amplitude_nHT_evaluation(:,1) = w/fs2;
h_amplitude_nHT_evaluation(:,2) = h_amplitude_notchHT;


%% フィルタ合成(ノッチフィルタ＋補正ヒルベルト変換)

%補正ヒルベルトフィルタの作成
h_coefficient_compensationHT=firpm(N_rearHT,f,fresp3,weight,'hilbert');
[h_amplitude_compensationCHT,w]=freqz(h_coefficient_compensationHT,1,F,fs);               %補正のフィルタの振幅特性

%フィルタ合成(ノッチフィルタ＋ヒルベルト変換器)
h_coefficient_notchCHT=conv(h_coefficient_notch,h_coefficient_compensationHT);
[h_amplitude_notchCHT,w]=freqz(h_coefficient_notchCHT,1,F,fs);
h_amplitude_nCHT_evaluation(:,1) = w/fs2;
h_amplitude_nCHT_evaluation(:,2) = h_amplitude_notchCHT;

%% nCHTのリプル

%全評価点を3バンドのインデックスに分ける．(左側阻止域　通過域　右側阻止域)
index_nCHT_pass = find((f_left_passband <= h_amplitude_nCHT_evaluation(:,1)) & (h_amplitude_nCHT_evaluation(:,1) <= f_right_passband));
index_nCHT_leftstop = find((0 <= h_amplitude_nCHT_evaluation(:,1)) & (h_amplitude_nCHT_evaluation(:,1) <= f_left_stopband));
index_nCHT_rightstop = find((f_right_stopband <= h_amplitude_nCHT_evaluation(:,1)) & (h_amplitude_nCHT_evaluation(:,1) <= 1));


%3バンドの振幅特性を代入
nCHT_pass_evaluation(:,2) = h_amplitude_nCHT_evaluation(index_nCHT_pass,2);
nCHT_leftstop_evaluation(:,2) = h_amplitude_nCHT_evaluation(index_nCHT_leftstop,2);
nCHT_rightstop_evaluation(:,2) = h_amplitude_nCHT_evaluation(index_nCHT_rightstop,2);


%理想振幅特性と実際の振幅特性の差を求めた後，誤差の最大値を抽出
max_rem_nCHT_pass(p,1)=N;
max_rem_nCHT_pass(p,2)=max(abs(abs(nCHT_pass_evaluation(:,2))-1));
max_rem_nCHT_leftstop=max(abs(nCHT_leftstop_evaluation(:,2))-0);
max_rem_nCHT_rightstop=max(abs(nCHT_rightstop_evaluation(:,2))-0);
max_rem_nCHT_stopall=[20*log(abs(max_rem_nCHT_leftstop)) 20*log(abs(max_rem_nCHT_rightstop))];
disp(max(max_rem_nCHT_stopall));

%% 設計したフィルタに対してSN比を計算

for j=(f_left_passband*fs2):2000:(f_right_passband*fs2)
    q=q+1;
    r=r+1;

    sig_f=j;
    normalized_sig_f = sig_f/fs2;
    sig_f_noize1 = notch_freq_1*fs2;
    sig_f_noize2 = notch_freq_2*fs2;
    
    % 正弦波入力信号
    sig_in = A*cos(2*pi*sig_f*t);
    
    % ノイズの信号
    sig_n1 = A_noize_rate_1*A*cos(2*pi*sig_f_noize1*t);
    sig_n2 = A_noize_rate_2*A*cos(2*pi*sig_f_noize2*t);
    
    %信号
    sig = sig_in + sig_n1;
%     sig = sig_in + sig_n2;
%     sig = sig_in + sig_n1 + sig_n2;

    sig_zeros = [sig_in zeros(1,N_rearHT/2)];
    sig_n_zeros = [sig zeros(1,N_rearHT/2)];

   %% ヒルベルト変換とSN比計算
   
   % 阻止域を有するFIRヒルベルト変換器とSN比
   sig_out = filter(h_coefficient_notchCHT,1,sig_n_zeros);
   SNR = snr(sig_out);
   snr_outcome(q,1) = normalized_sig_f;
   snr_outcome(q,2) = SNR;
   
end
    SNR_outcome(:,1) = max_rem_nCHT_pass(p,1);
    SNR_outcome(:,2) = mean(snr_outcome(:,2),1);
    disp(SNR_outcome)
%% 描画
% nCHT
figure1 = figure(2);
close 2;
figure1 = figure(2);
axes1 = axes('Parent',figure1);
hold(axes1,'on');
x = [notch_freq_1 notch_freq_1];
y = [y_low limy_high];
% ylim([y_low y_high]);
line(x,y,'Color','green','linestyle','--');%,"DisplayName","notch="+notch_freq_1);
x = [notch_freq_2 notch_freq_2];
y = [y_low limy_high];
% ylim([y_low y_high]);
line(x,y,'Color','green','linestyle','--');%,"DisplayName","notch="+notch_freq_2);
x = [f_left_passband f_left_passband];
y = [y_low limy_high];
%ylim([y_low y_high]);
line(x,y,'Color','black','linestyle','--');%,"DisplayName","notch="+notch_freq_1);
x = [f_right_passband f_right_passband];
y = [y_low limy_high];
%ylim([y_low y_high]);
line(x,y,'Color','black','linestyle','--');%,"DisplayName","notch="+notch_freq_2);
% legend('Location','southeast')
hold on
% plot(w/fs2, 20*log(abs(h_amplitude_notch)),"DisplayName","|H_{n}(\omega)|");
% plot(w/fs2, 20*log(abs(h_amplitude_compensationCHT)),"DisplayName","|H_{ht}(\omega)|");
plot(w/fs2, 20*log(abs(h_amplitude_notchCHT)));%,"DisplayName","|H_{nCHT}(\omega)|");
hold off
xlim([0 1]);
% ylim([limy_low limy_high]);
ylim([y_low y_high]);
% legend('Location','southeast')
xlabel("Normalized Angular Frequency   [rad/sample]");
ylabel("Magnitude   [dB] ");

box(axes1,'on');
set(axes1,'FontSize',16,'FontName','Times New Roman');

axes2 = axes('Parent',figure1,'Position',[0.4 0.5 0.25 0.25]);
hold(axes2, 'on');
hold on 
plot(w/fs2, 20*log(abs(h_amplitude_notchCHT)));%,"DisplayName","|H_{nCHT}(\omega)|");
hold off
set(axes2, 'FontSize',14,'FontName','Times New Roman');
set(axes2,'Color','w');
box(axes2,'on');
xlim([0.3 0.7]);
ylim([-2 2]);
% exportgraphics(figure1,strcat('..\卒論\fig\nCHT, N=',num2str(N),'f=',num2str(notch_freq_1),',',num2str(notch_freq_2),', w=',num2str(weight),'.pdf'),'ContentType','vector');




end

disp(SNR_outcome);

filename = 'outdata.xlsx';
table_name = {'N','delta'};
writecell(table_name,filename,'Sheet',1,'Range','D1');
writematrix(max_rem_nCHT_pass, filename,'Sheet',1,'Range','D2');
writematrix(max_rem_nCHT_stopall, filename, 'Sheet',1,'Range','F2');

% 重みを1 1 1にした時のdeltaが0.027
% ノッチの位置で確実に減衰させられれば，阻止域の減衰量もそんないらんくない？
% 阻止域全体で40dBから60dBくらいの減衰量を確保できたら，重みは通過域に回してよりフラットにするべき
% 通過域に重みを一番割くと一番フラットになってそのdeltaが0.005くらいになる
% 通過域deltaが0.027くらいまで許容とすれば，通過域重みマックスとして次数は18次まで落とせる．

