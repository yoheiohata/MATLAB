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
display('    N_rear　 N_HT　    nHT         nBPHT　   nCHT    dB(nHT-nBPHT)  dB(nBPHT-nCHT)');

%% フィルタ設計部

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%---------------------フィルタのパラメータ変更部分------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%各フィルタの次数部
N_rearHT=30;
N_rearHT2=36;

%正規化周波数帯域設定
f_left_stopband  = 0.2;
f_right_stopband = 0.8;
f_left_passband  = 0.3;
f_right_passband = 0.7;


%ノッチの場所
notch_freq_1 = 0.1;       
notch_freq_2 = 0.9;



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
f = [0 f_left_stopband f_left_passband f_right_passband f_right_stopband 1]; %振幅特性の通過域、遷移域、阻止域の周波数指定
a = [0 0 1 1 0 0];         %理想周波数特性[阻止域　通過域　阻止域]
%振幅特性0、振幅特性1、振幅特性0、
%繋がってないところは遷移域

% figureの設定
y_low = -100;
y_high = 20;
limy_high = 140;
limy_low = -150;
%% ノッチヒルベルト用ヒルベルト変換器(手法1と3に用いるやつ)

%ノッチヒルベルト用ヒルベルト変換器
[h_HT,err,res]=firpm(N_rearHT, f,a,{10000},Weight,'hilbert'); %フィルタ係数の算出
[H_HT,w]=freqz(h_HT,1,F,fs);

%阻止域を入れた時の振幅特性を比較する用
[h_HT2,err,res]=firpm(N_rearHT2, f,a,{10000},Weight,'hilbert'); %フィルタ係数の算出
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


%% ノッチ＋ヒルベルト変換器の最大誤差

%全評価点を3バンドのインデックスに分ける．(左側阻止域　通過域　右側阻止域)
index_nHT_pass = find((f_left_passband <= h_amplitude_nHT_evaluation(:,1)) & (h_amplitude_nHT_evaluation(:,1) <= f_right_passband));
index_nHT_leftstop = find((0 <= h_amplitude_nHT_evaluation(:,1)) & (h_amplitude_nHT_evaluation(:,1) <= f_left_stopband));
index_nHT_rightstop = find((f_right_stopband <= h_amplitude_nHT_evaluation(:,1)) & (h_amplitude_nHT_evaluation(:,1) <= 1));


%3バンドの振幅特性を代入
nHT_pass_evaluation(:,2) = h_amplitude_nHT_evaluation(index_nHT_pass,2);
nHT_leftstop_evaluation(:,2) = h_amplitude_nHT_evaluation(index_nHT_leftstop,2);
nHT_rightstop_evaluation(:,2) = h_amplitude_nHT_evaluation(index_nHT_rightstop,2);


%理想振幅特性と実際の振幅特性の差を求めた後，誤差の最大値を抽出
max_rem_nHT_pass=max(abs(abs(nHT_pass_evaluation(:,2))-1));
max_rem_nHT_leftstop=max(abs(nHT_leftstop_evaluation(:,2))-0);
max_rem_nHT_rightstop=max(abs(nHT_rightstop_evaluation(:,2))-0);
max_rem_nHT_all=[abs(max_rem_nHT_pass) abs(max_rem_nHT_leftstop) abs(max_rem_nHT_rightstop)];


%% フィルタ合成(ノッチフィルタ＋補正ヒルベルト変換)

%補正ヒルベルトフィルタの作成
h_coefficient_compensationHT=firpm(N_rearHT,f,fresp3,Weight,'hilbert');
[h_amplitude_compensationCHT,w]=freqz(h_coefficient_compensationHT,1,F,fs);               %補正のフィルタの振幅特性

%フィルタ合成(ノッチフィルタ＋ヒルベルト変換器)
h_coefficient_notchCHT=conv(h_coefficient_notch,h_coefficient_compensationHT);
[h_amplitude_notchCHT,w]=freqz(h_coefficient_notchCHT,1,F,fs);
h_amplitude_nCHT_evaluation(:,1) = w/fs2;
h_amplitude_nCHT_evaluation(:,2) = h_amplitude_notchCHT;


%% ノッチ＋補正ヒルベルト変換器の最大誤差

%全評価点を3バンドのインデックスに分ける．(左側阻止域　通過域　右側阻止域)
index_nCHT_pass = find((f_left_passband <= h_amplitude_nCHT_evaluation(:,1)) & (h_amplitude_nCHT_evaluation(:,1) <= f_right_passband));
index_nCHT_leftstop = find((0 <= h_amplitude_nCHT_evaluation(:,1)) & (h_amplitude_nCHT_evaluation(:,1) <= f_left_stopband));
index_nCHT_rightstop = find((f_right_stopband <= h_amplitude_nCHT_evaluation(:,1)) & (h_amplitude_nCHT_evaluation(:,1) <= 1));


%3バンドの振幅特性を代入
nCHT_pass_evaluation(:,2) = h_amplitude_nCHT_evaluation(index_nCHT_pass,2);
nCHT_leftstop_evaluation(:,2) = h_amplitude_nCHT_evaluation(index_nCHT_leftstop,2);
nCHT_rightstop_evaluation(:,2) = h_amplitude_nCHT_evaluation(index_nCHT_rightstop,2);


%理想振幅特性と実際の振幅特性の差を求めた後，誤差の最大値を抽出
max_rem_nCHT_pass=max(abs(abs(nCHT_pass_evaluation(:,2))-1));
max_rem_nCHT_leftstop=max(abs(nCHT_leftstop_evaluation(:,2))-0);
max_rem_nCHT_rightstop=max(abs(nCHT_rightstop_evaluation(:,2))-0);
max_rem_nCHT_all=[abs(max_rem_nCHT_pass) abs(max_rem_nCHT_leftstop) abs(max_rem_nCHT_rightstop)];



%%

for i=2:2:N_rearHT-6
    Ni_rear = i+2;
    Ni_HT = N_rearHT - Ni_rear;
    

%ノッチバンドパス用ヒルベルト変換器
[h_BPHT,err,res]=firpm(Ni_HT, f,a,{10000},Weight,'hilbert'); %フィルタ係数の算出
[H_BPHT,w]=freqz(h_BPHT,1,F,fs);


%% フィルタ合成(ノッチバンドパスフィルタ+ヒルベルト変換)


%補正フィルタ設計
h_coefficient_compensationBP=firpm(Ni_rear,f,fresp3,Weight);
[h_amplitude_compensationBPHT,w]=freqz(h_coefficient_compensationBP,1,F,fs);               %補正のフィルタの振幅特性

%フィルタ合成(ノッチフィルタ＋補正バンドパス)
h_coefficient_all=conv(h_coefficient_notch,h_coefficient_compensationBP);
[h_amplitude_all,w]=freqz(h_coefficient_all,1,F,fs);                 %ノッチバンドパスの振幅特性

%フィルタ合成(ノッチバンドパス＋ヒルベルト変換器)
h_coefficient_notchBPHT=conv(h_coefficient_all,h_BPHT);
[h_amplitude_notchBPHT,w]=freqz(h_coefficient_notchBPHT,1,F,fs);

%評価点と振幅特性のベクトル
h_amplitude_nBPHT_evaluation(:,1) = w/fs2;
h_amplitude_nBPHT_evaluation(:,2) = h_amplitude_notchBPHT;



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
max_rem_nBPHT_pass=max(abs(abs(nBPHT_pass_evaluation(:,2))-1));
max_rem_nBPHT_leftstop=max(abs(nBPHT_leftstop_evaluation(:,2))-0);
max_rem_nBPHT_rightstop=max(abs(nBPHT_rightstop_evaluation(:,2))-0);
max_rem_nBPHT_all=[abs(max_rem_nBPHT_pass) abs(max_rem_nBPHT_leftstop) abs(max_rem_nBPHT_rightstop)];


%%
% 描画

number_nCHT = numel(h_coefficient_notchCHT);
max_rem_nHT = max(max_rem_nHT_all);
max_rem_nBPHT = max(max_rem_nBPHT_all);
max_rem_nCHT = max(max_rem_nCHT_all);
dB_nHT_nBPHT = -20*log10(abs(max_rem_nHT-max_rem_nBPHT));
dB_nBPHT_nCHT = -20*log10(abs(max_rem_nCHT-max_rem_nBPHT));
outcome = [Ni_rear Ni_HT max_rem_nHT max_rem_nBPHT max_rem_nCHT dB_nHT_nBPHT dB_nBPHT_nCHT];
disp(outcome)


%notchBP
figure(1);
x = [notch_freq_1 notch_freq_1];
y = [-100 20];
line(x,y,'Color','green','linestyle','--');%,"DisplayName","notch="+notch_freq_1);
x = [notch_freq_2 notch_freq_2];
y = [-100 20];
line(x,y,'Color','green','linestyle','--');%,"DisplayName","notch="+notch_freq_2);
xlim([0 1]);
ylim([-100 20]);
hold on
plot(w/fs2,20*log10(abs(h_amplitude_notch)));%,"DisplayName","notch");
plot(w/fs2,20*log10(abs(h_amplitude_compensationBPHT)));%,"DisplayName","rear");
plot(w/fs2,20*log10(abs(h_amplitude_all)));%,"DisplayName","notch-bandpass");
hold off
xlabel("Normalized Angular Frequency   [rad/sample]");
ylabel("Magnitude       [dB] ");
% exportgraphics(gca,strcat('..\CAS\fig\bp','.pdf'),'ContentType','vector');
% exportgraphics(gca,strcat('..\CAS_final\fig\bp','.pdf'),'ContentType','vector');

% nBPHT
% figure100 = figure(10);
figure
x = [notch_freq_1 notch_freq_1];
% y = [y_low y_high];
y = [limy_low limy_high];
line(x,y,'Color','green','linestyle','--',"DisplayName","notch="+notch_freq_1);
x = [notch_freq_2 notch_freq_2];
% y = [y_low y_high];
y = [limy_low limy_high];
line(x,y,'Color','green','linestyle','--',"DisplayName","notch="+notch_freq_2);
legend()
xlim([0 1]);
% ylim([y_low y_high]);
ylim([limy_low limy_high]);
hold on
plot(w/fs2,20*log10(abs(h_amplitude_notch)),"DisplayName","|H_n(\omega)|");
plot(w/fs2,20*log10(abs(h_amplitude_compensationBPHT)),"DisplayName","|\^H_{bp}(\omega)|");
plot(w/fs2,20*log10(abs(H_HT)),"DisplayName","|H_{ht}(\omega)|");
plot(w/fs2,20*log10(abs(h_amplitude_notchBPHT)),"DisplayName","|H_{nBPHT}(\omega)|");
hold off
legend('Location','southeast')
% legend("|H_n(\omega)|","|\^H_{bp}(\omega)|","|H_{ht}(\omega)|","|H_{nBPHT}(\omega)|",'Location','southeast','$$\hat{x}$$','Interpreter','Latex');
xlabel("Normalized Angular Frequency(\times\pi)   [rad/sample]");
ylabel("Magnitude   [dB] ");
% exportgraphics(gca,strcat('..\CAS\fig\nBPHT\rear',num2str(Ni_rear),'ht',num2str(Ni_HT),'.pdf'),'ContentType','vector');
% exportgraphics(gca,strcat('..\CAS_final\fig\nBPHT\rear',num2str(Ni_rear),'ht',num2str(Ni_HT),'.pdf'),'ContentType','vector');


% % nHT vs nBPHT
% figure
% x = [notch_freq_1 notch_freq_1];
% y = [-200 140];
% line(x,y,'Color','green','linestyle','--',"DisplayName","notch="+notch_freq_1);
% x = [notch_freq_2 notch_freq_2];
% y = [-200 140];
% line(x,y,'Color','green','linestyle','--',"DisplayName","notch="+notch_freq_2);
% legend('Location','southeast')
% xlim([0 1]);
% ylim([-200 140]);
% hold on
% plot(w/fs2, 20*log(abs(h_amplitude_notchHT)),"DisplayName","nHT");
% plot(w/fs2, 20*log10(abs(h_amplitude_notchBPHT)),"DisplayName","nBPHT Nrear="+Ni_rear);
% hold off
% legend('Location','southeast')
% xlabel("Normalized Angular Frequency   [rad/sample]");
% ylabel("Magnitude   [dB] ");
% exportgraphics(gca,strcat('..\CAS\fig\nHT_nBPHT\nHT_nBPHT_rear',num2str(Ni_rear),'ht',num2str(Ni_HT),'.pdf'),'ContentType','vector');
% exportgraphics(gca,strcat('..\CAS_final\fig\nHT_nBPHT\nHT_nBPHT_rear',num2str(Ni_rear),'ht',num2str(Ni_HT),'.pdf'),'ContentType','vector');


%nBPHT vs nCHT
figure1 = figure(3);
close 3;
figure1 = figure(3);
axes1 = axes('Parent',figure1);
hold(axes1,'on');
x = [notch_freq_1 notch_freq_1];
y = [y_low y_high];
line(x,y,'Color','green','linestyle','--',"DisplayName","notch="+notch_freq_1);
x = [notch_freq_2 notch_freq_2];
y = [y_low y_high];
line(x,y,'Color','green','linestyle','--',"DisplayName","notch="+notch_freq_2);
legend('Location','southeast')
xlim([0 1]);
ylim([y_low y_high]);
hold on
plot(w/fs2, 20*log10(abs(h_amplitude_notchBPHT)),"DisplayName","|H_{nBPHT}(\omega)|");
plot(w/fs2, 20*log10(abs(h_amplitude_notchCHT)),"DisplayName", "|H_{nCHT}(\omega)|");
hold off
legend('Location','southeast')
xlabel("Normalized Angular Frequency   [rad/sample]");
ylabel("Magnitude   [dB] ");
box(axes1,'on');
set(axes1,'FontSize',16,'FontName','Times New Roman');

axes2 = axes('Parent',figure1,'Position',[0.375 0.4 0.3 0.3]);
hold(axes2,'on');
hold on
plot(w/fs2, 20*log10(abs(h_amplitude_notchBPHT)),"DisplayName","|H_{nBPHT}(\omega)| N_{bp}="+Ni_rear);
plot(w/fs2, 20*log10(abs(h_amplitude_notchCHT)),"DisplayName", "|H_{nCHT}(\omega)|");
hold off
set(axes2,'FontSize',14,'FontName','Times New Roman');
set(axes2,'Color','w');
box(axes2,'on');
xlim([0.3 0.7]);
ylim([-2 2]);
%exportgraphics(figure1,strcat('..\CAS\fig\nCHT\n_sta_rear',num2str(Ni_rear),'ht',num2str(Ni_HT),'.pdf'),'ContentType','vector');
%exportgraphics(figure1,strcat('..\CAS_final\fig\nCHT\n_sta_rear',num2str(Ni_rear),'ht',num2str(Ni_HT),'.pdf'),'ContentType','vector');
end

% nHT
figure;
x = [notch_freq_1 notch_freq_1];
y = [limy_low limy_high];
line(x,y,'Color','green','linestyle','--');%,"DisplayName","notch="+notch_freq_1);
x = [notch_freq_2 notch_freq_2];
y = [limy_low limy_high];
line(x,y,'Color','green','linestyle','--');%,"DisplayName","notch="+notch_freq_2);
% legend('Location','southeast')
hold on
plot(w/fs2, 20*log(abs(h_amplitude_notch)));%,"DisplayName","|H_{n}(\omega)|");
plot(w/fs2, 20*log(abs(H_HT)));%,"DisplayName","|H_{ht}(\omega)|");
plot(w/fs2, 20*log(abs(h_amplitude_notchHT)));%,"DisplayName","|H_{nHT}(\omega)|");
hold off
xlim([0 1]);
ylim([limy_low limy_high]);
% ylim([y_low y_high]);
% legend('Location','southeast')
xlabel("Normalized Angular Frequency   [rad/sample]");
ylabel("Magnitude   [dB] ");
% exportgraphics(gca,strcat('..\CAS\fig\nHT_1','.pdf'),'ContentType','vector');
% exportgraphics(gca,strcat('..\CAS_final\fig\nHT_1','.pdf'),'ContentType','vector');

% nCHT
figure;
x = [notch_freq_1 notch_freq_1];
y = [limy_low limy_high];
% ylim([y_low y_high]);
line(x,y,'Color','green','linestyle','--');%,"DisplayName","notch="+notch_freq_1);
x = [notch_freq_2 notch_freq_2];
y = [limy_low limy_high];
% ylim([y_low y_high]);
line(x,y,'Color','green','linestyle','--');%,"DisplayName","notch="+notch_freq_2);
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
% exportgraphics(gca,strcat('..\CAS\fig\nCHT_1','.pdf'),'ContentType','vector');
% exportgraphics(gca,strcat('..\CAS_final\fig\nCHT_1','.pdf'),'ContentType','vector');




% % nHT vs nCHT
% figure;
% limy = 140;
% x = [notch_freq_1 notch_freq_1];
% y = [y_low limy];
% line(x,y,'Color','green','linestyle','--',"DisplayName","notch="+notch_freq_1);
% x = [notch_freq_2 notch_freq_2];
% y = [y_low limy];
% line(x,y,'Color','green','linestyle','--',"DisplayName","notch="+notch_freq_2);
% legend('Location','southeast')
% hold on
% plot(w/fs2, 20*log(abs(h_amplitude_notchHT)),"DisplayName","nHT");
% plot(w/fs2, 20*log10(abs(h_amplitude_notchCHT)),"DisplayName","nCHT");
% hold off
% xlim([0 1]);
% ylim([y_low limy]);
% legend('Location','southeast')
% xlabel("Normalized Angular Frequency(\times\pi)   (rad/sample)");
% ylabel("Magnitude   [dB] ");
% exportgraphics(gca,strcat('..\CAS\fig\nCHT_nHT','.pdf'),'ContentType','vector');
% exportgraphics(gca,strcat('..\CAS_final\fig\nCHT_nHT','.pdf'),'ContentType','vector');

% notch filter だけ
figure;
x = [notch_freq_1 notch_freq_1];
y = [limy_low limy_high];
line(x,y,'Color','green','linestyle','--');%,"DisplayName","notch="+notch_freq_1);
x = [notch_freq_2 notch_freq_2];
y = [limy_low limy_high];
line(x,y,'Color','green','linestyle','--');%,"DisplayName","notch="+notch_freq_2);
% legend('Location','southeast')
hold on
plot(w/fs2, 20*log(abs(h_amplitude_notch)));%,"DisplayName","|H_{n}(\omega)|");
hold off
xlim([0 1]);
ylim([limy_low limy_high]);
% ylim([y_low y_high]);
% legend('Location','southeast')
xlabel("Normalized Angular Frequency   [rad/sample]");
ylabel("Magnitude   [dB] ");
% exportgraphics(gca,strcat('..\CAS\fig\notch','.pdf'),'ContentType','vector');
% exportgraphics(gca,strcat('..\CAS_final\fig\notch','.pdf'),'ContentType','vector');


% HTだけ
figure
hold on
% plot(w/fs2,20*log(abs(H_HT2)),"DisplayName","N="+N_rearHT2);
plot(w/fs2,20*log(abs(H_HT)));%,"DisplayName","N="+N_rearHT);
hold off
xlabel("Normalized Angular Frequency(\times\pi)   (rad/sample)");
ylabel("Magnitude    [dB] ");
%legend('Location','southeast')
xlim([0 1]);
ylim([-100 20]);
% exportgraphics(gca,strcat('..\CAS\fig\stopw_ht_resp','.pdf'),'ContentType','vector');
% exportgraphics(gca,strcat('..\CAS_final\fig\stopw_ht_resp','.pdf'),'ContentType','vector');



%理想特性の図
figure
hold on
grid off
xticks([0 0.1 0.2 0.25 0.3 0.7 0.8 1]);
xticklabels({'0','\omega_{b1}','\omega_{b2}','\omega_{s1}','\omega_{p1}','\omega_{p2}','\omega_{s2}','\pi'});
yticks(1);
yticklabels({'0'});
xlim([0 1]);
ylim([0 1.03]);
x = [0 0.01  0.09 0.1 0.11 0.19 0.2 0.21 0.25 0.25];
y = [0 0.1   0.1  0   0.1  0.1  0   0.1  0.1  0];
line(x,y,'Color','black');
x = [0.3 0.3 0.7 0.7];
y = [0 1 1 0];
line(x,y,'Color','black');
x = [0.8 0.8 0.9 0.99 1];
y = [0   0.1 0.1 0.1  0];
line(x,y,'Color','black');
xlabel("Normalized Angular Frequency   [rad/sample]");
ylabel("Magnitude   [dB] ");
set(gca,'TickLength',[0 0])                 %内向きの目盛りを消す
% exportgraphics(gca,strcat('..\CAS\fig\ideal_resp','.pdf'),'ContentType','vector');
% exportgraphics(gca,strcat('..\CAS_final\fig\ideal_resp','.pdf'),'ContentType','vector');

