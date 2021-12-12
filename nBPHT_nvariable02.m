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
N_initial=20;


%周波数帯域設定
f_left_stopband  = 0.2;
f_right_stopband = 0.8;
f_left_passband  = 0.3;
f_right_passband = 0.7;

%ノッチの場所
notch_freq_1 = 0.1;       
notch_freq_2 = 0.15;

%重み[左側阻止域　通過域　右側阻止域]
Weight = [1 1 1];  

%誤差の差
epsilon = 0.005;


%ヒルベルトフィルタ設計部
fs = 50000;          %Sampling frequency;
fs2 = fs/2;          %ナイキスト周波数
L=15000;             %信号長 15000
ts=1/fs;             %サンプリング周期
t=(0:L-1)*ts;        %サンプリング時間
F=0:0.1:fs2;         %周波数設定
f = [0 f_left_stopband f_left_passband f_right_passband f_right_stopband 1]; %振幅特性の通過域、遷移域、阻止域の周波数指定
a = [0 0 1 1 0 0];         %理想周波数特性[阻止域0　通過域1　阻止域0]
N_all = 0;



%% ノッチフィルタ設計部

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


%ノッチフィルタ
h_coefficient_notch = load('notch3.cof');
%ノッチフィルタの振幅特性
[h_amplitude_notch,w]=freqz(h_coefficient_notch,1,F,fs);             
%ノッチ1を考慮した設計
fresp3 = {'fresp_h3',a};      
%%

N_all = N_initial;

%% ノッチヒルベルト用ヒルベルト変換器

%フィルタ係数の算出
[h_HT,err,res]=firpm(N_all, f,a,{10000},Weight,'hilbert'); 
[H_HT,w]=freqz(h_HT,1,F,fs);



%% フィルタ合成(ノッチフィルタ＋ヒルベルト)


%補正ヒルベルトフィルタの作成
h_coefficient_compensationHT=firpm(N_all,f,fresp3,Weight,'hilbert');
[h_amplitude_compensation,w]=freqz(h_coefficient_compensationHT,1,F,fs);

%フィルタ合成(ノッチフィルタ＋ヒルベルト変換器)
h_coefficient_notchHT=conv(h_coefficient_notch,h_coefficient_compensationHT);
[h_amplitude_notchHT,w]=freqz(h_coefficient_notchHT,1,F,fs);

%評価点と振幅特性のベクトル
h_amplitude_nHT_evaluation(:,1) = w/fs2;
h_amplitude_nHT_evaluation(:,2) = h_amplitude_notchHT;



%% ノッチヒルベルト変換器の最大誤差

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
%%
while true    
j=0;
for i = 2:2:N_all-8      
                j=j+1;
                Ni_rear = i+2;
                Ni_HT = N_all - Ni_rear;
    
    

%ノッチバンドパス用ヒルベルト変換器
[h_BPHT,err,res]=firpm(Ni_HT, f,a,{10000},Weight,'hilbert');
[H_BPHT,w]=freqz(h_BPHT,1,F,fs);


%% フィルタ合成(ノッチバンドパスフィルタ+ヒルベルト変換)



%補正フィルタ
h_coefficient_compensation=firpm(Ni_rear,f,fresp3,Weight);
[h_amplitude_compensation,w]=freqz(h_coefficient_compensation,1,F,fs);


%フィルタ合成(ノッチフィルタ＋補正バンドパス)
h_coefficient_all=conv(h_coefficient_notch,h_coefficient_compensation);
[h_amplitude_all,w]=freqz(h_coefficient_all,1,F,fs);


%フィルタ合成(ノッチバンドパス＋ヒルベルト変換器)
h_coeff_notchBPHT=conv(h_coefficient_all,h_BPHT);
[h_amplitude_notchBPHT,w]=freqz(h_coeff_notchBPHT,1,F,fs);


%評価点と振幅特性のベクトル
h_amplitude_nBPHT_evaluation(:,1) = w/fs2;
h_amplitude_nBPHT_evaluation(:,2) = h_amplitude_notchBPHT;


%% ノッチバンドパス＋ヒルベルト変換器の最大誤差

%全評価点を3バンドのインデックスに分ける．(左側阻止域　通過域　右側阻止域)
index_nBPHT_pass = find((f_left_passband <= h_amplitude_nBPHT_evaluation(:,1)) & (h_amplitude_nBPHT_evaluation(:,1) <= f_right_passband));
index_nBPHT_leftstop = find((0 <= h_amplitude_nBPHT_evaluation(:,1)) & (h_amplitude_nBPHT_evaluation(:,1) <= f_left_stopband));
index_nBPHT_rightstop = find((f_right_stopband <= h_amplitude_nBPHT_evaluation(:,1)) & (h_amplitude_nBPHT_evaluation(:,1) <= 1));


%3バンドの振幅特性を代入
%nBPHT_pass_evaluation(:,2) = index_nBPHT_pass;
%nBPHT_pass_evaluation(:,1) = h_amplitude_nBPHT_evaluation(index_nBPHT_pass,1);
nBPHT_pass_evaluation(:,2) = h_amplitude_nBPHT_evaluation(index_nBPHT_pass,2);
nBPHT_leftstop_evaluation(:,2) = h_amplitude_nBPHT_evaluation(index_nBPHT_leftstop,2);
nBPHT_rightstop_evaluation(:,2) = h_amplitude_nBPHT_evaluation(index_nBPHT_rightstop,2);


%理想振幅特性と実際の振幅特性の差を求めた後，誤差の最大値を抽出
max_rem_nBPHT_pass=max(abs(abs(nBPHT_pass_evaluation(:,2))-1));
max_rem_nBPHT_leftstop=max(abs(nBPHT_leftstop_evaluation(:,2))-0);
max_rem_nBPHT_rightstop=max(abs(nBPHT_rightstop_evaluation(:,2))-0);
max_rem_nBPHT_all=[abs(max_rem_nBPHT_pass) abs(max_rem_nBPHT_leftstop) abs(max_rem_nBPHT_rightstop)];

max_rem_nBPHT(:,1) = Ni_rear;
max_rem_nBPHT(:,2) = max(max_rem_nBPHT_all);
% outcome1 = [Ni_rear Ni_HT max_rem_nBPHT];
% display('    N_rear　   N_HT      nBPHT　   ');
% disp(outcome1);



ref_nBPHT_max(j,2) = max_rem_nBPHT(:,2);
ref_nBPHT_coefficient(:,j) = h_amplitude_nBPHT_evaluation(:,2);
ref_nBPHT_N_rear(:,j) = max_rem_nBPHT(:,1);
[min_nBPHT, minInd] = min(ref_nBPHT_max(:,2));
N_rear_ref = ref_nBPHT_N_rear(:,minInd);
h_cor = ref_nBPHT_coefficient(:,minInd);


end

max_rem_nBPHT = min(ref_nBPHT_max(:,2));
max_rem_nHT = max(max_rem_nHT_all);
dB = 20*log10(abs(max_rem_nHT-max_rem_nBPHT));
outcome = [N_all N_rear_ref max_rem_nBPHT max_rem_nHT dB];
display('   N　        N_rear    nBPHT　    nHT     dB');
disp(outcome)


    if (max_rem_nBPHT) - (max_rem_nHT) > epsilon
        N_all = N_all + 2;  
    else (max_rem_nBPHT - max_rem_nHT) < epsilon
        break
    end
    

    
%% 描画


% figure
% %ノッチバンドパスフィルタ＋ヒルベルト変換　vs ノッチ補正ヒルベルト変換
% 
% 
% x = [notch_freq_1 notch_freq_1];
% y = [-200 20];
% line(x,y,'Color','green','linestyle','--',"DisplayName","notch="+notch_freq_1);
% x = [notch_freq_2 notch_freq_2];
% y = [-200 20];
% line(x,y,'Color','green','linestyle','--',"DisplayName","notch="+notch_freq_2);
% legend()
% xlim([0 1]);
% ylim([-200 20]);
% 
% hold on
% %ノッチバンドパスフィルタ＋ヒルベルト変換器
% plot(w/fs2, 20*log10(abs(h_cor)), "DisplayName","nBPHT Nrear="+N_rear_ref);
% %ノッチ補正ヒルベルト変換器
% plot(w/fs2, 20*log10(abs(h_amplitude_notchHT)), "DisplayName", "nCHT");
% hold off
% 
% legend()
% xlabel("Amplitude Normalized Frequency (\times\pi rad/sample)");
% ylabel("Magnitude");

         
    

end


 





