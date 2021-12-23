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
%% params

% filter params
N=30;
f_left_stop  = 0.2;
f_left_pass  = 0.3;
f_right_pass = 0.7;
f_right_stop = 0.8;
notch_f1 = 0.1;       
notch_f2 = 0.9;
Weight = [1 1 1];         
fs = 50000;
fs2 = fs/2;
L=15000;
ts=1/fs;
t=(0:L-1)*ts;
F=0:0.1:fs2;
f = [0 f_left_stop f_left_pass f_right_pass f_right_stop 1];
a = [0 0 1 1 0 0];

% signal params(signal frequency is below)
A=1;
A_noize_rate_1=0.7;
A_noize_rate_2=0.5;

% figure params
p=0;
y_low = -150;
y_high = 20;

%% design filter

% conventional
% bandpass hilbert
[h_BPHT,err,res]=firpm(N, f, a, {10000}, Weight,'hilbert'); %フィルタ係数の算出
[H_BPHT,w]=freqz(h_BPHT,1,F,fs);


% proposed
% notch
nf1 = 3600*180/fs2;    r1=1;    D2=1;
nf2 = 7200*180/fs2;    r2=1;    D4=1;
nf3 = 3600*180/fs2;    r3=1;
notch_weight1=[1 1 1];
notch_amplitude1=[0 0 1 1 0 0];
[notch1,notch2,notch3]=notch(nf1,nf2,nf3,D2,D4,fs2,r1,r2,r3);
notch1 = [1 -2*cos(notch_f1*pi) 1*1];
notch2 = [1 -2*cos(notch_f2*pi) 1*1];
notch3 = conv(notch1,notch2);
notch1 = notch1/sum(notch1);
notch2 = notch2/sum(notch2);
notch3 = notch3/sum(notch3);
save notch1.cof notch1 -ascii;
save notch2.cof notch2 -ascii;
save notch3.cof notch3 -ascii;
fresp3 = {'fresp_h3',a};
h_notch = load('notch3.cof');
[H_notch,w]=freqz(h_notch,1,F,fs);   

% hilbert
h_HT=firpm(N-2,f,fresp3,Weight,'hilbert');
[H_HT,w]=freqz(h_HT,1,F,fs);

% convolution
h_CHT=conv(h_notch,h_HT);
[H_CHT,w]=freqz(h_CHT,1,F,fs);

%% draw filter
figure;
x = [notch_f1 notch_f1];
y = [y_low y_high];
line(x,y,'Color','green','linestyle','--',"DisplayName","noize="+notch_f1);
x = [notch_f2 notch_f2];
y = [y_low y_high];
line(x,y,'Color','green','linestyle','--',"DisplayName","noize="+notch_f2);
hold on
plot(w/fs2, 20*log(abs(H_BPHT)),"DisplayName","BPHT")
plot(w/fs2, 20*log(abs(H_CHT)),"DisplayName","nCHT")
hold off
xlim([0 1]);
ylim([y_low y_high]);
xlabel("Normalized Angular Frequency(\times\pi)   [rad/sample]");
ylabel("Magnitude   [dB] ");
legend("location", "southwest")

%% filtering & calucurate SNR
for j=(f_left_pass*fs2):2000:(f_right_pass*fs2)
p=p+1;
    
%% define signals

% define frequency
sig_f=j;
normalized_sig_f = sig_f/fs2;
sig_f_noize1=notch_f1*fs2;
sig_f_noize2=notch_f2*fs2;

% signals
sig=A*cos(2*pi*sig_f*t);
% sig_n=A*cos(2*pi*sig_f*t)+A_noize_rate_1*A*cos(2*pi*sig_f_noize1*t); 
% sig_n=A*cos(2*pi*sig_f*t)+A_noize_rate_2*A*cos(2*pi*sig_f_noize2*t);
sig_n=A*cos(2*pi*sig_f*t)+A_noize_rate_1*A*cos(2*pi*sig_f_noize1*t)+...
      A_noize_rate_2*A*cos(2*pi*sig_f_noize2*t); 

%% filtering

sig_1=[sig zeros(1,N/2)];
sig_n_1=[sig_n zeros(1,N/2)];

% filtering BPHT (w/ w/o noize)
sig_Im_BPHT=filter(h_BPHT,1,sig_1);
sig_n_Im_BPHT=filter(h_BPHT,1,sig_n_1);

% filtering nCHT (w/ w/o noize)
sig_Im_nCHT=filter(h_CHT,1,sig_1);
sig_n_Im_nCHT=filter(h_CHT,1,sig_n_1);

%% calucurate SNR

SNR_BPHT_out=snr(sig_Im_BPHT, fs);
SNR_noize_BPHT_out=snr(sig_n_Im_BPHT, fs);
SNR_nCHT_out=snr(sig_Im_nCHT, fs);
SNR_noize_nCHT_out=snr(sig_n_Im_nCHT, fs);

snr_outcome(p,1) = normalized_sig_f;
snr_outcome(p,2) = SNR_BPHT_out;
snr_outcome(p,3) = SNR_noize_BPHT_out;
snr_outcome(p,4) = SNR_nCHT_out;
snr_outcome(p,5) = SNR_noize_nCHT_out;

% if figures are needed, insert here
end

%% display SNR outcome
disp(snr_outcome)
disp('------------------------------------------------------')
disp('avarage')
SNR_outcome(:,2)=mean(snr_outcome(:,2),1);
SNR_outcome(:,4)=mean(snr_outcome(:,4),1);
SNR_outcome(:,3)=mean(snr_outcome(:,3),1);
SNR_outcome(:,5)=mean(snr_outcome(:,5),1);
disp(SNR_outcome)
disp('------------------------------------------------------')

%% extras
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