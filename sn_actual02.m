clear;close all;clc

set(0,'defaultTextFontName','Times New Roman');
set(0,'defaultLegendFontName','Times New Roman');
set(0,'defaultAxesFontName','times new roman');
set(0,'defaultAxesGridLineStyle',':');
set(0,'defaultAxesBox','on');
set(0,'defaultAxesYGrid','on');
set(0,'defaultAxesXGrid','on');
set(0,'defaultAxesFontSize',16);
set(0,'DefaultAxesLineWidth', 1, 'DefaultLineLineWidth', 2); 

%% params
% filter params
N = 28;
f_left_stop  = 0.11;
f_left_pass  = 0.12;
f_right_pass = 0.8;
f_right_stop = 0.81;
notch_freq_1 = 0.07;       
Weight = [1 1 1];
fs = 4000;
fs2 = fs/2;
L=15000;
ts=1/fs;
t=(0:L-1)*ts;
F=0:0.1:fs2;
f = [0 f_left_stop f_left_pass f_right_pass f_right_stop 1];
a = [0 0 1 1 0 0];

% signal params
A=1;
A_noize_rate = 4.737639585013702 / 21.902957074846416;

p=0;
q=0;
r=0;
y_low = -100;
y_high = 20;

%% design filter

% notch filter
nf1 = 3600*180/fs2;    r1=1;    D2=1;
nf2 = 7200*180/fs2;    r2=1;    D4=1;
nf3 = 3600*180/fs2;    r3=1;
notch_weight1=[1 1 1];
notch_amplitude1=[0 0 1 1 0 0];
[notch1,notch2,notch3]=notch(nf1,nf2,nf3,D2,D4,fs2,r1,r2,r3);
notch1 = [1 -2*cos(notch_freq_1*pi) 1*1];
notch1 = notch1/sum(notch1);
save notch1.cof notch1 -ascii;
fresp3 = {'fresp_h1',a};
h_notch = load('notch1.cof');
[H_notch,w]=freqz(h_notch,1,F,fs);   

% hilbert filter 
h_HT=firpm(N,f,fresp3,Weight,'hilbert');
[H_HT,w]=freqz(h_HT,1,F,fs);

% convolution
h_CHT=conv(h_notch,h_HT);
[H_CHT,w]=freqz(h_CHT,1,F,fs);


% conventional
[h_HT1,err,res]=firpm(N+2,f,a,{10000},Weight);
[H_HT1,w]=freqz(h_HT1,1,F,fs);


% draw
figure; %nCHT
x = [notch_freq_1 notch_freq_1];
y = [y_low y_high];
line(x,y,'Color','green','linestyle','--',"DisplayName","noize="+notch_freq_1);

hold on
plot(w/fs2, 20*log(abs(H_notch)),"DisplayName","notch");
plot(w/fs2, 20*log(abs(H_HT)),"DisplayName","HT")
plot(w/fs2, 20*log(abs(H_CHT)),"DisplayName","all")
hold off
xlim([0 1]);
ylim([y_low y_high]);
legend("location", "southeast")
exportgraphics(gca,strcat('.\figure\nCHT\N=',num2str(N+2),'f=',num2str(f),'w=',num2str(Weight),'.pdf'),'ContentType','vector');


figure; % BPHT vs nCHT
x = [notch_freq_1 notch_freq_1];
y = [y_low y_high];
line(x,y,'Color','green','linestyle','--',"DisplayName","noize="+notch_freq_1);

hold on
plot(w/fs2, 20*log(abs(H_HT1)),"DisplayName", "conventional");
plot(w/fs2, 20*log(abs(H_CHT)),"DisplayName", "proposed");
hold off
xlim([0 1]);
ylim([y_low y_high]);
legend("location", "southeast")
exportgraphics(gca,strcat('.\figure\nCHTvsBPHT\N=',num2str(N+2),'w=',num2str(Weight),'.pdf'),'ContentType','vector');
%% assess signal

% define signal
sig_f = 296.875;
sig_f_n = 140.625;

sig = A*cos(2*pi*sig_f*t);
sig_n = sig + A_noize_rate*A*cos(2*pi*sig_f_n*t);

% out
sig_n_zeros = [sig_n zeros(1,N/2)];
sig_n_Re = sig_n_zeros(1:end-N/2);

sig_n_Im_nCHT = filter(h_CHT,1,sig_n_zeros);
sig_n_Im_BPHT = filter(h_HT1,1,sig_n_zeros);

sig_anly_CHT = sig_n_Re + 1i*sig_n_Im_nCHT((N/2+1):end);
sig_anly_BPHT = sig_n_Re + 1i*sig_n_Im_BPHT((N/2+1):end);

% cal SNrate
SN_BPHT = snr(filter(h_HT1,1,sig_n_zeros));
SN_nCHT = snr(filter(h_CHT,1,sig_n_zeros));


% figure;
% plot(t, sig_n);
% hold on
% plot(t, sig_n_Im_nCHT((N/2+1):end));
% plot(t, sig_n_Im_BPHT((N/2+1):end));
% hold off
% xlim([0.01 0.02]);
% ylim([-1.5 1.5])
% legend()
