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
N = 30;
f_left_stop  = 0.2;
f_left_pass  = 0.3;
f_right_pass = 0.7;
f_right_stop = 0.8;
notch_freq_1 = 0.04;
notch_freq_2 = 0.17;
Weight = [1 1 1];
fs = 4000;
fs2 = fs/2;
L=15000;
ts=1/fs;
t=(0:L-1)*ts;
F=0:0.1:fs2;
f = [0 f_left_stop f_left_pass f_right_pass f_right_stop 1];
a = [0 0 1 1 0 0];

y_low = -150;
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
notch2 = [1 -2*cos(notch_freq_2*pi) 1*1];
notch3 = conv(notch1, notch2);

notch1 = notch1/sum(notch1);
notch2 = notch2/sum(notch2);
notch3 = notch3/sum(notch3);
save notch1.cof notch1 -ascii;
save notch2.cof notch2 -ascii;
save notch3.cof notch3 -ascii;
fresp1 = {'fresp_h1',a};
fresp2 = {'fresp_h2',a};
fresp3 = {'fresp_h3',a};
h_notch = load('notch3.cof');
[H_notch,w]=freqz(h_notch,1,F,fs);   

% hilbert filter 
h_HT=firpm(N-4,f,fresp3,Weight,'hilbert');
[H_HT,w]=freqz(h_HT,1,F,fs);

% convolution
h_CHT=conv(h_notch,h_HT);
[H_CHT,w]=freqz(h_CHT,1,F,fs);

%% draw
figure;
x = [notch_freq_1 notch_freq_1];
y = [y_low y_high];
line(x,y,'Color','green','linestyle','--',"DisplayName","f(noise)="+notch_freq_1);
x = [notch_freq_2 notch_freq_2];
y = [y_low y_high];
line(x,y,'Color','green','linestyle','--',"DisplayName","f(noise)="+notch_freq_2);


hold on
% plot(w/fs2, 20*log(abs(H_notch)),"DisplayName","notch");
% plot(w/fs2, 20*log(abs(H_HT)),"DisplayName","ht")
plot(w/fs2, 20*log(abs(H_CHT)),"DisplayName","all")
hold off
xlim([0 1]);
ylim([y_low y_high]);
xlabel("Normalized Angular Frequency(\times\pi)   [rad/sample]");
ylabel("Magnitude   [dB] ");
legend("location", "northeast")
exportgraphics(gca,strcat('.\figure\amp_proposed_N=',num2str(N),'.pdf'),'ContentType','vector');

