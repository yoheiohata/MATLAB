clear;close all;clc

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


%ノッチの場所
notch_freq_1 = 0.1;       
notch_freq_2 = 0.3;

Weight = [1 1 1];         

fs = 50000;
fs2 = fs/2;
L=15000;
ts=1/fs;
t=(0:L-1)*ts;
F=0:0.1:fs2;
a = [0 0 1 1 0 0];
limy_low = -100;
limy_high = 70;

%% ノッチフィルタ

%　ノッチフィルタ
nf1 = 3600*180/fs2;    r1=1;        %notch1 & on,off
nf2 = 7200*180/fs2;    r2=1;        %notch2 & on,off
nf3 = 3600*180/fs2;    r3=1;        %notch3 & on,off

D2=1;
D4=1;

notch_weight1=[1 1 1];
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

h_notch = load('notch3.cof');
[H_notch,w]=freqz(h_notch,1,F,fs);   


% notch filter だけ
figure;
x = [notch_freq_1 notch_freq_1];
y = [limy_low limy_high];
line(x,y,'Color','green','linestyle','--',"DisplayName","f="+notch_freq_1);
x = [notch_freq_2 notch_freq_2];
y = [limy_low limy_high];
line(x,y,'Color','green','linestyle','--',"DisplayName","f="+notch_freq_2);
legend("location","southeast")
hold on
plot(w/fs2, 20*log(abs(H_notch)),"DisplayName","Amplitude");
hold off
xlim([0 1]);
ylim([limy_low limy_high]);
xlabel("Normalized Angular Frequency(\times\pi)   [rad/sample]");
ylabel("Magnitude   [dB] ");
exportgraphics(gca,strcat('.\figure\notch','.pdf'),'ContentType','vector');

