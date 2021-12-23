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

%% draw ideal amplitude characteristic
figure
hold on
grid off
xticks([0 0.1 0.2 0.25 0.3 0.7 0.8 1]);
xticklabels({'0','\omega_{b1}','\omega_{b2}','\omega_{ls}','\omega_{lp}','\omega_{rp}','\omega_{rs}','\pi'});
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
set(gca,'TickLength',[0 0]);
exportgraphics(gca,strcat('.\figure\ideal_amp','.pdf'),'ContentType','vector');