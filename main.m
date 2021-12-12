clear;close all;clc;

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

noise_flag=0;%0:�m�C�Y�� 1:�m�C�Y�L


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% �@�t�B���^�쐻�@�@%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs = 50000;          % Sampling frequency;
fs2 = fs/2;
N=42; %�q���x���g�ϊ���̎���
N_low=10; %�q���x���g�ϊ���̎���(�Ⴂ���)
L=15000*1; %�M���� 15000
ts=1/fs; %�T���v�����O����
t=(0:L-1)*ts; %�T���v�����O����
F=0:0.1:fs2;     %�@���g���ݒ�
p_area1=0.05; %�J�ڗ̈�(f��0)
p_area2=1-p_area1; %�J�ڗ̈�(f��1)
sig_f=1800; %�����g�M���̎��g��
sig=sin(2*pi*sig_f*t);
sig=cat(2,sin(2*pi*sig_f*t(1:end/2)),sin(2*pi*2700*t(end/2+1:end))); %�����g�M��
% sig=sig+sin(2*pi*1600*t);

sig_noise=awgn(sig,60,'measured');

figure;plot(sig_noise);
xlim([0 150]);
ylim([-1.1 1.1]);

figure;plot(sig_noise-sig);
xlim([0 150]);
ylim([0.8 1.1]);

noise=sig_noise-sig;

%�m�C�Y�m�F�p
bunsan_sig=var(sig)
bunsan_noise=var(noise)
ave_sig=mean(sig)
ave_noise=mean(noise)
snr=snr(sig,noise)
jikkouthi_noise=sqrt(mean(noise.^2))
jikkouthi_sig=sqrt(mean(sig.^2))
jikkouthi_noise=rms(noise)
jikkouthi_sig=rms(sig)
SNR=20*log10(jikkouthi_sig/jikkouthi_noise)

sig_no_noise=sig;
if noise_flag==1
    sig=sig_noise;
end

figure;
plot(t,sig);
hold on ;
plot(t,sig_no_noise);
grid;
xlim([0 0.003]);
ylim([-1.2 1.2]);




%%
% 1st �t�B���^ ------------------
N1=12; %�t�B���^����
f1=[   0.0        0.0001       ... % �� 1 �ш�̉������g�� ������g���i�ʉ߈�j
       3000/fs2     1             ... % �� 2 �ш�̉������g�� ������g���i�j�~��j
   ];

%  2nd �t�B���^ ------------------
N2=28; %�t�B���^����
f2=[   0.0             0.0001        ... % �� 1 �ш�̉������g�� ������g���i�ʉ߈�j
       3000/fs2     1              ... % �� 2 �ш�̉������g�� ������g���i�j�~��j
       ];

% 3rd �t�B���^ ----------------------
N3=8; %�t�B���^����
f3=[   0.0            0.0001         ... % �� 1 �ш�̉������g�� ������g���i�ʉ߈�j
       3000/fs2    1               ... % �� 2 �ш�̉������g�� ������g���i�j�~��j
       ];
  
%% �m�b�`

nf1 = 3600*180/fs2;    r1=1;        %notch1 & on,off
nf2 = 7200*180/fs2;    r2=1;        %notch2 & on,off
nf3 = 3600*180/fs2;    r3=1;        %notch3 & on,off

%%
w1=[   1 ...                      % �� 1 �ш�̃E�G�C�g
    1 ...                         % �� 2 �ш�̃E�G�C�g
   ];

a1=[  1 1 ...                     % �p�X�o���h���z���� 
     0 0 ...                      % �X�g�b�v�o���h���z����
   ];


D2=1;
D4=1;
[notch1,notch2,notch3]=notch(nf1,nf2,nf3,D2,D4,fs2,r1,r2,r3);
notch3=conv(notch1,notch2);

Notch3=freqz(notch3,1,F,fs);

save notch1.cof notch1 -ascii;
save notch2.cof notch2 -ascii;
save notch3.cof notch3 -ascii;

%% �t�B���^�[����

fresp1 = {'fresp_h1',a1};      %�m�b�`1���l�������݌v
h1=firpm(N1,f1,[1 1 0 0 ],w1);
h1=h1./sum(h1);
hh1=conv(h1,notch1);
hh1=h1;
hh1=hh1./sum(hh1);
H1=freqz(hh1,1,F,fs);
save hh1.cof hh1 -ascii;
% 
fresp2 = {'fresp_h2',a1};      %1�i�ڂƃm�b�`2���l�������݌v
h2=firpm(N2,f2,[1 1 0 0],w1);
h2=h2./sum(h2);
hh2=conv(h2,notch3);      
hh2=hh2./sum(hh2);
H_notch3=freqz(notch3,1,F,fs/D2);
H2=freqz(h2,1,F,fs/D2);
H22=freqz(hh2,1,F,fs/D2);
save hh2.cof hh2 -ascii;


fresp3 = {'fresp_h3',a1};      %1�i�ڂ̍�����
h3=firpm(N3,f3,fresp3,w1);
h3=h3./sum(h3);
hh3=conv(h3,notch3);      
hh3=hh3./sum(hh3);
H3=freqz(hh3,1,F,fs/D4);
H3_after=freqz(h3,1,F,fs/D4);
H3_before=freqz(notch3,1,F,fs/D4);
save hh3.cof hh3 -ascii;


%% �t�B���^�̐U������
figure;
plot(F,abs(H1));
grid;
xlabel('Frequency [Hz]');
ylabel('Amplitude')
hold on;
plot(F,abs(H3));

figure;
plot(F,abs(H1.*H2));
hold on;
plot(F,abs(H3));
grid;
xlabel('Frequency [Hz]');
ylabel('Amplitude')
legend('�c���ڑ�','�_�u���m�b�`')

figure;
plot(F,abs(H3));
grid;
xlabel('Frequency [Hz]');
ylabel('Amplitude')
title('�t�B���^�̐U������')

figure;
plot(F,abs(H3));
hold on;
plot(F,abs(H3_before));
plot(F,abs(H3_after));
grid;
ylim([0 1]);
xlabel('Frequency [Hz]');
ylabel('Amplitude')
legend('total','1st stage','2nd stage');


figure;
plot(F,20*log10(abs(H22)),'LineWidth',2);
hold on;
plot(F,20*log10(abs(H3_before)),'--','LineWidth',2);
plot(F,20*log10(abs(H2)),':','LineWidth',2.1);
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
grid on;
ylim([-60 40]);
legend('total(N=32)','1st step(N=4)','2nd step(N=28)')
set(gca,'FontSize',16);

figure;
plot(F,20*log10(abs(H3)),'LineWidth',2);
hold on;
plot(F,20*log10(abs(H3_before)),'LineWidth',2);
plot(F,20*log10(abs(H3_after)),'LineWidth',2);
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
grid on;
ylim([-60 40]);
legend('total(N=12)','1st step(N=4)','2nd step(N=8)')
set(gca,'FontSize',16);

figure;
plot(F,20*log10(abs(H3_before)),'color',color2,'LineWidth',2);
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
grid on;
ylim([-60 40]);
legend('total(N=12)','1st step(N=4)','2nd step(N=8)')
set(gca,'FontSize',16);

figure;
plot(F,20*log10(abs(H3)),'--','LineWidth',1.4);
plot(F,20*log10(abs(H3_before)),'--','color',color2,'LineWidth',1.4);
hold on;
plot(F,20*log10(abs(H3_after)),'color',color3,'LineWidth',2);
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
grid on;
ylim([-60 40]);
legend('1st step(N=4)','2nd step(N=8)')
set(gca,'FontSize',16);

figure;
plot(F,20*log10(abs(H3)),'LineWidth',2);
hold on;
plot(F,20*log10(abs(H3_before)),'--','LineWidth',1.4);
plot(F,20*log10(abs(H3_after)),'--','LineWidth',1.4);
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
grid on;
ylim([-60 40]);
legend('total(N=12)','1st step(N=4)','2nd step(N=8)')
set(gca,'FontSize',16);

figure;
plot(F,20*log10(abs(H3)),'--','LineWidth',1.4);
hold on;
plot(F,20*log10(abs(H3_before)),'LineWidth',2);
plot(F,20*log10(abs(H3_after)),'--','LineWidth',1.4);
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
grid on;
ylim([-60 40]);
legend('total(N=12)','1st step(N=4)','2nd step(N=8)')
set(gca,'FontSize',16);

figure;
plot(F,20*log10(abs(H3)),'LineWidth',2);
hold on;
plot(F,20*log10(abs(H3_before)),'-','LineWidth',2);
plot(F,20*log10(abs(H3_after)),':','LineWidth',2.1);
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
grid on;
ylim([-60 40]);
legend('total(N=12)','1st step(N=4)','2nd step(N=8)')
set(gca,'FontSize',16);


fvtool(hh3,1,hh2,1);
% ylabel ���쐬
 ylabel('Magnitude (dB)');
% xlabel ���쐬
 xlabel('Frequency (kHz)');
 legend('�O�i�l������','�O�i�l���Ȃ�');

%% �u�����g��
[h,err,res]=firpm(N,[p_area1 p_area2],[1 1],{10000},'Hilbert'); %�t�B���^�W���̎Z�o
[H,W]=freqz(h,1,F,fs);
[h60,err60,res60]=firpm(60,[p_area1 p_area2],[1 1],{10000},'Hilbert'); %�t�B���^�W���̎Z�o
[H60,W60]=freqz(h60,1,F,fs);

figure;
plot(W/fs2,abs(H),'LineWidth',2);
hold on;
plot(W60/fs2,abs(H60),'LineWidth',2);
grid;
ylabel('Magnitude');
xlabel('Normalized Frequency (\times\pi rad/sample)');
legend('N=42','N=60')
set(gca,'FontSize',16);
ylim([0 1.05]);

figure;
plot(W/fs2,abs(H),'LineWidth',2);
hold on;
plot(W60/fs2,abs(H60),'LineWidth',2);
grid;
ylabel('Magnitude');
xlabel('Normalized Frequency (\times\pi rad/sample)');
legend('N=42','N=60')
set(gca,'FontSize',16);
ylim([0.95 1.05]);
xlim([0.0 0.3]);


%%�x���␳
sig_1=[sig zeros(1,N/2)]; 
%%%%
sig_Im=filter(h,1,sig_1); %�q���x���g�ϊ�
sig_Re=sig_1(1:end-N/2); %����
sig_analysis=sig_Re+1i*sig_Im((N/2+1):end); %��͐M��
sig_analysis=sig_analysis(100:14000);
t=t(100:14000);

figure;
plot(sig);
hold on;
plot(sig_Im((N/2+1):end));
xlim([1 120]);

%%�u�����g���̓��o%%
phase=angle(sig_analysis); %��͐M���̈ʑ��p���v�Z
phase=unwrap(phase); %���炩�ɂȂ�悤�ʑ��p�̏C��
omega_tilde=diff(phase)/ts; %�ʑ��p�̔���
IF=omega_tilde/(2*pi); %�u�����g��


plotfft(IF,fs); %�u�����g���̃X�y�N�g����\��

%% �u�����g��Low
[h_low,err_low,res_low]=firpm(N_low,[p_area1 p_area2],[1 1],{10000},'Hilbert'); %�t�B���^�W���̎Z�o
% ylabel ���쐬
ylabel('Magnitude');
% xlabel ���쐬
xlabel('Normalized Frequency (\times\pi rad/sample)');

%%�x���␳
sig_low=[sig zeros(1,N_low/2)]; 
%%%%
sig_Im_low=filter(h_low,1,sig_low); %�q���x���g�ϊ�
sig_Re_low=sig_low(1:end-N_low/2); %����
sig_analysis_low=sig_Re_low+1i*sig_Im_low((N_low/2+1):end); %��͐M��
figure;
plot(sig);
hold on;
plot(sig_Im((N/2+1):end));
xlim([1 120]);

%%�u�����g���̓��o%%
phase_low=angle(sig_analysis_low); %��͐M���̈ʑ��p���v�Z
phase_low=unwrap(phase_low); %���炩�ɂȂ�悤�ʑ��p�̏C��
omega_tilde_low=diff(phase_low)/ts; %�ʑ��p�̔���
IF_low=omega_tilde_low/(2*pi); %�u�����g��
 plotfft(IF,fs); %�u�����g���̃X�y�N�g����\��

%% �t�B���^��ʂ�
IF_1=filter(hh1,1,IF); %�᎟��1�i��
IF_2=filter(hh2,1,IF_1); %�᎟��2�i��

IF_3=filter(hh3,1,IF); %������1�i��
IF_3=filter(notch3,1,IF_3); %������1�i��



%low
IF_1_low=filter(hh1,1,IF_low); %�᎟��1�i��
IF_2=filter(hh2,1,IF_1); %�᎟��2�i��
IF_3_low=filter(hh3,1,IF_low); %������1�i��


%% ��r�pNF
%%�U��������(NF)%%
r2=0.4;
w20=0; %NF��w0�̏����l(w(n))
w21=0; %NF��w1�̏����l(w(n-1))
w22=0; %NF��w2�̏����l(w(n-2))
e_2(1:2)=IF(1:2); %�����pNF�̏o��
e_3(1:2)=e_2(1:2);
 for i=2:length(t)-1
    a2=-(1+r2)*cos(2*pi*sig_f*2/fs); %NF�̃t�B���^�W������
    x2=IF(i); %NF�̓���
    w22=w21;
    w21=w20; %��ԕϐ��̍X�V
    w20=x2-(a2*w21+r2*w22);
    e_2(i)=0.5*(x2+r2*w20+a2*w21+w22); %NF�̏o��
 end
 
  for i=2:length(t)-1
    a2=-(1+r2)*cos(2*pi*sig_f*4/fs); %NF�̃t�B���^�W������
    x3=e_2(i); %NF�̓���
    w22=w21;
    w21=w20; %��ԕϐ��̍X�V
    w20=x3-(a2*w21+r2*w22);
    e_3(i)=0.5*(x3+r2*w20+a2*w21+w22); %NF�̏o��
  end
 
%NF�̐U�������m�F  
w=linspace(0,pi,L-1);
f_w=(w/(2*pi))*fs;
H_z_a=[1+r2 a2*2 1+r2]; %���q
H_z_b=2.*[1 a2 r2]; %����
H_z=freqz(H_z_a,H_z_b,L-1);
mag_H_z=abs(H_z);
figure;
plot(f_w,20*log10(mag_H_z),'LineWidth',2)
set(gca,'FontSize',16);
grid;
xlabel('Frequency[Hz]')
ylabel('Magunitude[dB]')

%NF�Ƃ̐��\��r
figure;
set(gca,'FontSize',16);
plot(t(1:end-1),IF,'LineWidth',1);
hold on;
plot(t(1:end-1),IF_3,'LineWidth',2);
ylim([1740 1860]);
xlim([0 0.002]);
xlabel("Time[s]");ylabel("Frequency[Hz]");
grid;
plot(t(1:end-1),e_2,'--','LineWidth',2);
plot(t(1:end-1),e_3,':','LineWidth',2);
legend('IF','FIR','NF','NF2');

figure;
set(gca,'FontSize',16);
plot(IF,'LineWidth',1);
hold on;
plot(IF_3,'LineWidth',2);
ylim([1799 1801]);
xlim([0 55]);
xlabel("Time[s]");ylabel("Frequency[Hz]");
grid;
plot(e_2,'--','LineWidth',2);
plot(e_3,':','LineWidth',2);
legend('IF','FIR','NF','NF2');

%% �����x�̐��x�̂��߂̃q���x���g�ϊ��펟������
middle_sig=IF_3(5000:10000);
delta_max=max(middle_sig)-sig_f;
delta_min=min(middle_sig)-sig_f;
figure;
plot(middle_sig)
delta=max(abs(delta_max),abs(delta_min));

delta_now=delta+1;

N=200-2;
while delta_now > delta
    N=N+2

    sig=sin(2*pi*sig_f*t); %�����g�M��
  
    h=firpm(N,[p_area1 p_area2],[1 1],'Hilbert'); %�t�B���^�W���̎Z�o
    fvtool(h,1); %�U�������m�F

    %%���p�r�s��%%(�x���̂��߁H)
    sig_3=[sig zeros(1,N/2)]; 
    %%%%
    sig_Im=filter(h,1,sig_3); %�q���x���g�ϊ�
    sig_Re=sig_3(1:end-N/2); %����
    sig_analysis=sig_Re+1i*sig_Im((N/2+1):end); %��͐M��
    %%�������̑O�����������ƂŒx�������l�����Ă���(?)

    %%�u�����g���̓��o%%
    phase=angle(sig_analysis); %��͐M���̈ʑ��p���v�Z
    phase=unwrap(phase); %���炩�ɂȂ�悤�ʑ��p�̏C��
    delta_phase=diff(phase); %�ʑ��p�̔���
    IF_loop=delta_phase/(2*pi)*fs; %�u�����g��

    middle_sig=IF_loop(5000:10000);
    delta_max=max(middle_sig)-sig_f;
    delta_min=min(middle_sig)-sig_f;
    delta_now=max(abs(delta_max),abs(delta_min));
end


%% ���ʏo��
figure;
plot(t(1:end-1),IF_3,'LineWidth',2,'color',color2); %�u�����g���̊m�F
hold on;
plot(t(1:end-1),IF_loop,'--','LineWidth',2.4,'color',color3);
plot(t(1:end-1),e_3,':','LineWidth',1.7,'color',color4);
plot(t(1:end-1),e_2,'-.','LineWidth',2);
ylim([1798 1802]);
ylim([1780 1830]);
xlim([0.0 0.003]);
legend('proposed','instantaneous frequency(N=200)','notch filter(f_n=3,600and7,200)')
xlabel("Time[s]");ylabel("Frequency[Hz]");
pbaspect([3 2 1]);
grid;
set(gca,'FontSize',16);


%% �G������
figure;
plot(t(1:end-1),IF,'LineWidth',2,'color',color1); %�u�����g���̊m�F
hold on;
plot(t(1:end-1),IF_3,'LineWidth',2,'color',color2); %�u�����g���̊m�F
plot(t(1:end-1),e_3,':','LineWidth',1.7,'color',color3);
ylim([1799.9995 1800.001]);
xlim([0.0 0.002]);
legend('IF','proposed','notch filter(f_n=3,600and7,200)')
xlabel("Time[s]");ylabel("Frequency[Hz]");
grid;
set(gca,'FontSize',16);

figure;
plot(t(1:end-1),IF,'LineWidth',2,'color',color1); %�u�����g���̊m�F
hold on;
plot(t(1:end-1),IF_3,'LineWidth',2,'color',color2); %�u�����g���̊m�F
plot(t(1:end-1),IF_loop,'-','LineWidth',2,'color',color3);
ylim([1799.9999 1800.0002]);
xlim([0.0 0.003]);
% pbaspect([3 2 1]);
legend('without filter(N=42)','proposed','without filter(N=200)')
xlabel("Time[s]");ylabel("Frequency[Hz]");
grid on;
set(gca,'FontSize',16);
 
figure;
plot(t(1:end-1),IF,'LineWidth',2,'color',color1); %�u�����g���̊m�F
hold on;
plot(t(1:end-1),IF_3,'-','LineWidth',2,'color',color2); %�u�����g���̊m�F
plot(t(1:end-1),IF_loop,'-','LineWidth',2,'color',color3);
ylim([1780 1830]);
xlim([0.0 0.003]);
grid on
legend('without filter(N=42)','proposed','without filter(N=200)')
xlabel("Time[s]");ylabel("Frequency[Hz]");
set(gca,'FontSize',16);

 %% �G������
figure;
plot(t(1:end-1),IF,'LineWidth',2,'color',color1); %�u�����g���̊m�F
hold on;
plot(t(1:end-1),IF_3,'LineWidth',2,'color',color2); %�u�����g���̊m�F
plot(t(1:end-1),e_3,':','LineWidth',1.7,'color',color3);
ylim([1790 1820]);
xlim([0.0 0.002]);
legend('IF','proposed','notch filter(f_n=3,600and7,200)')
xlabel("Time[s]");ylabel("Frequency[Hz]");
grid;
set(gca,'FontSize',16);

 
figure;
h1=plot(t(1:end-1),IF,'LineWidth',2,'color',color1); %�u�����g���̊m�F
hold on;
h3=plot(t(1:end-1),IF_loop,'-','LineWidth',2,'color',color3);
h2=plot(t(1:end-1),IF_3,'-','LineWidth',2,'color',color2); %�u�����g���̊m�F
ylim([1780 1830]);
xlim([0.0 0.003]);
grid on
legend([h1,h2,h3],{'HT only(N=42)','proposed(N=42+12)','HT only(N=200)'});
xlabel("Time[s]");ylabel("Frequency[Hz]");
set(gca,'FontSize',16);

figure;
h1=plot(t(1:end-1),IF,'LineWidth',2,'color',color1); %�u�����g���̊m�F
hold on;
h3=plot(t(1:end-1),IF_loop,'-','LineWidth',2,'color',color3);
h2=plot(t(1:end-1),IF_3,'-','LineWidth',2,'color',color2); %�u�����g���̊m�F
ylim([1799.9999 1800.0002]);
xlim([0.0 0.003]);
grid on
legend([h1,h2,h3],{'HT only(N=42)','proposed(N=42+12)','HT only(N=200)'});
xlabel("Time[s]");ylabel("Frequency[Hz]");
set(gca,'FontSize',16);


plotfft(IF_3,fs);
plotfft(IF_3_low,fs);

figure;
grpdelay(hh3,1);

figure;
grpdelay(notch3,1)

figure;
grpdelay(h3,1)


Fs = 500;
N = 500;
rng default
xn = ecg(N)+0.1*randn([1 N]);
xn=zeros(1,N);
xn(1)=1;

Nfir = 70;
Fst = 75;

firf = designfilt('lowpassfir','FilterOrder',Nfir, ...
    'CutoffFrequency',Fst,'SampleRate',Fs);
xf = filter(firf,xn);

figure;
plot(xn);
hold on;
plot(xf);
title 'Electrocardiogram'
xlabel 'Time (s)'
legend('Original','FIR Filtered')
grid;
xlim([0 100]);


test1=filter(h3,1,xn);
test2=filter(notch3,1,xn);
test3=filter(hh3,1,xn);

figure;
plot(test1);
hold on;
plot(test2);
plot(test3);
title 'Electrocardiogram'
xlabel 'Time (s)'
legend('Original','FIR Filtered')
grid;
xlim([0 100]);
ylim([-1 1]);
