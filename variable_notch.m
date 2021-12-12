clear;
run('C:\Users\kt33_\Box\相川研\個人\30期生\高尾\MATLAB\my_init.m')
% close all;
% Minimize stopband ripple of a linear phase lowpass FIR filter
% "Filter design" lecture notes (EE364) by S. Boyd
% (figures are generated)
%
% Designs a linear phase FIR lowpass filter such that it:
% - minimizes the maximum passband ripple
% - has a constraint on the maximum stopband attenuation
%
% This is a convex problem.
%
%   minimize   delta
%       s.t.   1/delta <= H(w) <= delta     for w in the passband
%              |H(w)| <= atten_level        for w in the stopband
%
% where H is the frequency response function and variables are
% delta and h (the filter impulse response).
%
% Written for CVX by Almir Mutapcic 02/02/06

%********************************************************************
% user's filter specifications
%********************************************************************
% filter order is 2n+1 (symmetric around the half-point)
n_all = 20;

wpass = 0.2*pi;        % passband cutoff freq (in radians)
wstop = 0.35*pi;        % stopband start freq (in radians)
% atten_level = -30;      % stopband attenuation level in dB

N_of_Notch=2;    %ノッチの数
alpha_range=linspace(0.2,0.4,20); % 設計時の可変範囲
% beta_range=linspace(0,1,length(alpha_range));
L1=5;%多項式の次数


n=(n_all-2*N_of_Notch)/2;

%********************************************************************
% create optimization parameters
%********************************************************************
N = 30*n+1;                            % freq samples (rule-of-thumb)
w = linspace(0,pi,N);
A = [cos(kron(w',repelem([0:n],L1+1)))]; % matrix of cosines

Ap=[];
As=[];
L_range=[0:L1];

for q=1:length(alpha_range)
    %重み(ノッチフィルタ部分)の作成
    weight_notch_filter=ones(N,1);
%     notch_coef=[1];
    for i=1:N_of_Notch;
        my_theta=alpha_range(q)*pi*i*2;
        % my_theta=calculation_theta(my_theta,mode_flag);
        weight_notch_filter = weight_notch_filter.*[2*(cos(w')-cos(my_theta))];
%         notch_coef_temp = [1 -2*cos(my_theta) 1];
%         notch_coef=conv(notch_coef_temp,notch_coef);
    end
    %   figure;
    %   plot(w/pi,abs(weight_notch_filter));
    
    beta=((alpha_range(q)-alpha_range(1))/(alpha_range(end)-alpha_range(1)));
    beta_box=beta.^L_range;
    beta_box=repmat(beta_box,1,n+1);
    
    A_q=A.*weight_notch_filter;
    A_q=A_q.*beta_box;
    
    % passband 0 <= w <= w_pass
    ind = find((0 <= w) & (w <= wpass));   % passband
    Ap  = [Ap;A_q(ind,:)];
    
    % transition band is not constrained (w_pass <= w <= w_stop)
    
    % stopband (w_stop <= w)
    ind = find((wstop <= w) & (w <= pi));  % stopband
    %   Us  = 10^(atten_level/20)*ones(length(ind),1);
    As  = [As;A_q(ind,:)];
end

%********************************************************************
% optimization
%********************************************************************
% formulate and solve the linear-phase lowpass filter design
cvx_begin
    variable delta
    variable g((n+1)*(L1+1));

    minimize( delta )
    subject to
    % passband bounds
    Ap*g <= delta+1;
    inv_pos(Ap*g) <= delta+1;

    % stopband bounds
    abs( As*g ) <= delta;

cvx_end

% check if problem was successfully solved
disp(['Problem is ' cvx_status])
if ~strfind(cvx_status,'Solved')
    return
else
    % construct the full impulse response
    %   h = [flipud(h(2:end)); h];
    fprintf(1,'The optimal minimum passband ripple is %4.3f dB.\n\n',...
        20*log10(delta));
end



%********************************************************************
% fripm
%********************************************************************

h1_LP=firpm(n_all,[0 wpass/pi wstop/pi 1],[1 1 0 0],[1 1]);
% h1_LP_comp=h1_LP_comp./sum(h1_LP_comp);
H1_LP=freqz(h1_LP,1,w);



%********************************************************************
% plots
%********************************************************************
% figure(1)
% % FIR impulse response
% plot([0:2*n],h','o',[0:2*n],h','b:')
% xlabel('t'), ylabel('h(t)')
figure(1000)    
hold off;
plot([0 wpass/pi],[0 0],'k--');
hold on;
plot([wstop/pi 1],[20*log10(max(delta)) 20*log10(max(delta))],'k--');
axis([0,1,-50,5])
xlabel('w'), ylabel('mag H(w) in dB')

alpha_check=linspace(0.2,0.4,111);
for i=1:length(alpha_check)
    beta=((alpha_check(i)-alpha_range(1))/(alpha_range(end)-alpha_range(1)));
    h=get_filter_coef(g,beta,n,L1);
    
    h_notch=[1];
    for j=1:N_of_Notch;
        my_theta=alpha_check(i)*pi*j*2;
        % my_theta=calculation_theta(my_theta,mode_flag);
%         weight_notch_filter = weight_notch_filter.*[2*(cos(w')-cos(my_theta))];
        notch_coef_temp = [1 -2*cos(my_theta) 1];
        h_notch=conv(notch_coef_temp,h_notch);
    end
    H_notch=freqz(h_notch,1,w);
    H_temp=freqz(h,1,w);
    
    h_all=conv(h,h_notch);
    H=freqz(h_all,1,w);
  
    
    figure(50+i)
    hold off;
    
    % magnitude
    plot(w/pi,20*log10(abs(H)),'DisplayName',['VMX']);
    hold on;
    plot(w/pi,20*log10(abs(H1_LP)),'DisplayName',['Remez']);
    plot(w/pi,20*log10(abs(H_notch)),'DisplayName',['notch']);
    plot([0 wpass/pi],[0 0],'k--');
    plot([wstop/pi 1],[20*log10(delta) 20*log10(delta)],'k--');
    axis([0,1,-50,5])
    xlabel('w'), ylabel('mag H(w) in dB')
    legend;
    
%         figure(150+i)
%         hold off;
%         % magnitude
%         plot(w/pi,(abs(H)),'DisplayName',['VMX']);
%         hold on;
%         plot(w/pi,(abs(H1_LP)),'DisplayName',['Remez']);
%         axis([0,1,0,1.1])
%         legend;
        % xlabel('w'), ylabel('mag H(w) in dB')
    
    % figure(2)
    % % frequency response
    % H = exp(-j*kron(w',[0:2*n]))*h;
    % % magnitude
    % subplot(2,1,1)
    % plot(w,20*log10(abs(H)),[wstop pi],[atten_level atten_level],'r--');
    % axis([0,pi,-40,10])
    % xlabel('w'), ylabel('mag H(w) in dB')
    % % phase
    % subplot(2,1,2)
    % plot(w,angle(H))
    % axis([0,pi,-pi,pi])
    % xlabel('w'), ylabel('phase H(w)')
    
    figure(1000)    
    % magnitude
        hold on;
    plot(w/pi,20*log10(abs(H)),'DisplayName',['']);

%     legend;
end