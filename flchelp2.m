%function [h,edge] = flchelp2(M,L,ce)
% [h,edge] = flchelp2(M,L,es)
% A program for the design of lowpass FIR filters with maximally 
% flat passband and equiripple stopband characteristics.
% This program can treat Type1 and Type2 filters.
%  Type1 : filter length is  odd, impluse response is symmetric
%  Type2 : filter length is even, impluse response is symmetric
%
% require subprograms : local_max.m,taylor1.m
%
% <parameters>
%  M  : filter order
%  L  : degree of flatness at w=0  ( L < M/2 )
%  ce : Chebyshev error in stopband
%
% Author: Jun Hashimoto, Tokyo Metropolitan University
% Update: Naoyuki Aikawa, Tokyo Metropolitan University
% Original program author: Ivan W.Selesnick, Rice University
%
% 論文
%   相川　直幸，佐藤　正光
% 　通過域平たんかつ阻止域等リプル特性を有する直線位相FIRディジタルフィルタの設計
%   電子情報通信学会論文誌(A) , Vol. J83-A No. 6.，pp744-749.，(2000,6)
%
%
% <example>
     M  =  21;
     L  =  6;
     ce =  0.05;
     wp = 0;                                     % passband center frequency
% ---------------------------------------------------------------------------------------
PF = 0;                                          % flag : Plot Figures
if L >= M/2
   fprintf(1,'\n    L must be less than M/2 \n')
   return
end
if rem(M,2) == 0
   ftype = 1;
else
   ftype = 2;
   M = M - 1;
end
N = M + 1;                                       % filter length
g = 2^ceil(log2(20*N));                          % number of grid points
SN = 1e-9;                                       % SN : SMALL NUMBER
q  = (N-2*L+1)/2;                                % number of frequency points
% ---------- initialize some constants --------------------------------------------------
Y  = [-ce*(1-(-1).^(1:q))/2 + ce*((-1).^(1:q)+1)/2]';
w  = [0:g]'*pi/g;                                % frequency axis
wo = pi*L/(L+q);
wo = wo - 2*ce;
if ftype == 1
   rs = [0:q - 1]'*(pi-wo)/(q-1) + wo;
else
   rs = [0:q]'*(pi-wo)/q + wo;
   rs(length(rs)) = [];
end
Z  = zeros(2*(g-q)+1,1);
n  = 0:q-1;
if ftype == 1
   B1 = 1;
   B1r = 1;
   B2 = 1;
   B2r = 1;
else
   B1 = cos(w/2);
   B1r = cos(rs/2);
   B2 = 0;
   B2r = 0;
   for i = 0:L-1
      B2 = B2 + taylor1(i,wp)*((1-cos(w))/2).^i;
      B2r = B2r + taylor1(i,wp)*((1-cos(rs))/2).^i;
   end
end
A1 = (-1)^L * sin(w/2).^(2*L);
A1r = (-1)^L * sin(rs/2).^(2*L);
it = 0;
while 1 & (it < 10)
   % ---------- calculate interpolation values ------------------------------------------
   a2 = cos(rs*n)\[((Y./B1r)-B2r)./A1r];
   A2 = real(fft([a2(1);a2(2:q)/2;Z;a2(q:-1:2)/2])); A2 = A2(1:g+1);
   A = B1.*(B2 + A1.*A2);
   % ---------- plot figures if PF == 1 -------------------------------------------------
   if PF
      figure(2)
      subplot(2,1,1), plot(w/pi,A), axis([0 1 -0.2 1.05]),
         hold on, plot(rs/pi,Y,'o'), hold off,
         xlabel('Normalized Angular Frequency'), ylabel('Amplitude')
      subplot(2,1,2), plot(w/pi,20*log10(abs(A))), axis([0 1 -100 5]),
         hold on, plot(rs/pi,20*log10(abs(Y)),'o'), hold off,
         xlabel('Normalized Angular Frequency'), ylabel('Amplitude [dB]'), 
      pause(1)
   end
   % ---------- determine new reference set ---------------------------------------------
   ri = sort([local_max(A); local_max(-A)]);
   if ftype == 1
      ri(1:length(ri)-q) = [];
      rs = (ri-1)*pi/g;
      A1r = (-1)^L * sin(rs/2).^(2*L);
      Ar = 1 + (cos(rs*n)*a2) .* A1r;
   else
      ri(1:length(ri)-(q+1)) = [];
      ri(length(ri)) = [];
      rs = (ri-1)*pi/g;
      A1r = (-1)^L * sin(rs/2).^(2*L);
      B1r = cos(rs/2);
      B2r = 1;
      for i = 1:L-1                              % calculate Taylor's series coeffifients
         B2r = B2r + taylor1(i,wp)*((1-cos(rs))/2).^i;
      end
      Ar = B1r.*(B2r + (cos(rs*n)*a2) .* A1r);
   end
   % ---------- calculate maximum constraint violation ----------------------------------
   Err = max([max(Ar)-ce, -ce-min(Ar)]);
   if PF
      fprintf(1,'    Err = %18.15f\n',Err);
   end
   if Err < SN
      break
   end
   it = it + 1;
end
% ---------- calculate filter coefficients and roots-------------------------------------
ha2 = [a2(q:-1:2)/2; a2(1); a2(2:q)/2];
h = ha2;
for k = 1:L
   h = conv(h,[1 -2 1]')/4;
end
if ftype == 1
   h((N+1)/2) = h((N+1)/2) + 1;
else
   hb2 = zeros(N,1); 
   for j = 0:L-1;
      x = [-0.25 0.5 -0.25];
      k = 1;
      for i = 1:j
         k = conv(x,k);
      end
      c = taylor1(j,wp) * k';
      hb2 = hb2 + [zeros((N-((j*2)+1))/2,1); c; zeros((N-((j*2)+1))/2,1)];
   end
   h = h + hb2;
   h = conv(h,[0.5 0.5]);
end
% ---------- calculate band edge --------------------------------------------------------
k = 1;
while 1
   j = A(ri(1)-k) - ce;
   if j > 0
      ew = ri(1) - k;
      break
   end
   k = k + 1;
end
edge = w(ew)/pi;
fprintf(1,'    Stopband edge             =  %8.5f\n',edge)
c = roots(h);                                    % calculate roots
figure(1) 
subplot(2,2,1), plot(w/pi,A), axis([0 1 -0.2 1.05]),
   hold on, plot(rs/pi,Y,'o'), hold off,
   xlabel('Normalized Angular Frequency'), ylabel('Amplitude')
subplot(2,2,3), plot(w/pi,20*log10(abs(A))), axis([0 1 -100 5]),
   hold on, plot(rs/pi,20*log10(abs(Y)),'o'), hold off,
   xlabel('Normalized Angular Frequency'), ylabel('Amplitude [dB]'), 
subplot(2,2,2), impz(h,1)
   xlabel('Impulse Response')
subplot(2,2,4), zplane(c)