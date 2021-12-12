%function [h,edge] = flchebp2(M,L,wp,ce1,ce2,itype)
% 等リプルに補正する方
% [h,edge] = flchebp2(M,L,wp,ce1,ce2,itype)
% A program for the design of bandpass FIR filters with maximally 
% flat passband and equiripple stopband characteristics.
% This program can treat all types of the linear phase FIR filters.
%  Type1 : filter length is  odd, impluse response is symmetric
%  Type2 : filter length is even, impluse response is symmetric
%  Type3 : filter length is  odd, impluse response is anti-symmetric
%  Type4 : filter length is even, impluse response is anti-symmetric
% 
% require subprograms : local_max.m,taylor1.m,taylor2
%
% <parameters>
%  M  : filter order
%  L  : degree of flatness at w=wo  ( L < M/2 )
%  wp : passband center frequency  ( 0.0 - 1.0 )
%  ce1 : Chebyshev error in 1st stopband
%  ce2 : Chebyshev error in 2nd stopband
%
%  itype : type of impluse response sequence
%          ( symmetric --- 1, anti-symmetric --- 2 )
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
% Original program author: Ivan W.Selesnick, Rice University
%
% <example>
     M  =  34;
     L  =  6;
     wp  = 0.4;
     ce1 =  0.1;
     ce2 =  0.1;
     itype = 1; 
% ----------------------------------------------------------------------
PF = 0;                              % flag : Plot Figures
if (L >= M/2) | (rem(L,2) == 1)
   fprintf(1,'\n    L must be even and be less than M/2 \n')
   return
end
if (itype ~= 1) & (itype ~= 2)
   fprintf(1,'\n    itype must be 1 or 2 \n')
   return
end
if rem(M,2) == 0
   if itype == 1
      ftype = 1;
   else
      ftype = 3;
      M = M - 2;
   end
else
   if itype == 1
      ftype = 2;
      M = M - 1;
   else
      ftype = 4;
      M = M - 1;
   end
end
N = M + 1;                              % filter length
g = 2^ceil(log2(20*N));                 % number of grid points
wp  = wp * pi;
SN = 1e-9;                              % SN : SMALL NUMBER
q  = (N-2*L+1)/2;                         % number of frequency points = q+1
% ---------- initialize some constants ------------------------------
w  = [0:g]'*pi/g;                    % frequency axis
ws1 = wp*0.8;
a = 5; b = 1;
ws2 = (a*wp+b*pi)/(a+b);
d = ws1/(pi-ws2);
q1 = round(q/(1+1/d));       % q1 : number of ref. set freq. in 1st stopband
if q1 == 0
   q1 = 1;
elseif q1 == q
   q1 = q - 1;
end
q2 = q - q1;                 % q2 : number of ref. set freq. in 2nd stopband
if q1 == 1;
   rs1 = ws1;
else
   if (ftype == 1) | (ftype == 2)
      rs1 = [0:q1-1]'*(ws1/(q1-1));
   elseif (ftype == 3) | (ftype == 4)
      rs1 = [0:q1]'*(ws1/q1);
      rs1(1) = [];
   end
end
if q2 == 1
   rs2 = ws2;
else
   if (ftype == 1) | (ftype == 4)
      rs2 = [0:q2-1]'*(pi-ws2)/(q2-1)+ws2;
   elseif (ftype == 2) | (ftype == 3)
      rs2 = [0:q2]'*(pi-ws2)/q2+ws2;
      rs2(length(rs2)) = [];
   end
end
rs = [rs1; rs2];
Y1 = [-ce1*(1-(-1).^(q1:-1:1))/2 + ce1*((-1).^(q1:-1:1)+1)/2]';
Y2 = [-ce2*(1-(-1).^(1:q2))/2 + ce2*((-1).^(1:q2)+1)/2]';
Y  = [Y1; Y2];
Z  = zeros(2*(g-q)+1,1);
n  = 0:q-1;
if ftype == 1
   B1 = 1;
   B1r = 1;
   B2 = 1;
   B2r = 1;
elseif ftype == 2
   B1 = cos(w/2);
   B1r = cos(rs/2);
   B2 = 0;
   B2r = 0;
   for i = 0:L-1
      B2 = B2 + taylor1(i,wp) * (((1-cos(w))/2)-((1-cos(wp))/2)).^i;
      B2r = B2r + taylor1(i,wp) * (((1-cos(rs))/2)-((1-cos(wp))/2)).^i;
   end
elseif ftype == 3
   B1 = sin(w);
   B1r = sin(rs);
   B21 = 0;
   B2r1 = 0;
   for i = 0:L/2
      B21 = B21 + taylor1(i,wp) * (((1-cos(w))/2)-((1-cos(wp))/2)).^i;
      B2r1 = B2r1 + taylor1(i,wp) * (((1-cos(rs))/2)-((1-cos(wp))/2)).^i;
   end
   B22 = 0;
   B2r2 = 0;
   for i = 0:L/2
      B22 = B22 + taylor2(i,wp) * (((1-cos(w))/2)-((1-cos(wp))/2)).^i;
      B2r2 = B2r2 + taylor2(i,wp) * (((1-cos(rs))/2)-((1-cos(wp))/2)).^i;
   end
   B2 = (B21 .* B22)/2;
   B2r = (B2r1 .* B2r2)/2;
else
   B1 = sin(w/2);
   B1r = sin(rs/2);
   B2 = 0;
   B2r = 0;
   for i = 0:L-1
      B2 = B2 + taylor2(i,wp) * (((1-cos(w))/2)-((1-cos(wp))/2)).^i;
      B2r = B2r + taylor2(i,wp) * (((1-cos(rs))/2)-((1-cos(wp))/2)).^i;
   end
end
A1  = (-1)^L * ((((1-cos(w))/2)-((1-cos(wp))/2)).^L);
A1r = (-1)^L * ((((1-cos(rs))/2)-((1-cos(wp))/2)).^L);
it = 0;
while 1 & (it < 15)
   % ---------- calculate interpolation values ----------
   a2 = cos(rs*n)\[((Y./B1r)-B2r)./A1r];
   A2 = real(fft([a2(1);a2(2:q)/2;Z;a2(q:-1:2)/2])); A2 = A2(1:g+1);
   A = B1.*(B2 + A1.*A2);
   % ---------- plot figures if PF == 1 ----------
   if PF
      figure(2)
      subplot(2,1,1), plot(w/pi,A), axis([0 1 -0.5 1.05]),
         hold on, plot(rs/pi,Y,'o'), hold off,
         xlabel('Normalized Angular Frequency'), ylabel('Amplitude')
         for i = 1:length(A)
            if A(i) == 0
               A(i) = SN;
            end
         end
      subplot(2,1,2), plot(w/pi,20*log10(abs(A))), axis([0 1 -120 5]),
         hold on, plot(rs/pi,20*log10(abs(Y)),'o'), hold off,
         xlabel('Normalized Angular Frequency'), ylabel('Amplitude [dB]'),
         for i = 1:length(A)
            if A(i) == SN
               A(i) = 0;
            end
         end
      pause(1)
   end
   % ---------- determine new reference set ----------
   ri = sort([local_max(A); local_max(-A)]);
   lri = length(ri);
   if (ftype == 2) | (ftype == 3)
      ri(lri) = [];
      lri = lri - 1;
   end
   if (ftype == 3) | (ftype == 4)
      ri(1) = [];
      lri = lri - 1;
   end
   if lri ~= q+1
      if abs(A(ri(lri))-A(ri(lri-1))) < abs(A(ri(1))-A(ri(2)))
         ri(lri) = [];
      else
         ri(1) = [];
      end
   end
   rs = (ri-1)*pi/g;
   [temp, k] = min(abs(rs - wp)); rs(k) = [];
   q1 = sum(rs < wp);
   q2 = sum(rs > wp);
   Y1 = [-ce1*(1-(-1).^(q1:-1:1))/2 + ce1*((-1).^(q1:-1:1)+1)/2]';
   Y2 = [-ce2*(1-(-1).^(1:q2))/2 + ce2*((-1).^(1:q2)+1)/2]';
   Y = [Y1; Y2];
   A1r = (-1)^L * ((cos(wp)-cos(rs))/2).^L;
   if ftype == 1
      B1r = 1;
      B2r = 1;
   elseif ftype == 2
      B1r = cos(rs/2);
      B2r = 0;
      for i = 0:L-1
         B2r = B2r + taylor1(i,wp) * (((1-cos(rs))/2)-((1-cos(wp))/2)).^i;
      end
   elseif ftype == 3
      B1r = sin(rs);
      B2r=1;
      B2r1 = 0;
      for i = 0:L/2
         B2r1 = B2r1 + taylor1(i,wp) * (((1-cos(rs))/2)-((1-cos(wp))/2)).^i;
      end
      B2r2 = 0;
      for i = 0:L/2
         B2r2 = B2r2 + taylor2(i,wp) * (((1-cos(rs))/2)-((1-cos(wp))/2)).^i;
      end
      B2r = (B2r1 .* B2r2)/2;
   else
      B1r = sin(rs/2);
      B2r = 0;
      for i = 0:L-1
         B2r = B2r + taylor2(i,wp) * (((1-cos(rs))/2)-((1-cos(wp))/2)).^i;
      end
   end
   Ar = B1r.*(B2r + (cos(rs*n)*a2) .* A1r);
   % ---------- calculate maximum constraint violation ----------
   for i = 1:length(Ar)
      if i <= length(Y1)
         Err(i) = abs(Ar(i)) - ce1;
      else
         Err(i) = abs(Ar(i)) - ce2;
      end
   end
   Err = max(Err);
   if PF
      fprintf(1,'    Err = %18.15f\n',Err);
   end
   if Err < SN
      break
   end
   it = it + 1;
end
% ---------- calculate filter coefficients and roots--------------------
ha2 = [a2(q:-1:2)/2; a2(1); a2(2:q)/2];
h = ha2;
for k = 1:L
   h = conv(h,[1 -2*cos(wp) 1]')/4;
end
if ftype == 1
   h((N+1)/2) = h((N+1)/2) + 1;
elseif ftype ~= 3
   hb2 = zeros(N,1); 
   for j = 0:L-1;
      x = [-0.25 cos(wp)/2 -0.25];
      k = 1;
      for i = 1:j
         k = conv(x,k);
      end
      if ftype == 2
         c = taylor1(j,wp) * k';
      else
         c = taylor2(j,wp) * k';
      end
      hb2 = hb2 + [zeros((N-((j*2)+1))/2,1); c; zeros((N-((j*2)+1))/2,1)];
   end
   h = h + hb2;
else
   hc1 = zeros(L+1,1);
   for j = 0:L/2;
      x = [-0.25 cos(wp)/2 -0.25];
      k = 1;
      for i = 1:j
         k = conv(x,k);
      end
      c1 = taylor1(j,wp) * k';
      hc1 = hc1 + [zeros(L/2-j,1); c1; zeros(L/2-j,1)];
   end
   hc2 = zeros(L+1,1);
   for j = 0:L/2;
      x = [-0.25 cos(wp)/2 -0.25];
      k = 1;
      for i = 1:j
         k = conv(x,k);
      end
      c2 = taylor2(j,wp) * k';
      hc2 = hc2 + [zeros(L/2-j,1); c2; zeros(L/2-j,1)];
   end
   c = conv(hc1,hc2)/2;
   hb2 = [zeros((N-length(c))/2,1); c; zeros((N-length(c))/2,1)];
   h = h + hb2;
end
% ---------- filter coefficients --------------------------------------------------------
if ftype == 2
   h = conv(h,[0.5 0.5]);
elseif ftype == 3
   h = conv(h,[0.5 0 -0.5]);
elseif ftype == 4
   h = conv(h,[0.5 -0.5]);
end
% ---------- calculate band edges -------------------------------------------------------
k = 1;
while 1
   j = A(round(g*wp/pi)-k) - ce1;
   if j < 0
      ew1 = round(g*wp/pi) - k;
      break
   end
   k = k + 1;
end
edge1 = w(ew1)/pi;
k = 1;
while 1
   j = A(round(g*wp/pi)+k) - ce2;
   if j < 0
      ew2 = round(g*wp/pi) + k;
      break
   end
   k = k + 1;
end
edge2 = w(ew2)/pi;
fprintf(1,'    1st Stopband edge         =  %8.5f\n',edge1)
fprintf(1,'    2nd Stopband edge         =  %8.5f\n',edge2)
c = roots(h);                           % calculate roots
figure(1) 
subplot(2,2,1), plot(w/pi,A), axis([0 1 -0.2 1.05]),
   hold on, plot(rs/pi,Y,'o'), hold off,
   xlabel('Normalized Angular Frequency'), ylabel('Amplitude')
   for i = 1:length(A)
      if A(i) == 0;
          A(i) = SN;
       end
   end
subplot(2,2,3), plot(w/pi,20*log10(abs(A))), axis([0 1 -100 5]),
   hold on, plot(rs/pi,20*log10(abs(Y)),'o'), hold off,
   xlabel('Normalized Angular Frequency'), ylabel('Amplitude [dB]'), 
subplot(2,2,2), impz(h,1)
   xlabel('Impulse Response')
subplot(2,2,4), zplane(c)