%function [h,del] = flchebp1(M,L,wp,ws,wt,itype)
%　遷移域の幅をきめて減衰量を決める
% 
% [h,del] = flchebp1(M,L,wp,ws,wt,itype)
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
%  M  : filter order　フィルタの次数
%  L  : degree of flatness at w=wp  ( L < M/2 )　多項式の次数
%  wp : passband center frequency  ( 0.0 - 1.0 )　通過域中心周波数
%  ws : distance from wp to stopband edge　通過域中心周波数wpから阻止域端までの距離
%        ( 0 < wp-ws < wp < wp+ws < 1 )
%  wt  : weight ratio　重みづけの割合(分数になってる↓)
%        ( weight in 2nd stopband / weight in 1st stopband )
%          二番目(右側)の阻止域の重み/一番目(左側)の阻止域の重み
%  itype : type of impluse response sequence
%          ( symmetric --- 1, anti-symmetric --- 2 )
%           インパルス応答列の種類（偶対称なら1、奇対称なら2）
% Author: Jun Hashimoto, Tokyo Metropolitan University
% Update: Naoyuki Aikawa, Tokyo Metropolitan University
% Original program author: Ivan W.Selesnick, Rice University
%
% 論文
%   相川　直幸，佐藤　正光
% 　通過域平たんかつ阻止域等リプル特性を有する直線位相FIRディジタルフィルタの設計
%   電子情報通信学会論文誌(A) , Vol. J83-A No. 6.，pp744-749.，(2000,6)
%
% <example>
     M   =  28;     %フィルタの次数M=28 
     L   =  8;      %フィルタの多項式の次数
     wp  = 0.595;   %通過域中心周波数wp=0.595
     ws  = 0.35;    %通過域中心周波数wpから阻止域端までの距離  ws=0.35
     wt  = 1;       %阻止域の割合は同じ
     itype = 1;     %インパルス応答は偶対称
% ----------------------------------------------------------------------
%
PF = 0;                              % PF : flag : Plot Figures
if (L >= M/2) | (rem(L,2) == 1)      %多項式次数Lがフィルタの次数の半分　または　多項式の次数が奇数　
   fprintf(1,'\n    L must be even and be less than M/2 \n')
   return
elseif (0 >= wp-ws) | (wp-ws >= wp) | (wp >= wp+ws) | (wp+ws >= 1)      %通過域と遷移域の幅を指定
   fprintf(1,'\n    need : 0 < wp-ws < wp < wp+ws < 1 \n')
   return
elseif (itype ~= 1) & (itype ~= 2)
   fprintf(1,'\n    itype must be 1 or 2 \n')                           %インパルス応答のタイプを指定
   return
end
if rem(M,2) == 0
   if itype == 1            %type1フィルタ
      ftype = 1;
   else
      ftype = 3;            %type3フィルタ
      M = M - 2;
   end
else
   if itype == 1            %type2フィルタ
      ftype = 2;
      M = M - 1;
   else
      ftype = 4;            %type4フィルタ
      M = M - 1;
   end
end
N = M + 1;                              % filter length
g = 2^ceil(log2(20*N));                 % 格子点の個数？ ceil関数...正の無限大方向に丸めるnumber of grid points
ws1 = (wp - ws) * pi;
ws2 = (wp + ws) * pi;
wp  = wp * pi;
SN = 1e-9;                              % SN : SMALL NUMBER
q  = (N-2*L+1)/2;                         % number of frequency points = q+1
% ---------- initialize some constants (初期条件)------------------------------
w = [0:g]'*pi/g;                       % frequency axis
d = ws1/(pi-ws2);
q1 = round((q+1)/(1+1/d));   % q1 : number of ref. set freq. in 1st stopband
if q1 == 0
   q1 = 1;
elseif q1 == q+1
   q1 = q;
end
q2 = q + 1 - q1;             % q2 : number of ref. set freq. in 2nd stopband
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
Z = zeros(2*(g+1-q)-1,1);
si = (-1).^[0:q]';
for i=1:length(rs1)           % weight
   si(i) = si(i) * wt;
end
s1n = length(rs1);
n = 0:q-1;
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
   x = [cos(rs*n), si./(B1r.*A1r)]\[B2r./A1r];
   a = -x(1:q);
   del = x(q+1);
   A2 = real(fft([a(1);a(2:q)/2;Z;a(q:-1:2)/2])); A2 = A2(1:g+1);
   A = B1.*(B2 + A1.*A2);
   Y  = si*del;
   % ---------- plot figures if PF == 1 ----------
   if PF
      figure(2)
      subplot(2,1,1), plot(w/pi,A), axis([0 1 -0.2 1.05]),
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
      pause(2)
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
   for i=1:s1n
      si(i) = si(i) / wt;
   end
   rs = (ri-1)*pi/g;
   [temp, k] = min(abs(rs - wp));
   rs(k) = [];
   Aws1 = 1 + (cos(ws1*n)*a).*((-1)^L * (sin(ws1/2-wp/2).*sin(ws1/2+wp/2)).^L);
   Aws2 = 1 + (cos(ws2*n)*a).*((-1)^L * (sin(ws2/2-wp/2).*sin(ws2/2+wp/2)).^L);
   if (Aws1 > Aws2) | any((rs > wp) & (rs < ws2))
      rs = sort([ws1; rs]);
   else
      rs = sort([ws2; rs]);
   end
   s1n = 0;
   for i = 1:length(rs)
      if rs(i) < wp
         s1n = s1n + 1;
      end
   end
   for i=1:s1n
      si(i) = si(i) * wt;
   end
   A1r = (-1)^L * (sin(rs/2-wp/2).*sin(rs/2+wp/2)).^L;
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
   Ar = B1r.*(B2r + (cos(rs*n)*a) .* A1r);
   % ---------- calculate maximum constraint violation ----------
   for i = 1:s1n
      Ar(i) = Ar(i) / wt;
   end
   Err = max([max(Ar)-abs(del), abs(del)-abs(min(Ar))]);
   if PF
      fprintf(1,'    Err = %18.15f\n',Err);
   end
   if Err < SN
      break
   end
   it = it + 1;
end
% ---------- calculate filter coefficients and roots----------
ha = [a(q:-1:2)/2; a(1); a(2:q)/2];
h = ha;
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

fprintf(1,'    1st Stopband attenuation  =  %8.5f (%8.5f[dB])\n',abs(del*wt),20*log10(del*wt))
fprintf(1,'    2nd Stopband attenuation  =  %8.5f (%8.5f[dB])\n',abs(del),20*log10(del))
c = roots(h);                           % calculate roots
figure(1) 
subplot(2,2,1), plot(w/pi,A), axis([0 1 -0.2 1.05]),
   hold on, plot(rs/pi,Y,'o'), hold off,
   xlabel('Normalized Angular Frequency'), ylabel('Amplitude')
   for i = 1:length(A)
      if A(i) == 0
         A(i) = SN;
      end
   end
subplot(2,2,3), plot(w/pi,20*log10(abs(A))), axis([0 1 -110 5]),
   hold on, plot(rs/pi,20*log10(abs(Y)),'o'), hold off,
   xlabel('Normalized Angular Frequency'), ylabel('Amplitude [dB]'), 
subplot(2,2,2), impz(h,1),
   xlabel('Impulse Response')
subplot(2,2,4), zplane(c)
fid=fopen('foddbp','w');
for i=1:size(A)
   fprintf(fid,'%12.8f , %12.8f \r',w(i)/pi,20*log10(abs(A(i))));
end;
fclose(fid);
fid=fopen('roddbp','w');
[mm nn]=size(c);
for i=1:mm
   fprintf(fid,'%12.8f , %12.8f \r',real(c(i)),imag(c(i)));
end;
fclose(fid);
