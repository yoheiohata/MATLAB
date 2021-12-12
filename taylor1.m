function Ta1 = taylor1(N,wo)
% Ta1 = taylor1(N,wo)
% This program returns the coefficients of Nth order of 
% Taylor's series of (1-x)^-1/2 at x=xo.
% 0.0 =< xo < 1.0 that is, 0.0 =< wo < pi
%
% Author: Jun Hashimoto, Tokyo Metropolitan University
% ----------------------------------------------------------
i = N - 1;
if i == -1
   T=1;
elseif i==0
   T=1/2;
else
   a = 1;
   b = 2;
   for k = 1:i
      a = a * (2*k+1);
      b = b * (2*k+2);
   end
   T = a / b;
end
xo = (1-cos(wo))/2;
Ta1 = ((1-xo)^(-1/2)) * T * ((1/(1-xo))^N);
