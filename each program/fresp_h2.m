function [des,wt] = fresp_h3(order, ff, grid, wtx, aa, varargin)
%	d‚ÝŠÖ”•ÏX
%	Makoto Nakatani 2006/08/05
%   Modify N. Aikawa 2006/08/22
% fs=2000;
nbands = length(ff)/2;

  h = load('notch2.cof');
  h = h';
  des=grid;
  wt=grid;

for i=1:nbands
    k = find(grid >= ff(2*i-1) & grid <= ff(2*i));
    ww = pi.*grid(k)'*linspace(0,length(h)-1,length(h));
    pref(k)=abs(cos(ww)*h+j*sin(ww)*h);
    if (aa(2*i-1)==1)
        des(k)=aa(2*i-1)./pref(k);
    else
        des(k)=0;
    end;
    wt(k)=wtx(i)*pref(k); 
end
