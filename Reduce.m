function [ t,f,ind ] = Reduce( t,f,N )
%REDUCE Reduce the resolution of a signal for plotting
%   [ t,f ] = Reduce( t,f,N )
%   t = time, f = function
%   Samples at default N=16x the fundamental frequency as determined by 
%       zero crossings

if nargin < 3
    N = 16;
end

if N < 4
	warning('REDUCE:MIN','Reducing to < 4x the fundamental frequency is not recommended');
end

if size(f,1) == 1 || size(f,2) == 1
	f = f(:)';
end

delta_min = Inf;

for ii = 1:size(f,1)
	f_c = f(ii,:);
	f_c = f_c - mean(f_c);

	if length(t) ~= length(f_c);
		error('REDUCE:LENGTH','Lengths must match');
	end

	f1 = f_c(2:end);
	f2 = f_c(1:end-1);

	Ncrossing = length(find(f2 >= 0 & f1 <= 0));

	delta = round(length(f)/Ncrossing/N);
	if(delta < delta_min)
		delta_min = delta;
	end

end

ind = 1:delta_min:length(t);

t = t(ind);
f = f(:,ind);

end

