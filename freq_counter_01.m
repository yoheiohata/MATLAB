clear;close all;clc
fs = 1000;
t = (0:1/fs:100)';
x = sin(30*pi*t) + 0.67*sin(40*pi*t) + 0.33*sin(50*pi*t) + 1;

dcblker = dsp.DCBlocker('Algorithm','FIR','Length',100);


hsa = dsp.SpectrumAnalyzer('SampleRate',fs, ...
    'PowerUnits','dBW','FrequencySpan','Start and stop frequencies',...
    'StartFrequency',-30,'StopFrequency',30,'YLimits',[-200 20],...
    'Title','Signal Spectrum');


hsb = clone(hsa);
hsb.Title = 'Signal Spectrum After DC Blocker';

y = dcblker(x);

hsa(x);

hsb(y);
mydata1=xlsread('freq_analysis.xlsx','A4:B259');
x=mydata1(:,1);
y=mydata1(:,2);

plot(x,y,'k');
% M = readmatrix('freq_analysis.xlsx');
% function [y,tOut] = removeDC(x,tIn)
% % Remove the DC value of a signal by subtracting its mean
%    y = x - mean(x);
%    tOut = tIn;
% end
% function [y,tOut] = timealign(x,tIn,startTime)
% % Change the starting time of a signal
%    y = x;
%    t = tIn;
%    if ~isempty(t)
%        t = t - t(1) + startTime;
%    end
%    tOut = t;
% end

        
