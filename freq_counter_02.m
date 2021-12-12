clear;close all;clc

mydata1=xlsread('freq_analysis.xlsx','A4:B259');
x=mydata1(:,1);
y_each=mydata1(:,2);
y_mean=mean(y_each);   %　DC成分(オフセット)を平均値として計算

y1=y_each-y_mean;
y=y1(:,1);

zeroCross = y(1:end-1).*y(2:end) < 0;  % 符号が変化しているところを読み取る

mydata=[zeroCross x(1:end-1)];

[x_1 y_1]=size(mydata);
n=0;
for m=1:x_1
    if  mydata(m,1)==1
        n=n+1;
        z(n,:)=mydata(m,2);   %2列目を配列に入れて代入
        
    end
end

[z_1]=length(z);
p=0;
for q=1:(z_1)-1
    p=p+1;
    z_2(p,:)=z(q+1)-z(q);
end

T_half=mean(z_2);
freqency=1/(2*T_half);



figure
plot(x,y,'k');




