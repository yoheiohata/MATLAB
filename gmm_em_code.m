clear;close all;clc

%平均
mu=[4; -3];
%分散 
sigma=cat(3,4,7);
%混合比
p = ones(1,2)/2;

%1次元混合ガウスモデル
GMM = gmdistribution(mu,sigma,p);

x=-10:0.1:10;
x = transpose(x);
y = pdf(GMM,x);


figure
plot(x,y, '-b');

mu1 = [1 2];
Sigma1 = [2 0; 0 0.5];
mu2 = [-1 -2];
Sigma2 = [1 0;0 1];
rng(1); % For reproducibility
X = [mvnrnd(mu1,Sigma1,1000);mvnrnd(mu2,Sigma2,1000)];

%混合ガウス モデルを近似します。成分数を 2 に指定します。

GMModel = fitgmdist(X,2);

%近似された混合ガウス モデルの等高線上にデータをプロットします。

figure
y = [zeros(1000,1);ones(1000,1)];
h = gscatter(X(:,1),X(:,2),y);
hold on
ezcontour(@(x1,x2)pdf(GMModel,[x1 x2]),get(gca,{'XLim','YLim'}))
title('{\bf Scatter Plot and Fitted Gaussian Mixture Contours}')
legend(h,'Model 0','Model1')
hold off


