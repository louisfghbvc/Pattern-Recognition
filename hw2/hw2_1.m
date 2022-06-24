clear all;close all; clc;

N = 1000;
% 資料1
P = 0.25;
% 產生白努利分布資料
data = binornd(1, P, [1 N]);
% ML0估算
est1 = Bernoulli_MLEstimatior(P, data);

% 資料2
P = 0.5;
% 產生白努利分布資料
data2 = binornd(1, P, [1 N]);
% ML估算結果
est2 = Bernoulli_MLEstimatior(P, data2);

% 白努利 maximum likelihood公式 
% P - 機率
% X - 資料
function [est] = Bernoulli_MLEstimatior(P,X)
[tmp,N] = size(X);
syms P
% 估算式子
loss_p = log(P)*sum(X) + log(1-P)*(N-sum(X));

% 對P微分之後 計算等於0的值
% est = sum(X) / N;
est = double( solve(diff(loss_p, P) == 0) );
end 