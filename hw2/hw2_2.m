clear all;close all; clc;

N = 1000;
% 資料1
m1 = [1;1]; % 平均
s1 = [5 3; 3 4]; % 共變異數
data1 = mvnrnd(m1, s1, N); % 產生高斯分佈資料

% 資料2
m2 = [10;5]; % 平均
s2 = [7 4; 4 5]; % 共變異數
data2 = mvnrnd(m2, s2, N); % 產生高斯分佈資料

% 顯示資料1
figure(1);
plot(data1(:,1),data1(:, 2),'b+');

% 顯示資料2
figure(2);
plot(data2(:,1),data2(:, 2),'g+');

% 資料1 ML估算結果
[est_mean1, est_var1] = Normal_MLEstimatior(data1);
% 資料2 ML估算結果
[est_mean2, est_var2] = Normal_MLEstimatior(data2);

% 常態分佈 maximum likelihood 
% X - 資料
function [est_mean, est_var] = Normal_MLEstimatior(X)
[N,tmp] = size(X);
% 常態分佈 平均估算
est_mean = sum(X) / N;

% 常態分佈 變異數估算
est_var = zeros(2);
for i = 1:N
    est_var = est_var + (X(i,:) - est_mean)' * (X(i,:) - est_mean);
end
est_var = est_var / N;
end 