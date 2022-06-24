clear all;close all; clc;
% randn('seed', 0);

N = 100;
%  -------------(a 小題)--------------
% 資料1, 2.
% 平均
a_m1 = [2;4]; 
a_m2 = [2.5;10]; 
a_S = [1 0; 0 1];

% 產生高斯分佈資料
a_data1 = mvnrnd(a_m1, a_S, N); 
a_data2 = mvnrnd(a_m2, a_S, N); 

figure(1)
% 顯示資料1
plot(a_data1(:,1),a_data1(:, 2),'r+');
% 顯示資料2
hold on;
plot(a_data2(:,1),a_data2(:, 2),'g+');

% 計算 FDR
a_FDR = FDR(a_data1, a_data2);

%  -------------(b 小題)--------------
% 資料1, 2
% 平均
b_m1 = [2;4];
b_m2 = [2.5;10]; 
b_S = [0.25 0; 0 0.25];

% 產生高斯分佈資料
b_data1 = mvnrnd(b_m1, b_S, N); 
b_data2 = mvnrnd(b_m2, b_S, N); 

figure(2)
% 顯示資料1
plot(b_data1(:,1),b_data1(:, 2),'r+');
% 顯示資料2
hold on;
plot(b_data2(:,1),b_data2(:, 2),'g+');

% 計算 FDR
b_FDR = FDR(b_data1, b_data2);

% 計算 FDR 
% a - 資料 a
% b - 資料 b
% res - 回傳結果 2個feature 的FDR(1x2)
function [res] = FDR(a, b)
a_m = mean(a);
b_m = mean(b);
a_var = var(a);
b_var = var(b);

res1 = (a_m(:, 1) - b_m(:, 1)) ^2 / (a_var(:, 1) + b_var(:, 1));
res2 = (a_m(:, 2) - b_m(:, 2)) ^2 / (a_var(:, 2) + b_var(:, 2));
res = [res1, res2];
end