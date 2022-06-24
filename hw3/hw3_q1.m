clear all;close all; clc;

N = 100;

%  -------------(a 小題)--------------
% 資料1, 2, 3, 4
a_m1 = [-10;-10]; % 平均
a_m2 = [-10;10]; 
a_m3 = [10;-10]; 
a_m4 = [10;10]; 
a_S = [0.2 0; 0 0.2]; % 共變異數

a_data1 = mvnrnd(a_m1, a_S, N); % 產生高斯分佈資料
a_data2 = mvnrnd(a_m2, a_S, N); 
a_data3 = mvnrnd(a_m3, a_S, N); 
a_data4 = mvnrnd(a_m4, a_S, N); 

figure(1)
% 顯示資料1
plot(a_data1(:,1),a_data1(:, 2),'r+');
% 顯示資料2
hold on;
plot(a_data2(:,1),a_data2(:, 2),'g+');
% 顯示資料3
hold on;
plot(a_data3(:,1),a_data3(:, 2),'b+');
% 顯示資料4
hold on;
plot(a_data4(:,1),a_data4(:, 2),'black+');

% 合併資料
a_data = [a_data1 a_data2 a_data3 a_data4];
a_mean = [a_m1 a_m2 a_m3 a_m4];
% 計算 Sw, Sb, Sm
a_Sw = CalSw(a_data, a_mean);
a_Sb = CalSb(a_mean);
a_Sm = a_Sw + a_Sb;
% 計算 J3
a_J3 = trace(inv(a_Sw)*a_Sm);

%  -------------(b 小題)-------------- 
% 資料1, 2, 3, 4
b_m1 = [-1;-1]; % 平均
b_m2 = [-1;1]; 
b_m3 = [1;-1]; 
b_m4 = [1;1]; 
b_S = [0.2 0; 0 0.2]; % 共變異數

b_data1 = mvnrnd(b_m1, b_S, N); % 產生高斯分佈資料
b_data2 = mvnrnd(b_m2, b_S, N); 
b_data3 = mvnrnd(b_m3, b_S, N); 
b_data4 = mvnrnd(b_m4, b_S, N); 

figure(2)
% 顯示資料1
plot(b_data1(:,1),b_data1(:, 2),'r+');
% 顯示資料2
hold on;
plot(b_data2(:,1),b_data2(:, 2),'g+');
% 顯示資料3
hold on;
plot(b_data3(:,1),b_data3(:, 2),'b+');
% 顯示資料4
hold on;
plot(b_data4(:,1),b_data4(:, 2),'black+');

% 合併資料
b_data = [b_data1 b_data2 b_data3 b_data4];
b_mean = [b_m1 b_m2 b_m3 b_m4];
% 計算 Sw, Sb, Sm
b_Sw = CalSw(b_data, b_mean);
b_Sb = CalSb(b_mean);
b_Sm = b_Sw + b_Sb;
% 計算 J3
b_J3 = trace(inv(b_Sw)*b_Sm);

%  -------------(c 小題)-------------- 
% 資料1, 2, 3, 4
c_S = [3 0; 0 3]; % 共變異數

c_data1 = mvnrnd(a_m1, c_S, N); % 產生高斯分佈資料
c_data2 = mvnrnd(a_m2, c_S, N); 
c_data3 = mvnrnd(a_m3, c_S, N); 
c_data4 = mvnrnd(a_m4, c_S, N); 

figure(3)
% 顯示資料1
plot(c_data1(:,1),c_data1(:, 2),'r+');
% 顯示資料2
hold on;
plot(c_data2(:,1),c_data2(:, 2),'g+');
% 顯示資料3
hold on;
plot(c_data3(:,1),c_data3(:, 2),'b+');
% 顯示資料4
hold on;
plot(c_data4(:,1),c_data4(:, 2),'black+');

% 合併資料
c_data = [c_data1 c_data2 c_data3 c_data4];
c_mean = [a_m1 a_m2 a_m3 a_m4];
% 計算 Sw, Sb, Sm
c_Sw = CalSw(c_data, c_mean);
c_Sb = CalSb(c_mean);
c_Sm = c_Sw + c_Sb;
% 計算 J3
c_J3 = trace(inv(c_Sw)*c_Sm);


% 計算 Sw = (1/N) * sum of 1:k class Si
% X - 資料 Nx8
% M - 平均 2x4
% res - Sw
function [res] = CalSw(X, M)
[N, d] = size(X);
res = zeros(2,2);
for i = 1:d/2
    for j = 1:N
        res = res + (X(j, 2*i-1:2*i)' - M(:, i)) * ( X(j, 2*i-1:2*i)' - M(:, i))';
    end
end
res = res / (N*d/2);
end

% 計算 Sb = sigma Pi * (mi - M) * (mi - M)'
% M - 平均 2x4
% res - Sb
function [res] = CalSb(M)
globalM = zeros(2,1);
res = 0;
for i = 1:4
    globalM = globalM + (1/4) * M(:, i);
end
for i = 1:4
    res = res + (1/4) * (M(:, i) - globalM)*(M(:, i) - globalM)';
end
end