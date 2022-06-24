clear all;close all; clc;

%  -------------(a 小題)--------------
l = 2;
N = 1000;
w = [1;1];
w0 = 0;
a = 10;
e = 1;
sed = 0;
X = generate_hyper(w, w0, a, e, N, sed);

%  -------------(b 小題)--------------
Y = PCA(X, 2);
hold on;
figure(1);
plot(Y(1, :), Y(2, :), '.g');

% 依照 題目要求 產生hyper plane
% w - 權重
% w0 - 常數
% a - uniform distribute points[-a, a]
% e - uniform distribute noise [-a, a]
% N - 點的數目
% sed - 亂數種子
% X - 回傳結果
function [X] = generate_hyper(w, w0, a, e, N, sed)
    randn('seed', sed);
    [l, tmp] = size(w);
    t = (rand(l-1, N) - 0.5) * 2 * a;
    t_last = -(w(1:l-1) / w(1))' * t + 2 * e * (rand(1,N) - 0.5) - (w0 / w(1));
    X = [t; t_last];
    if(l == 2)
        figure(1);
        plot(X(1, :), X(2, :), '.b');
    end
    figure(1);
    axis equal;
end

% 計算 PCA
% X - 資料
% m - 想要的維度 
% Y - 轉換X的資料
function [Y] = PCA(X, m)
[l,N]=size(X);

% 計算 mean 並轉置, 減去平均
mean_vec = mean(X')';
X_zero= X - mean_vec * ones(1,N);

% 計算 eigen值
R = cov(X_zero');
[V, D] = eig(R);

eigenval = diag(D);
[tmp,ind] = sort(eigenval, 1, 'descend');
eigenvec = V(:,ind);

eigenvec = eigenvec(:,1:m);

% 計算轉置矩陣
A = eigenvec(:,1:m)';

% 轉換 Y
Y = A * X;
end