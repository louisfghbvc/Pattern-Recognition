clear all;close all; clc;

N = 500;

% ---------------------------- (a) 小題 ------------------------------
% 平均
a_m1 = [1;1];
a_m2 = [5;5];
a_m3 = [9;1];

% 共變異數
s1 = [1 0.4; 0.4 1];
s2 = [1 -0.6; -0.6 1];
s3 = [1 0; 0 1];

% 產生 dataset
a_data = GenerateData(a_m1, a_m2, a_m3, s1, s2, s3, N);

k = 3;
q = 2;
% 實驗 1, k=3, q=2
[a_mean1, a_U1, ~] = fuzzy_k_means(a_data, k, q);
q = 3;
% 實驗 1, k=3, q=3
[a_mean2, a_U2, ~] = fuzzy_k_means(a_data, k, q);

% 取最大
[~, a_id1] = max(a_U1);
[~, a_id2] = max(a_U2);

% 畫圖
figure(1)
subplot(3,1,1);
DrawData(a_data);
subplot(3,1,2);
gscatter(a_data(:, 1), a_data(:, 2), a_id1);
subplot(3,1,3);
gscatter(a_data(:, 1), a_data(:, 2), a_id2);

% ----------------------- (b) 小題 ------------------------
% 平均
b_m1 = [1;1];
b_m2 = [3.5;3.5];
b_m3 = [6;1];

% 產生 dataset
b_data = GenerateData(b_m1, b_m2, b_m3, s1, s2, s3, N);

k = 3;
q = 2;
% 實驗 1, k=3, q=2
[b_mean1, b_U1, ~] = fuzzy_k_means(b_data, k, q);
q = 3;
% 實驗 1, k=3, q=3
[b_mean2, b_U2, ~] = fuzzy_k_means(b_data, k, q);

% 取最大
[~, b_id1] = max(b_U1);
[~, b_id2] = max(b_U2);

% 畫圖
figure(2)
subplot(3,1,1);
DrawData(b_data);
subplot(3,1,2);
gscatter(b_data(:, 1), b_data(:, 2), b_id1);
subplot(3,1,3);
gscatter(b_data(:, 1), b_data(:, 2), b_id2);

% --------------------- (c) 小題 --------------------------
% 平均
c_m1 = [1;1];
c_m2 = [2;2];
c_m3 = [3;1];

% 產生 dataset
c_data = GenerateData(c_m1, c_m2, c_m3, s1, s2, s3, N);

k = 3;
q = 2;
% 實驗 1, k=3, q=2
[c_mean1, c_U1, ~] = fuzzy_k_means(c_data, k, q);
q = 3;
% 實驗 1, k=3, q=3
[c_mean2, c_U2, ~] = fuzzy_k_means(c_data, k, q);

% 取最大
[~, c_id1] = max(c_U1);
[~, c_id2] = max(c_U2);

% 畫圖
figure(3)
subplot(3,1,1);
DrawData(c_data);
subplot(3,1,2);
gscatter(c_data(:, 1), c_data(:, 2), c_id1);
subplot(3,1,3);
gscatter(c_data(:, 1), c_data(:, 2), c_id2);



% fuzzy_k_means
% X - data N x 2
% m - numbers of center
% q - The criterion allows each pattern to belong to multiple clusters
function [theta, U, obj_fun] = fuzzy_k_means(X, m, q)
    options(1) = q;
    [theta, U, obj_fun] = fcm(X, m, options);
    theta = theta';
end

% Generate Dataset
% m - mean
% s - varience
% N - size
function [data] = GenerateData(m1, m2, m3, s1, s2, s3, N)
    % 產生 dataset
    data = zeros(N, 2);
    for i = 1:N
       if mod(i, 4) == 1 || mod(i, 4) == 2
           data(i,:) = mvnrnd(m2, s2, 1);
       elseif mod(i, 4) == 3
           data(i,:) = mvnrnd(m1, s1, 1);
       else
           data(i,:) = mvnrnd(m3, s3, 1);
       end
    end
end

% Draw Data
% data - N x 2, input data
function DrawData(data)
    [N, ~] = size(data);
    for i = 1 : N
       if mod(i, 4) == 1 || mod(i, 4) == 2
            plot(data(i, 1), data(i, 2),'r+');
            hold on;
       elseif mod(i, 4) == 3
            plot(data(i, 1), data(i, 2),'g+');
            hold on;
       else
            plot(data(i, 1), data(i, 2),'b+');
            hold on;
       end
    end
end