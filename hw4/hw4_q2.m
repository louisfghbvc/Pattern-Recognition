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

% run k means
[a_est_means1, a_mnI1] = k_means(a_data, 2, 1000);
[a_est_means2, a_mnI2] = k_means(a_data, 3, 1000);
[a_est_means3, a_mnI3] = k_means(a_data, 4, 1000);

% 畫圖
figure(1)
subplot(2,2,1);
DrawData(a_data);

subplot(2,2,2);
gscatter(a_data(:, 1), a_data(:, 2), a_mnI1);
subplot(2,2,3);
gscatter(a_data(:, 1), a_data(:, 2), a_mnI2);
subplot(2,2,4);
gscatter(a_data(:, 1), a_data(:, 2), a_mnI3);

% ----------------------- (b) 小題 ------------------------
% 平均
b_m1 = [1;1];
b_m2 = [3.5;3.5];
b_m3 = [6;1];

% 產生 dataset
b_data = GenerateData(b_m1, b_m2, b_m3, s1, s2, s3, N);

% run k means
[b_est_means1, b_mnI1] = k_means(b_data, 2, 1000);
[b_est_means2, b_mnI2] = k_means(b_data, 3, 1000);
[b_est_means3, b_mnI3] = k_means(b_data, 4, 1000);

% 畫圖
figure(2)
subplot(2,2,1);
DrawData(b_data);

subplot(2,2,2);
gscatter(b_data(:, 1), b_data(:, 2), b_mnI1);
subplot(2,2,3);
gscatter(b_data(:, 1), b_data(:, 2), b_mnI2);
subplot(2,2,4);
gscatter(b_data(:, 1), b_data(:, 2), b_mnI3);


% --------------------- (c) 小題 --------------------------
% 平均
c_m1 = [1;1];
c_m2 = [2;2];
c_m3 = [3;1];

% 產生 dataset
c_data = GenerateData(c_m1, c_m2, c_m3, s1, s2, s3, N);

% run k means
[c_est_means1, c_mnI1] = k_means(c_data, 2, 1000);
[c_est_means2, c_mnI2] = k_means(c_data, 3, 1000);
[c_est_means3, c_mnI3] = k_means(c_data, 4, 1000);

% 畫圖
figure(3)
subplot(2,2,1);
DrawData(c_data);

subplot(2,2,2);
gscatter(c_data(:, 1), c_data(:, 2), c_mnI1);
subplot(2,2,3);
gscatter(b_data(:, 1), c_data(:, 2), c_mnI2);
subplot(2,2,4);
gscatter(c_data(:, 1), c_data(:, 2), c_mnI3);

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

% K means algorithms
% X - data N x 2
% k - numbers of center
% max_iter - 最大疊代數
function [means, mnI] = k_means(X, k, max_iter)
    [N, l] = size(X);
    
    % dim x numbers of center
    means = rand(l, k);

    e = 1;
    iter = 0;
    while e ~= 0 && iter < max_iter

        means_old = means;
        dist_all = [];
        
        for j = 1 : k
            extand_mean = ones(N, 1) * means(:, j)';
            % sum the distance, no need sqrt, we just compare.
            dist = sum( ((extand_mean - X).^2)');
            
            % append on matrix
            dist_all = [dist_all; dist];
        end
        
        [~, mnI] = min(dist_all);

        for j = 1 : k
            % calculate same label j
            num = sum(mnI == j);
            if( num ~= 0)
                % convert to same dim with data                 
                extand_ = (mnI==j)'* ones(1,l);
                means(:, j) = sum( X.* extand_) / num;
            end
        end

        e = sum( sum(abs(means - means_old)) );
        iter = iter + 1;
    end
end