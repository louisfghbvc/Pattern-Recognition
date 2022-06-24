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

% 畫圖
figure(1)
subplot(2,1,1);
% 顯示資料a
DrawData(a_data);

% Running em algorithm
[a_em, a_es, a_p] = ExpectationMaximization(a_data, 500, 0);

% Compare Result after EM
% 顯示EM結果 a
subplot(2,1,2);
str = ["r+", "g+", "b+"];
for i = 1:3
    tmp = mvnrnd(a_em(2*i-1 : 2*i, :), a_es(2*i-1 : 2*i, :), fix(N * a_p(i)));
    plot(tmp(:, 1), tmp(:, 2), str(i));
    hold on;
end

% ----------------------- (b) 小題 ------------------------
% 平均
b_m1 = [1;1];
b_m2 = [3.5;3.5];
b_m3 = [6;1];

% 產生 dataset
b_data = GenerateData(b_m1, b_m2, b_m3, s1, s2, s3, N);

% 畫圖
figure(2)
subplot(2,1,1);
% 顯示資料b
DrawData(b_data);

% Running em algorithm
[b_em, b_es, b_p] = ExpectationMaximization(b_data, 500, 0);

% Compare Result after EM
% 顯示EM結果 b
subplot(2,1,2);
for i = 1:3
    tmp = mvnrnd(b_em(2*i-1 : 2*i, :), b_es(2*i-1 : 2*i, :), fix(N * b_p(i)));
    plot(tmp(:, 1), tmp(:, 2), str(i));
    hold on;
end

% --------------------- (c) 小題 --------------------------
% 平均
c_m1 = [1;1];
c_m2 = [2;2];
c_m3 = [3;1];

% 產生 dataset
c_data = GenerateData(c_m1, c_m2, c_m3, s1, s2, s3, N);

% 畫圖
figure(3)
subplot(2,1,1);
% 顯示資料c
DrawData(c_data);

% Running em algorithm
[c_em, c_es, c_p] = ExpectationMaximization(c_data, 500, 0);

% Compare Result after EM
% 顯示EM結果 c
subplot(2,1,2);
for i = 1:3
    tmp = mvnrnd(c_em(2*i-1 : 2*i, :), c_es(2*i-1 : 2*i, :), fix(N * c_p(i)));
    plot(tmp(:, 1), tmp(:, 2), str(i));
    hold on;
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

% EM algorithm
% data - input data N x 2
% T - 最大疊代數
% seed - Random seed
function [em, es, p] = ExpectationMaximization(data, T, seed)
    randn('seed', seed);
    
    [N, ~] = size(data);
    % initial values 
    em = rand(6, 1); % mean
    
    % random a positive det
    A = randn(2);
    A = A'*A;
    B = randn(2);
    B = B'*B;
    C = randn(2);
    C = C'*C;
    es = [A; B; C;];
    
    p = [0.3; 0.3; 0.4];
    
    % error
    e = 1;
    t = 0;
    % each iteration, two condition
    while e ~= 0 && t < T
        next_p = zeros(3, 1);
        next_mean = zeros(6, 1);
        next_var = zeros(6, 2);
        
        % sum of data, each posterior j
        sigma_p = zeros(3, 1);
        
        for j = 1:N
            % E-step
            posterior = zeros(3, 1);
            posterior_sum = 0;
            for i = 1:3
                cur_es = es(2*i-1:2*i,:);
                cur_em = em(2*i-1:2*i,:);
                posterior(i,:) = sqrt(det(inv(cur_es))) * exp( (-0.5) * (data(j, :)' - cur_em)' * inv(cur_es) * (data(j, :)' - cur_em) ) * p(i, :);
                posterior_sum = posterior_sum + posterior(i,:);
            end
            for i = 1:3
                posterior(i,:) = posterior(i,:) / posterior_sum;
            end
            
            % M-step
            for i = 1:3
                % sum mean
                next_mean(2*i-1:2*i, :) = next_mean(2*i-1:2*i, :) + posterior(i, :) * data(j,:)';
                sigma_p(i, :) = sigma_p(i, :) + posterior(i, :);
                
                % sum varience
                cur_em = em(2*i-1:2*i,:);
                next_var(2*i-1:2*i, :) = next_var(2*i-1:2*i, :) + posterior(i, :) * (data(j,:)' - cur_em) * (data(j,:)' - cur_em)';
            end
        end
        
        % M-step
        % determine the next value
        for i = 1:3
            next_mean(2*i-1:2*i, :) = next_mean(2*i-1:2*i, :) / sigma_p(i, :);
            next_var(2*i-1:2*i, :) = next_var(2*i-1:2*i, :) /  sigma_p(i, :);
            next_p(i, :) = (1/N) * sigma_p(i, :);
        end
        
        % 與前一個的差異總和
        e = sum(sum(abs(next_mean - em))) + sum(sum(abs(next_var - es))) + sum(abs(next_p - p));
        
        em = next_mean;
        es = next_var;
        p = next_p;
        
        t = t+1;
    end
    
end