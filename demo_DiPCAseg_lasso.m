%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script demonstrates the time series segmentation algoirthm 
% introduced in the following paper:
%   S. Lu and S. Huang, "Segmentation of multivariate industrial time 
%   series data based on dynamic latent variable predictability," 
%   IEEE Access, vol. 8, pp. 112092¨C112103, 2020.
%
% Written by : Shaowen Lu, Shuyu Huang @ Northeastern University, 2019.
%

global workDir;
workDir = pwd;
addpath(genpath(workDir));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_path = './data/';
data_name = 'data4seq';
data_file = strcat(data_path, data_name, '.mat'); 
data = load(data_file);
X = data.Xgen;
groundTruth_file = strcat(data_path, data_name, 'GroundTruth', '.csv');
posGrandTruth = readmatrix(groundTruth_file); % The grand truth positions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = 3; % the number of time lags (dynamic order)
L = 3; % the dimension of latent vector
lambda = 0.003; % lasso regularization coefficient
iniLen = 50; % initial number of samples in each subsequence
maxNoS = 4; % maximum number of subsequences

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Segmentation based on DiPCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
costPara = {s; L; lambda};
[pos, C, obj] = bottomUp(X, iniLen, maxNoS, @costDiPCA, costPara, 'max');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
strMethod = 'DiPCA';
save_path = './results/';
save_name = strcat('results_', data_name, '_NoS=', int2str(maxNoS), '_',...
    'method=', strMethod, '_', ...
    string(datetime('now','TimeZone','local','Format','dd-MM-yy-HHmm')));
save_file = strcat(save_path, save_name, '.mat');
save(save_file, 'pos', 'C', 'obj', 's', 'L', 'lambda', 'iniLen', 'maxNoS');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(save_file);
XtickStep = 500;
posSet = {pos, 'DiPCASeg'};
plotSegResults(X, posGrandTruth, posSet, XtickStep);
savefig(gcf, strcat(save_path, save_name, '.fig'));


% %% load data
% % defined in Section III.A.
% clc; clear;
% load('data4seq.mat')
% 
% %% initialize
% m = 5; % the dimension of sample vector
% n = 1000; % the number of samples in each sequence 
% s = 3; % the number of time lags (dynamic order)
% L = 3; % the dimension of latent vector
% lambda = 0.05;
% 
% var_v = 1; % variance of innovation
% var_e = 0.1; % variance of observation error
% 
% vkp = cell(2);
% vkp{1} = zeros(1, L); % mu of v_k
% vkp{2} = var_v*eye(L); % sigma of v_k
% 
% ekp = cell(2);
% ekp{1} = zeros(1, m); % mu of e_k
% ekp{2} = var_e*eye(m); % sigma of e_k
% 
% c = zeros(1, L);
% 
% P = [0.4315 0.1723 -0.0574;
%     0.1202 -0.1463 0.5348;
%     0.2483 0.1982 0.4797;
%     0.1151 0.1557 0.3739;
%     0.2258 0.5461 -0.0424];
% 
% iniLen = 50; % initial number of samples in each subsequence
% maxNoS = 4; % maximum number of subsequences
% 
% 
% %% segmentation
% NN = 6; % repeatation times
% 
% costPara = cell(2,1);
% costPara{1} = s;
% costPara{2} = L;
% costPara{3} = lambda;
% seg = zeros(NN, maxNoS-1);
% 
% parfor i=1:NN
%     i
%     [pos, C, obj] = bottomUp(Xgen, iniLen, maxNoS, @costDiPCA, costPara, 'max');
%     seg(i,:) = pos;    
% end
% save('data4seq_results_NoS=4_NN=6_lambda=0.05.mat', 'seg');

