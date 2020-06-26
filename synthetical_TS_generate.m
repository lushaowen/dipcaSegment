%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This routine generate a synthetical time series for segementation based
% on dynamical predictability
%
% Written by : Shaowen Lu, Shuyu Huang @ Northeastern University, 2019.
%


%% A time series with 4 sequences 'data4seq.mat'
% Defined in Secction III.A. in the paper
%   S. Lu and S. Huang, "Segmentation of multivariate industrial time 
%   series data based on dynamic latent variable predictability," 
%   IEEE Access, vol. 8, pp. 112092¨C112103, 2020.

clc; clear;

data_path = './data/';
data_name = 'data4seq';
ground_truth_name = strcat(data_name, 'GroundTruth');
mat_file = strcat(data_path, data_name, '.mat');
csv_file = strcat(data_path, data_name, '.csv');
groundTruth_file = strcat(data_path, ground_truth_name, '.csv');

m = 5;       % the dimension of sample vector
n = 1000;    % the number of samples in each sequence 
s = 3;       % the number of time lags (dynamic order)
L = 3;       % the dimension of latent vector
var_v = 1;   % variance of innovation
var_e = 0.1; % variance of observation error
ground_truth = [1000; 2000; 3000]; % the true segmentation posistions

vkp = cell(2);
vkp{1} = zeros(1, L); % mu of v_k
vkp{2} = var_v*eye(L); % sigma of v_k

ekp = cell(2);
ekp{1} = zeros(1, m); % mu of e_k
ekp{2} = var_e*eye(m); % sigma of e_k

c = zeros(1, L);

P = [0.4315 0.1723 -0.0574;
    0.1202 -0.1463 0.5348;
    0.2483 0.1982 0.4797;
    0.1151 0.1557 0.3739;
    0.2258 0.5461 -0.0424];
A1 = [0.1395 0.0958 0.0256;
    0.0248 0.1307 0.0466;
    0.327 0.1005 0.0431];
A2 = [0.1488 0.0194 0.0364;
    0.0002 0.0067 0.0440;
    0.0622 0.0410 0.1628];
A3 = [0.0674 0.0126 0.1594;
    0.1015 0.0967 0.1291;
    0.0638 0.1064 0.1193];
A4 = [0.1132 0.0005 0.1198;
    0.0531 0.1595 0.0756;
    0.1170 0.0114 0.0388];
A5 = [0.0205 0.0239 0.0843;
    0.0580 0.1317 0.1511;
    0.1158 0.1260 0.0573];
A6 = [0.1273 0.0933 0.1589;
    0.1024 0.0429 0.0385;
    0.1128 0.1599 0.0853];

% segment #1
S1 = cell(s);
S1{1} = A1;
S1{2} = A2;
S1{3} = A3;

% segment #2
S2 = cell(s);
S2{1} = -A1;
S2{2} = -A2;
S2{3} = -A3;

% segment #3
S3 = cell(s);
S3{1} = A4;
S3{2} = A5;
S3{3} = A6;

% segment #4
S4 = cell(s);
S4{1} = -A4;
S4{2} = -A5;
S4{3} = -A6;

% note that the number of samples is n+s
[Tgen1,Xgen1,Vgen1,Egen1] = genVarData( n+s, s, S1, c, vkp, P, ekp ); 
[Tgen2,Xgen2,Vgen2,Egen2] = genVarData( n+s, s, S2, c, vkp, P, ekp ); 
[Tgen3,Xgen3,Vgen3,Egen3] = genVarData( n+s, s, S3, c, vkp, P, ekp ); 
[Tgen4,Xgen4,Vgen4,Egen4] = genVarData( n+s, s, S4, c, vkp, P, ekp ); 

Tgen = [Tgen1; Tgen2; Tgen3; Tgen4];
Xgen = [Xgen1; Xgen2; Xgen3; Xgen4];
Vgen = [Vgen1; Vgen2; Vgen3; Vgen4];
Egen = [Egen1; Egen2; Egen3; Egen4];

save(mat_file, 'Tgen', 'Xgen', 'Vgen', 'Egen', 'ground_truth', 's', 'L', 'var_v', 'var_e');
writematrix(Xgen, csv_file);
writematrix(ground_truth, groundTruth_file);



%% Data with 4 sequences 'data4seq2.mat' with different observation matrix

clc; clear;

data_path = './data/';
data_name = 'data4seq2';
ground_truth_name = strcat(data_name, 'GroundTruth');
mat_file = strcat(data_path, data_name, '.mat');
csv_file = strcat(data_path, data_name, '.csv');
groundTruth_file = strcat(data_path, ground_truth_name, '.csv');

m = 5;       % the dimension of sample vector
n = 1000;    % the number of samples in each sequence 
s = 3;       % the number of time lags (dynamic order)
L = 3;       % the dimension of latent vector
var_v = 1;   % variance of innovation
var_e = 0.1; % variance of observation error
ground_truth = [1000; 2000; 3000]; % the true segmentation posistions

vkp = cell(2);
vkp{1} = zeros(1, L); % mu of v_k
vkp{2} = var_v*eye(L); % sigma of v_k

ekp = cell(2);
ekp{1} = zeros(1, m); % mu of e_k
ekp{2} = var_e*eye(m); % sigma of e_k

c = zeros(1, L);

P1 = [0.4315 0.1723 -0.0574;
    0.1202 -0.1463 0.5348;
    0.2483 0.1982 0.4797;
    0.1151 0.1557 0.3739;
    0.2258 0.5461 -0.0424];
P2 = [0.1695 -0.0631 0.1065;
    -0.0450 0.6723 0.6273;
    -0.5631 -0.5148 -1.0799;
    0.8026 -0.2585 0.5440;
    -0.0900 0.4606 0.3706];
A1 = [0.1395 0.0958 0.0256;
    0.0248 0.1307 0.0466;
    0.327 0.1005 0.0431];
A2 = [0.1488 0.0194 0.0364;
    0.0002 0.0067 0.0440;
    0.0622 0.0410 0.1628];
A3 = [0.0674 0.0126 0.1594;
    0.1015 0.0967 0.1291;
    0.0638 0.1064 0.1193];
A4 = [0.1132 0.0005 0.1198;
    0.0531 0.1595 0.0756;
    0.1170 0.0114 0.0388];
A5 = [0.0205 0.0239 0.0843;
    0.0580 0.1317 0.1511;
    0.1158 0.1260 0.0573];
A6 = [0.1273 0.0933 0.1589;
    0.1024 0.0429 0.0385;
    0.1128 0.1599 0.0853];

% segment #1
S1 = cell(s);
S1{1} = A1;
S1{2} = A2;
S1{3} = A3;

% segment #2
S2 = cell(s);
S2{1} = A1;
S2{2} = A2;
S2{3} = A3;

% segment #3
S3 = cell(s);
S3{1} = -A1;
S3{2} = -A2;
S3{3} = -A3;

% segment #4
S4 = cell(s);
S4{1} = -A1;
S4{2} = -A2;
S4{3} = -A3;

% note that the number of samples is n+s
[Tgen1,Xgen1,Vgen1,Egen1] = genVarData( n+s, s, S1, c, vkp, P1, ekp ); 
[Tgen2,Xgen2,Vgen2,Egen2] = genVarData( n+s, s, S2, c, vkp, P1, ekp ); 
[Tgen3,Xgen3,Vgen3,Egen3] = genVarData( n+s, s, S3, c, vkp, P2, ekp ); 
[Tgen4,Xgen4,Vgen4,Egen4] = genVarData( n+s, s, S4, c, vkp, P2, ekp ); 

Tgen = [Tgen1; Tgen2; Tgen3; Tgen4];
Xgen = [Xgen1; Xgen2; Xgen3; Xgen4];
Vgen = [Vgen1; Vgen2; Vgen3; Vgen4];
Egen = [Egen1; Egen2; Egen3; Egen4];

save(mat_file, 'Tgen', 'Xgen', 'Vgen', 'Egen', 'ground_truth', 's', 'L', 'var_v', 'var_e');
writematrix(Xgen, csv_file);
writematrix(ground_truth, groundTruth_file);


