function [ cost ] = costDiPCA( X1, X2, parameters )
% To calculate the merging cost of two subsequences X1 and X2
% The cost is defined as the prediction accuracy of using the reference model
% built using the longer sequence to predict the other sequence
%
% Written by : Shaowen Lu, Shuyu Huang @ Northeastern University, 2019.
%
% To cite: 
%   S. Lu and S. Huang, "Segmentation of multivariate industrial time 
%   series data based on dynamic latent variable predictability," 
%   IEEE Access, vol. 8, pp. 112092¨C112103, 2020.
%
%  INPUTS
%   X1, X2      -   the two data sequence of 'nxm', each row repressenting a 
%                   sample of X=[x1;x2;...;xn], where 'xi' is sample vector 
%                   of m variables
%   parameters  -   passing the following parameters
%                   s:  the number of time lags (dynamic order)
%                   L:  the number of dynamic latent variables
%                   lambba: regularization parameter
%
%  OUTPUTS
%   cost        -   caulcated cost
%

s = parameters{1};
L = parameters{2};
lambda = parameters{3};

% dipca model
[T, P1, W1, ~, Theta_hat1, ~] = dipca(X1, s, L, lambda);
[~, P2, W2, ~, Theta_hat2, ~] = dipca(X2, s, L, lambda);
reg = 0.01*eye(L);
R1 = W1/(P1'*W1+reg);
R2 = W2/(P2'*W2+reg);

% calculate the score of X1 on model_2 and X2 on model_1
T_X1M2 = X1*R2;
T_X2M1 = X2*R1;

% form s-lagged matrix Ts from T
[Ts1, Tmat1] = formLagMat(T_X1M2, s); 
[Ts2, Tmat2] = formLagMat(T_X2M1, s); 

% predict T
That_X1M2 = Ts1*Theta_hat2; %estimate/predict T1 by model 2
That_X2M1 = Ts2*Theta_hat1; %estimate/predict T2 by model 1

% prediction accuracy based on the inner product of two matrices
acc1 = dotMat(Tmat1{s+1}, That_X1M2);
acc2 = dotMat(Tmat2{s+1}, That_X2M1);

% choose the longer sequence as a reference model to predict the other
if size(X1,1) > size(X2,1)
    cost = acc2;
elseif size(X1,1) < size(X2,1)
    cost = acc1;
else
    cost = (acc1+acc2)/2;
end

% When X1 is equal to X2, we compute the reconstruction accuracy.
if isequal(X1, X2)
    [Ts, Tmat] = formLagMat(T, s); 
    T_splus1 = Tmat{s+1}; 
    T_hat = Ts*Theta_hat1;
    cost = dotMat(T_splus1, T_hat);
end

end

