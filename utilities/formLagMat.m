function [ Zs, Xmat ] = formLagMat( X, s )
% To form the stacked s-lag matrices from the given historical data matrix.
%
% Written by : Shaowen Lu, Shuyu Huang @ Northeastern University, 2019.
%
% INPUTS
%   X       -   histroical data matrix 
%   s       -   time lag
%
% Output
%   Zs      -   [X1 X2 ... Xs] is block matrix formed by Xi
%   Xmat    -   [x_i x_{i+1} ... x_{R+i-1}]' for i=1,2,...,s+1
%

[N, M] = size(X);
N = N-s;
Xmat = cell(1,s+1); 
Zs = zeros(N,M*s); 
for i=1:s
    Xmat{i} = X(i:N+i-1,:);
    Zs(:,((i-1)*M+1):(i*M)) = Xmat{i}; 
end
Xmat{s+1} = X(s+1:N+s,:);
end

