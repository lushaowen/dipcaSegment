function [ t, x, v_k, e_k ] = genVarData( n, s, A, c, vkp, P, ekp)
% To generate data using a VAR process according to the following:
%     t_k = c + A^1*t_{k-1} + A^2*t_{k-2} + ... + A^L*t_{k-L} + v_k
%     x_k = P*t_k + e_k
%
% Written by : Shaowen Lu, Shuyu Huang @ Northeastern University, 2019.
%
% INPUTS
%   n   -   the number of samples
%   s   -   the time lag
%   A   -   is a cell array, each element contains A^{i} for the ith time lag
%   c   -   bias 
%   P   -   loading matrix
%   vkp -   parameters of the MA-dimensional gaussian process [mu, sigma]
%   ekp -   parameters of the MP-dimensional gaussian process [mu, sigma]
%
% OUTPUTS
%   t   -  the time series of latent variable t_k
%   x   -  the time series data x_k
%


% the dimension of each sample
L = size(A{1}, 1); % the dimension of latent vector
m = size(P, 1); % the dimension of data vector

% initialization
t = zeros(n, L);
t(1:s, :) = randn(s, L);
v_k = mvnrnd(vkp{1}, vkp{2}, n);

x = zeros(n, m);
e_k = mvnrnd(ekp{1}, ekp{2}, n);

for i=s+1:n
    Sigma_At = zeros(1, L);
    for j=1:s
        Sigma_At = Sigma_At + t(i-j,:)*A{j}';
    end
    t(i,:) = c + Sigma_At + v_k(i,:);
    x(i,:) = t(i,:)*P' + e_k(i,:);
end

t = t(s+1:end,:);
x = x(s+1:end,:);
e_k = e_k(s+1:end,:);
v_k = v_k(s+1:end,:);

end

