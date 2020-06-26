function [ T, P, W, Beta, Theta_hat, J ] = dipca( X, s, L, lambda )
% This function implements the DiPCA algorithm introduced in paper: 
%   "Yining Dong and S. Joe Qin, A novel dynamic PCA algorithm for dynamic
%   data modeling and process monitoring, JPC, 2018"
%   
% Written by : Shaowen Lu, Shuyu Huang @ Northeastern University, 2019.
%
% INPUTS
%   X       -   a data matrix of 'nxm', each row repressenting a sample
%               X=[x1;x2;...;xn], where 'xi' is a sample vector of m variables
%   s       -   the number of time lags (dynamic order)
%   L       -   the number of dynamic latent variables
%   lambda  -   L1 regularization coefficient
%
% OUTPUTS
%   T       -   each row stores one of the L latent score vectors
%   P       -   each row stores one of the L loading vectors
%   W       -   each row stores one of the L w vectors
%   Beta    -   each row stores one of the L beta vectors
%   Theta_hat - least square estimates of the coefficients of T
%   J       -   a row vector stores the maximum objective values of each LV
%


%% main algorithm
%   beta: is a vector [beta_1 beta_2 ... beta_s]', where t_k =
%         beta_1*t_{k-1} +...+ beta_s*t_{k-s} + r_k, see eq.(4), where r_k
%         is the residual
%
%   w:    an 'mx1' vector denoting the weight for t_k = x_k'*w,  see eq.(5),
%         where x_k is the sample vector at time k, t_k is the latent variable
%         expressed as a linear combination of the past LV
%

% constants
THRESHOLD_DLV = 1e-4; % threshold for DLV iteration
MAX_ITERATION_DLV = 10; % maximum iteration number for DLV iteration

%%% step 1: 
% Scale X to zero-mean and unit-variance.
[n,m] = size(X); % X has 'n' samples, and each sample has dimension 'm'
X = zscore(X, 0, 1);

% Initialize
T = zeros(n,L); 
P = zeros(m,L); 
W = zeros(m,L); 
Beta = zeros(s,L); 
J = zeros(1,L);


%%% step 2-4: extract L latent variables

for i=1:L
    
    [Zs, Xmat] = formLagMat(X, s); % Form s-lagged matrix Zs from X
    
    % step 2: extract one LV
    Jmax = 0; tmax=zeros(n,1); pmax=zeros(m,1); wmax=zeros(m,1); beta_max=zeros(s,1);
    for j=1:m % try m times, and return w and beta which maximize J
       
        % initialize w and beta
        w = rand_unit_vector(m,1); % Initialize w to a random unit vector.
        %w = zeros(m,1); 
        w(j)=1; % each time assign a direction along coloumn to 1
        beta = zeros(s,1);
        
        % Calculate w and beta until convergence
        delta=1; iter=0;
        while(delta>THRESHOLD_DLV && iter<MAX_ITERATION_DLV ) % extracting one latent variable        
            t = X*w; % calculate the latent scores
            [ts, tmat] = formLagMat(t, s); % form s-lagged matrix ts from t            
            beta = ts'* tmat{s+1};
            beta = beta/norm(beta);
            w_new = zeros(m,1);
            for k=1:s
                X_term = Xmat{s+1}'*tmat{k} + Xmat{k}'*tmat{s+1};
                w_new = w_new + beta(k)*X_term;
            end
            w_new = w_new/norm(w_new);
            delta = norm(w-w_new);
            w = w_new;
            iter = iter+1;
        end % end while        
        
        t = X*w; % update the latent scores
        p = X'*t/(t'*t); % calculate the loading vector, see eq.(18)
        
        % Calcuate the objective J
        L1_reg = lambda*norm(beta);
        obj = abs(w'*Xmat{s+1}'*Zs* kron(beta,w)) / (n-s) - L1_reg; % objective, see eq.(9) 
        
        % cache t, p, w, beta which maximizing J
        if obj>Jmax
            Jmax = obj; tmax = t; pmax = p; wmax = w; beta_max = beta;            
        end             
    end % end for j=1:m
    
    % step 3: deflation
    X = X - t*p';
    
    % step 4: update and extract the next LV
    T(:,i) = tmax; 
    P(:,i) = pmax; 
    W(:,i) = wmax; 
    Beta(:,i) = beta_max;
    J(i) = Jmax;    
    
end % end for i=1:L


%%% step 5: Dynamic inner modeling. 
% Build a VAR model for latent scores based on (21) and (22).
[Ts_bar, Tmat] = formLagMat(T, s); % Form s-lagged matrix Zs from X

% The least squares estimate for Theta, see eq.(21)
Theta_hat = pinv(Ts_bar'*Ts_bar)*Ts_bar'*Tmat{s+1}; 

% Predicted T_{s+1}, see eq.(22)
T_splus1_hat = Ts_bar*Theta_hat; 

    
end


