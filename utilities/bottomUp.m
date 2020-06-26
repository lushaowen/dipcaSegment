function [ pos, C, obj ] = bottomUp( X, iniLen, maxNoS, costFun, costPara, mode )
% This function implements the bottom up segmentation algorithm.
%
% Written by : Shaowen Lu, Shuyu Huang @ Northeastern University, 2019.
%
%  INPUTS
%   X           -   a data matrix of 'nxm', each row repressenting a sample
%                   of X=[x1;x2;...;xn], where 'xi' is a sample vector of m 
%                   variables
%   iniLen      -   the initial number of samples in each subsequence
%   maxNoS      -   the maximum number of sequences
%   costFun     -   the handle of merging cost calculation function
%   costPara    -   parameters of costFun which includes
%                   - s:  the number of time lags (dynamic order)
%                   - L:  the number of dynamic latent variables
%   mode        -   the mode of greedy search
%                   - 'min': residual oriented 
%                   - 'max': accuracy oriented
%
%  OUTPUTS
%   pos         -   locations of segmentation
%   C           -   subsequences
%   obj         -   objective
%

[n, ~] = size(X);

% Partition X into equal length subsequences
NoS = ceil(n/iniLen); % number of subsequences
rowDist = ones(NoS,1)*iniLen; 
if mod(n, iniLen) > 0
    rowDist(end) = mod(n, iniLen);
end
C = mat2cell(X, rowDist); % initial segments

obj = zeros(NoS, 1); % store the average cost

% Bottom up segmentation
while NoS>maxNoS
    % Calculate all the merging cost
    cost_merge = zeros(NoS-1,1);
    for i=1:NoS-1        
        cost_merge(i) = costFun(C{i}, C{i+1}, costPara);
    end
    
    switch mode
        case 'min'  % Find the minial merging cost (reconstruction error)
            index = find(cost_merge==min(cost_merge), 1, 'first');
        case 'max'  % Find the maximal merging cost (prediction accuracy)
            index = find(cost_merge==max(cost_merge), 1, 'first');
        otherwise
            index = find(cost_merge==max(cost_merge), 1, 'first');
    end
    % index = find(cost_merge==max(cost_merge), 1, 'first');
    
    % Merge C{index} and C{index+1}
    C{index} = vertcat(C{index}, C{index+1});
    C(index+1) = [];
    
    % Update
    NoS = length(C);
    
    % Average cost
    temp = zeros(NoS,1);
    for ic = 1:length(C)
        temp(ic) = costFun(C{ic}, C{ic}, costPara);
    end
    obj(NoS) = mean(temp);   
end % while 

% Segmentation locations
pos = zeros(NoS-1,1);
pos(1) = length(C{1});
for i=2:NoS-1
    pos(i) = pos(i-1) + length(C{i});
end

end

