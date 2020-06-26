function ret = dotMat( A, B )
% To calcuate the inner product of two matrices <A, B>
% The definition of matrices inner product is seen in 
%   Matrix Algebra: Theory, Computations, and Applications in Statistics by
%   Gentle, 2007, pp.132
%
% Written by : Shaowen Lu, Shuyu Huang @ Northeastern University, 2019.
%
%  INPUTS
%   A, B    -   the input matrices 
%
%  OUTPUTS
%   ret     -   the inner product of A and B
%

ret = trace(A' * B) / (norm(A,'fro')*norm(B,'fro')); 
end

