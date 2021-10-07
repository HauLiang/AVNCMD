function x = Estimate_NCS(A, D, g, K, N)
%
% This code implements the Section 3.1 "Estimating the Nonlinear Chirp Signal" 
%
% The optimization problem may be expressed as
%    minimize   alpha*|| g - Ax ||_2^2 + || Dx ||_2^2
%
%
% Inputs:
%    A:  dictionary
%    D:  2K block second-order difference matrix
%    g:  sampled signal
%    K:  number of the modes
%    N:  number of the samples
% Outputs:
%    x:  solution of the above optimization problem
%
% Author: Hao Liang
% Last modified by: 21/10/02
%

% Find a matrix M whose rows are orthogonal to those in D
M = zeros(4*K,2*K*N); ii = 1;
for i = 1:2*K
    M(ii,i*N-1) = 1;
    M(ii+1,i*N) = 1;
    ii = ii+2;
end

% Construct a full rank matrix D_tilde
D_tilde = [D;M];

% Parameter setting
Lambda = A/D_tilde;
Lambda1 = Lambda(:,1:2*K*(N-2));
Lambda2 = Lambda(:,2*K*(N-2)+1:end);
theta_matrix = (Lambda2'*Lambda2)\Lambda2';
w_matrix = eye(N,N) - Lambda2*theta_matrix;
y = w_matrix*g;
Phi = w_matrix*Lambda1;

% Update w
w = Bayesian_strategy(Phi, y);

% Update theta
theta = theta_matrix*(g-Lambda1*w);

% Reconstruct x
x = D_tilde\[w;theta];

end