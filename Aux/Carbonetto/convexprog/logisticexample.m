% This little script demonstrates the use the primal-dual interior-point
% solver to compute a logistic regression model for predicting the binary
% {0,1} outputs of a input vectors. It computes the set of parameters that
% maximizes the likelihood (or minimizes the sum of squared errors), but
% subject to a penalization on (otherwise known as the "Lasso" or "Basis
% pursuit denoising"). The best way to understand what is going on here
% is to go and read:
%
%    * Trevor Hastie, Robert Tibshirani and Jerome Friedman (2001). 
%      The Elements of Statistical Learning. Springer.
%
%    * Scott S. Chen, David L. Donoho and Michael A. Saunders (2001). 
%      Atomic Decomposition by Basis Pursuit. SIAM Review, Vol. 43, 
%      No. 1. (2001), pp. 129-159.
%
% The computed solution should be fairly close to the "true" regression
% coefficients (beta). Note that the Hessian in this approach is intensely
% ill-conditioned (due to transformation of the regression coefficients into
% their positive and negative components), so in general this may not be the
% best approach for logistic regression with L1 regularization. The steepest
% descent direction actually works well descent despite a Hessian with a
% very large condition number. There is probably a good reason why, but
% at this point I don't know. (Sorry.)
%
%                                         Peter Carbonetto
%                                         Dept. of Computer Science
%                                         University of British Columbia
%                                         Copyright 2008
clear

% CREATE DATA SET.
% Generate the input vectors from the standard normal, and generate the
% binary responses from the regression with some additional noise, and then
% transform the results using the logistic function. The variable "beta" is
% the set of true regression coefficients of length m.
n       = 100;   % The number of training examples.
epsilon = 0.25;  % Standard deviation in noise of outputs.
beta    = [ 0  0 2 -4 0 0 -1 3 ]';         % True regression coefficients.
sigma   = [ 10 1 1  1 1 1  1 1 ]';         % Standard deviation of coord's.
m       = length(beta);                    % Number of dimensions/features.
A       = repmat(sigma',n,1).*randn(n,m);  % The n x m matrix of examples.
noise   = epsilon*randn(n,1);              % Noise in outputs.
y       = rand(n,1) < logit(A*beta+noise); % The binary outputs.

% COMPUTE SOLUTION WITH INTERIOR-POINT METHOD.
% Compute the L1-regularized maximum likelihood estimator.
lambda = 1/2;          % Level of L1 regularization (sparsity).
P      = [A -A];
x      = ones(2*m,1);  % The inital point.
z      = ones(2*m,1);  % A dummy variable; does nothing.
data   = { P y lambda };
x      = ipsolver(x,@(x)logisticl1(x,z,data,'objective'),...
		  @(x)logisticl1(x,z,data,'gradient'),...
		  @(x)logisticl1(x,z,data,'constraints'),...
		  @(x,z)logisticl1(x,z,data,'jacobian'),...
		  'steepest',1e-4,100,true);
w      = x(1:m) - x(m+1:end);
fprintf('\nSolution:\n');
disp(w);
