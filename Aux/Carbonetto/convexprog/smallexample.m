% This short script demonstrates the use of the interior-point solver to
% compute the solution to a quadratic program with convex objective
% (i.e. positive-definite Hessian) and convex, quadratic inequality
% constraints. More precisely, it finds the solution to the following
% optimization problem:
%
%   minimize    (1/2)x'Hx + q'x
%   subject to  ci(x) < b
%
% where the inequality constraints are quadratic functions
%
%   ci(x) = (1/2)x'P{i}x + r{i}'x
%
% and the quantities {H,q,P,r,b} are all specified below. Note that this
% code is not a particular efficient way to optimize a constrained quadratic
% program, and should not be used for solving large optimization problems.
%
% This particular example originally comes from the book: H. P. Schwefel
% (1995) Evolution and Optimum Seeking. The minimium occurs at (0,1,2,-1).
%
%                                         Peter Carbonetto
%                                         Dept. of Computer Science
%                                         University of British Columbia
%                                         Copyright 2008
clear

% These specify the quadratic objective function.
H = diag([ 2 2 4 2 ]);
q = [ -5 -5 -21 7 ]';

% These specify the quadratic inequality constraints.
P{1} = diag([ 4 2 2 0 ]);
P{2} = diag([ 2 2 2 2 ]);
P{3} = diag([ 2 4 2 4 ]);

r{1} = [  2 -1  0 -1 ]';
r{2} = [  1 -1  1 -1 ]';
r{3} = [ -1  0  0 -1 ]';

b = [ 5 8 10 ]';

% Set the starting point.
x0 = [ 0 0 0 0 ]';

% Solve the optimization problem with the primal-dual interior point solver.
fprintf('Running exact interior point method to obtain solution.\n');
data = { H q P r b };   % All the information about the quadratic program.
z    = zeros(size(b));  % A dummy variable; does nothing.
x    = ipsolver(x0,@(x)qprog(x,z,data,'objective'),...
		@(x)qprog(x,z,data,'gradient'),...
		@(x)qprog(x,z,data,'constraints'),...
		@(x,z)qprog(x,z,data,'jacobian'),...
		'bfgs',1e-6,100,true);
fprintf('\nSolution:\n');
disp(x);
