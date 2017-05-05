function varargout = logisticl1 (x, z, data, what)
  [ P y lambda ] = deal(data{:});
  m = length(x);   % The number of primal variables.
  u = logit(P*x);  % Responses of the logistic function.
    
  if strcmp(what,'objective')
    
    % Compute the response of the logistic loss function with a penalty
    % on the L1 norm of the regression coefficients.
    f = -sum(y.*log(u) + (1-y).*log(1-u)) + lambda*sum(x);
    varargout = { f };
  elseif strcmp(what,'gradient')
    
    % Compute the gradient of the objective function. And, if requested,
    % output the Hessian of the objective.
    g = -P'*(y - u) + lambda;
    varargout{1} = g;
    if nargout == 2
      U = diag(sparse(u.*(1-u)));
      H = P'*U*P;
      varargout{2} = H;
    end
  elseif strcmp(what,'constraints')
    
    % Compute the responses of the lower bounds (x > 0) on the optimization
    % variables.
    c = -x;
    varargout = { c };
  elseif strcmp(what,'jacobian')
    
    % Compute the Jacobian of the bound constraints. Since the
    % constraints are linear, the Hessian of the Lagrangian is zero.
    J = -speye(m);
    W = spzeros(m,m);
    varargout = { J W };
  end

% ------------------------------------------------------------------
% Return a sparse n x m matrix of zeros.
function A = spzeros (n, m)
  A = sparse([],[],[],n,m);
  