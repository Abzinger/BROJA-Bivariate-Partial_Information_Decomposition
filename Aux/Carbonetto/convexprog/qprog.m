function varargout = qprog (x, z, data, what)
  [ H q P r b] = deal(data{:});
  n            = length(x);  % The number of (primal) variables.
  m            = length(b);  % The number of inequality constraints.
  
  if strcmp(what,'objective')
    
    % Compute the response of the quadratic objective function at x.
    f = x'*H*x/2 + q'*x;
    varargout = { f };
  elseif strcmp(what,'gradient')
    
    % Compute the gradient of the quadratic objective function at x. And, if
    % requested, output the Hessian of the quadratic objective.
    g = H*x + q;  
    varargout{1} = g;
    if nargout == 2
      varargout{2} = H;
    end
  elseif strcmp(what,'constraints')
    
    % Compute the response of the vector-valued inequality constraint
    % function. Repeat for each inequality constraint.
    c = zeros(m,1);
    for i = 1:m
      c(i) = x'*P{i}*x/2 + r{i}'*x - b(i);
    end
    varargout = { c };
  elseif strcmp(what,'jacobian')
    
    % Compute the m x n Jacobian of the vector-valued constraint function,
    % and the the n x n Hessian of the Lagrangian (minus the Hessian of the
    % objective). Repeat for each inequality constraint.
    J = zeros(m,n);
    W = zeros(n);
    for i = 1:m
      J(i,:) = (P{i}*x + r{i})';
      W      = W + z(i)*P{i};
    end
    varargout = { J W };
  end