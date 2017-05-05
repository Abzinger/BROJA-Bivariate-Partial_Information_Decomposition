% Return the response of the logistic function at x.
function u = logit (x)
  u = 1./(1 + exp(-x));