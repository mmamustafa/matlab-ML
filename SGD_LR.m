% Stochastic Gradient Descent for Logistic regression

function [a, loss, iterations] = SGD_LR(X,Y,a,learning_rate,batch_size,eps)
  if nargin < 2
    error('Minimum number of input parameter is 2!')
  end
  X = [X; ones(1,size(X,2))];
  [n_parameters, n_samples] = size(X);
  if nargin < 6
    eps = 1e-2;
    if nargin < 5
      batch_size = n_samples;
      if nargin < 4
        learning_rate = 1e-2;
        if nargin < 3
          a = randn(n_parameters,1);
        end
      end
    end
  end
  n0 = ceil(n_samples/batch_size);
  if n0 > n_samples
    n0 = n_samples;
  end
  ind_tmp = batch_size*[1:n0];
  ind_start = [1, ind_tmp(1:end-1)];
  ind_finish = ind_tmp - 1;
  ind_finish(end) = n_samples;
  old_loss = 1e6;
  for i = 1:5e3
    for j = 1:n0
      X_sub = X(:,ind_start(j):ind_finish(j));
      Y_sub = Y(:,ind_start(j):ind_finish(j));
      y_hat = logistic(X_sub,a);
      % Jacobian
      G = X_sub*(y_hat - Y_sub)';
      % Hessian
      %R = diag(y_hat .* (1 - y_hat));
      %H = X_sub * R * X_sub';

      a = a - learning_rate * G;
    end
    % stopping condition
    loss = cross_entropy_loss(X,Y,a);
    delta = abs(old_loss - loss);
    old_loss = loss;
    if delta < 1e-4
      break
    end
  end
  iterations = i;
return


function y_hat = logistic(X,a)
  y_hat = 1 ./ (1 + exp(-a'*X));
return

function loss = cross_entropy_loss(X,Y,a)
  y_hat = logistic(X,a);
  t0 = Y.*log(y_hat) + (1-Y).*log(1-y_hat);
  t0(~isfinite(t0)) = 0;
  loss = -sum(t0);
return