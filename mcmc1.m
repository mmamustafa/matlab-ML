clear, clc

%rand('seed', 10);
a = 2;        % greater than 0
n_dist = 3;   % number of distributions
%{
pis = rand(1, n_dist);  % mixing parameters
pis = pis / sum(pis);   % normalize
mus = -a/2 + a*rand(1, n_dist);     % meansq
stds = abs(a/3*(-a + 2*a*rand(1, n_dist)));   % std
%}

pis = [0.2, 0.5, 0.3];
mus = [-1.5, 0, 1.5];
stds = [0.1, 0.5, 0.2];

% generate samples
x = -a:1e-2:a;
y = zeros(size(x));   % initialize
for i=1:n_dist
  y = y + pis(i) * normpdf(x, mus(i), stds(i));
end

sum(y)*1e-2
figure(1); clf
subplot(2,1,1)
plot(x, y)


% MCMC
n = 1e4;
x1 = zeros(1, n);
x1(1) = -a + 2*a*rand;
for i=1:n-1
  x_c = normrnd(x1(i), 0.1);    % Monte Carlo
  num = 0;
  den = 0;
  for j=1:n_dist
    num = num + pis(j) * normpdf(x_c, mus(j), stds(j));
    den = den + pis(j) * normpdf(x1(i), mus(j), stds(j));
  end
  % Markov Chain
  if rand < min(1, num/den)
    % accept
    x1(i+1) = x_c;
  else
    % reject
    x1(i+1) = x1(i);
  end
end

subplot(2,1,2)
h = histogram(x1,50);
t1 = diff(h.BinEdges);
subplot(2,1,1)
hold on
plot(h.BinEdges(1:end-1), h.Values/sum(h.Values)/t1(1), 'r--')
hold off
legend('true', 'estimate')