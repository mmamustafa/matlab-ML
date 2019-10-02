% Binary variable: x
% ------------------
clc, clear all
% (1) Define parameters
u = 0.2;    % p(x=1|u)
N = 1e3;    % number of data points
a0 = 20;    % beta distribution parameter
b0 = 2;    % beta distribution parameter
u_prior = a0/(a0+b0);
var_prior = a0*b0/((a0+b0)^2*(a0+b0+1));  % variance of the paramter u not the data

% (2) Generate x observations based on parameters
x = rand(1,N)<u;

% (3) MLE for u
m = sum(x);
u_ML = m/N;      % p(D)

%{
% (4a) MAP for u
u_MAP = zeros(1,N);
for n=1:N
  x_hat = x(1:n);
  m = sum(x_hat);
  N_hat = n;
  u_MAP(n) = (m+a0)/(a0+b0+N_hat);
end
% >>>>> u_MAP is the same as the one obtained below <<<<<
%}

% (4b) Sequential MAP for u
u_MAP = zeros(1,N);
var_MAP = zeros(1,N);
a = a0;   b = b0;
for n=1:N
  % update beta parameters
  a = a+x(n);
  b = b+1-x(n);
  u_MAP(n) = a/(a+b);
  var_MAP(n) = a*b/((a+b)^2*(a+b+1));
end
figure(1), clf
subplot(1,2,1)
plot(1:N,u_MAP,'b-', 'LineWidth',2)
xlabel('n')
title('\mu_{MAP} (Sequential)')
subplot(1,2,2)
plot(1:N,var_MAP,'b-','LineWidth',2)
xlabel('n')
title('var[\mu_{MAP}] (Sequential)')

% (5) Plot beta prior and posterior
figure(2), clf
u1 = 0:1e-3:1;
plot(u1, betapdf(u1, a0, b0), 'b--', 'LineWidth',2)
xlabel('\mu')
hold on
plot(u1, betapdf(u1, a, b), 'b-', 'LineWidth',2)
hold off
title(['N = ' num2str(N)])

format = 2;
legend({['prior (E[\mu] = ' num2str(u_prior, format) ', var[\mu] = ' num2str(var_prior, format) ')'],...
         ['posterior (E[\mu] = ' num2str(u_MAP(end), format) ', var[\mu] = ' num2str(var_MAP(end), format) ')']})

