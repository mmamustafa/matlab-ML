% multinomial variable: x
% -----------------------
clc, clear all
% (1) Define parameters
K = 6;      % dimensions (min 2)
N = 2e3;    % number of data points
t0 = rand(K,1);
u = t0/sum(t0);    % p(x=1|u)
L = [1, 10];      % limits of alpha
a0 = L(1) + diff(L)*rand(K,1);  % Dirichlet distribution parameters

% prior parameters 
ao = sum(a0);
u_prior = a0./ao;
var_prior = u_prior.*(1-u_prior)./(ao+1);   % variance of the paramter u not the data

% (2) Generate x observations based on parameters
%     https://en.wikipedia.org/wiki/Multinomial_distribution#Sampling_from_a_multinomial_distribution
x = zeros(K,N);
[u_hat,ind0] = sort(u,'descend');
[~,ind0] = sort(ind0);  % for later to put in the same order
ind1 = K - sum((repmat(cumsum(u_hat), 1, N) - repmat(rand(1, N), K, 1)) >= 0) + 1;
x(sub2ind(size(x),ind1,1:N)) = 1;
% put in the original order of u
x = x(ind0,:);

% (3) MLE for u
m = sum(x,2);
u_ML = m/N;      % p(D)   D here is data

% (4b) Sequential MAP for u
u_MAP = zeros(K,N);
var_MAP = zeros(K,N);
a = a0;
for n=1:N
  % update dirichlet parameters
  a = a+x(:,n);
  
  ao = sum(a);
  u_MAP(:,n) = a./ao;
  var_MAP(:,n) = u_MAP(:,n).*(1-u_MAP(:,n))./(ao+1);
end
figure(1), clf
subplot(1,2,1)
h=plot(1:N,u_MAP, 'LineWidth',2);
xlabel('n')
title('\mu_{MAP} (Sequential)')
subplot(1,2,2)
plot(1:N,var_MAP,'LineWidth',2)
xlabel('n')
title('var[\mu_{MAP}] (Sequential)')

% (5) Plot beta prior and posterior
figure(2), clf
u1 = 0:1e-3:1;
colors = get(h,'Color');
for k=1:K
  % >>> Note that this should be Dirichlet pdf not beta pdf <<<<<<<<
  % >>> Parameters here are decoupled for visualization purposes <<<
  plot(u1, betapdf(u1, a0(k), sum(a0)-a0(k)), '--', 'LineWidth',2,'Color',colors{k})
  if k==1, hold on, end
  plot(u1, betapdf(u1, a(k), sum(a)-a(k)), '-', 'LineWidth',2,'Color',colors{k})
end
xlabel('\mu')
title(['N = ' num2str(N)])
