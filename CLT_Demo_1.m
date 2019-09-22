function CLT_Demo_1()
% Central Limit Theorem CLT Demo (Probabilities).
% Sum of many (nv --> inf) independent random variables is normal.
% In this test, we add nv random variables generated from different
% distributions as shown in Fig.1.
clc
fid = 0;
fn = 'Times';   fz = 15;    col = 'k-';     w = 1.5;

ns = 1e4;       % number of samples to generate for each random variable
nv = 1e3;         % number of random variables per distribution


% (1) Plot the distribution of random variables.
fid = fid+1;    figure(fid),    clf,    set(gcf, 'color', [1 1 1])
subplot(1,3,1),hist(rand_quad(1,ns),50)
title(['power distribution histogram'],'Interpreter','latex','FontSize',fz)
subplot(1,3,2),hist(rand(1,ns),50)
title(['uniform distribution histogram'],'Interpreter','latex','FontSize',fz)
subplot(1,3,3),hist(double(rand(1,ns)<=0.3),50)
title(['bernoulli p=0.3 distribution histogram'],'Interpreter','latex','FontSize',fz)


% (2) Add nv samples generated from each random variables defined by
%     distributions shown in Fig.1
sum1 = zeros(1,ns);     sum2 = sum1;        sum3 = sum1;
% >>> comment any combination of below to see the effect of each.
sum1 = sum(rand(nv,ns)<=0.3, 1);	% bernoulli distribution
sum2 = sum(rand(nv,ns), 1);         % uniform distibution
sum3 = sum(rand_quad(nv,ns), 1);    % power distribution

% (3) Plot the distribution of the SUM of random variables.
fid = fid+1;    figure(fid),    clf,    set(gcf, 'color', [1 1 1])
hist((sum1+sum2+sum3)/(3*nv),50)
title(['Histogram of adding nv random variables'],'Interpreter','latex','FontSize',fz)
return
% -------------------------------------------------------------------------
function s = rand_quad(ns_row,ns_col)
% generates ns_row x ns_col samples from power distribution
s = (125*rand(ns_row,ns_col)).^(1/3)-1;     % inverse transform sampling
return