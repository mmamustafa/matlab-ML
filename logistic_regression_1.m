%{
Logistic regression binary classification Example.

x \in R^n
y \in [0, 1]

y_prob = Pr(y=1|x;a) = 1 / (1 + exp(-(a'x + a0)))

logistic regression is to find a \in R^n, a0 \in R.

Cost (loss) function w.r.t a:
cost(a; x1, x2, ..., y1, y2, ...) = -sum_{i=1} yi*log(y_prob) + (1-yi)*log(1-y_prob)
%}

clc
% (1) Generate binary dataset
%rand('seed',15);   randn('seed',1);

N = 1e2;              % number of points to generate
x_domain = [2, 10];   % domain of x values
x_0 = -(x_domain(1) + diff(x_domain) * rand);   % center value where y=0.5
x_1 = rand*1e1;     % paramters controlling the variance of bernoulli distribution (the higher the steeper)
a_true = [1; x_0] * x_1 * sign(rand - 0.5);          % logistic parameters
x = x_domain(1) + diff(x_domain) * rand(1,N);
y_true_prob = 1 ./ (1 + exp(-a_true'*[x; ones(1,size(x,2))]));
y_true = rand(1,N) < y_true_prob;

% (x) Consider N-d...


% (x) Split to train/test/validate dataset...


% (x) Estimate parameters using SGD
[a_pred, loss, iterations] = SGD_LR(x,y_true,randn(size(a_true)),1e-2,30)

% (x) Prediction after training
y_pred_prob = 1 ./ (1 + exp(-a_pred'*[x; ones(1,size(x,2))]));

% (x) Imbalanced dataset... Consider 0.5 threshold
y_pred = y_pred_prob >= 0.5;
[accuracy, precision, recall, f1] = binary_classifier_performance(y_true,y_pred)

% (x) Plots
figure(1), clf
subplot(2,1,1)
plot(x,y_true_prob,'b.', 'MarkerSize',10)

hold on
plot(x,y_pred_prob,'r.', 'MarkerSize',10)
hold off

grid on
xlabel('x'), ylabel('Pr(y)')
subplot(2,1,2)
plot(x,y_true,'x', 'MarkerSize',10)
xlabel('x'), ylabel('y')
title('Generate dataset')

figure(2), clf
precision_recall_curve(y_true, y_true_prob);
figure(3), clf
precision_recall_curve(y_true, y_pred_prob,'r-');
