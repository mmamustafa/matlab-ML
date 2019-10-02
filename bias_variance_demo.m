clc, close all
%rand('seed',0), randn('seed',0)

% generate the data
n = 100;
lims1 = [-1,1];
x = lims1(1):diff(lims1)/n:lims1(2);
% generate random polynomial with known order
sys_ord = 5;    % system order
lims2 = [-1,1];
coeff = lims2(1) + diff(lims2)*rand(1, sys_ord+1);
t = evaluate_poly(x, coeff);
sigma = 0.1;
noise = sigma * randn(size(x));
y = t + noise;

% split the data to test and train
p_train = 0.5;    % probability of split
n_train = round(n*p_train);
t1 = randperm(n);
ind_train = t1(1:n_train);
ind_test = setdiff(1:n, ind_train);
x_train = x(ind_train);
y_train = y(ind_train);
x_test = x(ind_test);
y_test = y(ind_test);

% generate error curves
n_orders = 20;
err_train_curve = zeros(1, n_orders);
err_test_curve = err_train_curve;
for o=1:n_orders
  coeff_hat = polyfit(x_train,y_train,o);
  y_train_hat = evaluate_poly(x_train, coeff_hat);
  err_train_curve(o) = sqrt(mean((y_train_hat - y_train).^2));
  y_test_hat = evaluate_poly(x_test, coeff_hat);
  err_test_curve(o) = sqrt(mean((y_test_hat - y_test).^2));
end


var_curve = abs(err_test_curve - err_train_curve);

% plot
subplot(121)
plot(x,t,'k-')
hold on
plot(x_train,y_train,'r.','MarkerSize',12)
plot(x_test,y_test,'g.','MarkerSize',12)
hold off

subplot(122)
plot(1:n_orders, err_train_curve, 'r-', 'LineWidth', 2)
axis([0,n_orders,0,1.2*sigma^0.5])
hold on
plot(1:n_orders, err_test_curve, 'b-', 'LineWidth', 2)
plot(1:n_orders, var_curve, 'g--', 'LineWidth', 2)
hold off
xlabel('polynomial order')
ylabel('error')
legend('training','testing','varaince')