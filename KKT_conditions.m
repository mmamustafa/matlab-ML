np = 25;
x_lim = [-1, 4];
lambda_lim = [0, 5];

x = linspace(x_lim(1), x_lim(2), np);
lambda = linspace(lambda_lim(1), lambda_lim(2), np);

[x1, l1] = meshgrid(x, lambda);
L = (x1-1).^2 -x1 + l1.*(x1 - 1);


%figure(1)
%surf(x1, l1, L), xlabel('x'), ylabel('\lambda'), zlabel('L'), axis equal

x_min = x(1);
x_max = x(end);
l_min = lambda(1);
l_max = lambda(end);
y_min = min(L(:));
y_max = max(L(:));
figure(2)
for i=2:length(x)
  subplot(1,3,1)
  plot(x, L(i,:), '-')
  if i==2
    hold on
    axis equal
    axis([x_min, x_max, y_min, y_max])
    ylabel('L')
    grid on
  end
  xlabel(['x = ' num2str(x(i))])
  
  subplot(1,3,2)
  plot(lambda, L(:,i), '-')
  if i==2
    hold on
    axis equal
    axis([l_min, l_max, y_min, y_max])
    ylabel('L')
    grid on
  end
  xlabel(['\lambda = ' num2str(lambda(i))])
  
  subplot(1,3,3), cla
  surf(x1(1:i, 1:i), l1(1:i, 1:i), L(1:i, 1:i))
  xlabel('x'), ylabel('\lambda'), zlabel('L'),
  axis equal
  axis([x_min, x_max, l_min, l_max, y_min, y_max])
  pause(0.2)
end
subplot(1,3,1)
hold off
xlabel('x')
subplot(1,3,2)
hold off
xlabel('\lambda')
