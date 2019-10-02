function auc = precision_recall_curve(y_true, y_prob, cs)
 if nargin<3
      cs = 'b-';
  end
  N = length(y_true);
  thresholds = sort(y_prob);    % all thresholds
  accuracies = zeros(1,N);
  precisions = zeros(1,N);
  recalls = zeros(1,N);
  f1s = zeros(1,N);
  for i=1:N
    y_pred = y_prob >= thresholds(i);
    [accuracies(i), precisions(i), recalls(i), f1s(i)] = binary_classifier_performance(y_true,y_pred);
  end
  % make sure they start at (1,0) and finish at (0,1)
  recalls = [1, recalls, 0];
  precisions = [0, precisions, 1];
  f1s = [0, f1s, 0];
  thresholds = [0, thresholds, 1];
  % Area under the curve (trapezpoid method)
  w = -diff(recalls);
  h = (precisions(1:end-1) + precisions(2:end)) / 2;
  auc = sum(h.*w);
  % Optimal threshold that maximizes F1 score
  [f1_max,ind] = max(f1s);
  threshold_opt = thresholds(ind);  
  % plot
  subplot(2,1,1)
  plot(recalls, precisions, cs, 'LineWidth', 2)
  axis equal
  axis([0 1 0 1])
  xlabel('recall'), ylabel('precision')
  title({'Precision-recall curve', ['AUC = ' num2str(auc)]})
  legend
  subplot(2,1,2)
  plot(thresholds, f1s, cs, 'LineWidth', 2)
  hold on
  plot([threshold_opt, threshold_opt], [0, f1_max], 'k--', 'LineWidth', 1.5)
  hold off
  xlabel('threshold'), ylabel('F1 score')
  title({'Optimal threshold maximizing F1 score', ...
         ['(threshold = ' num2str(threshold_opt) ',  f1 = ' num2str(f1_max) ')']})
return