function [accuracy, precision, recall, f1_score] = binary_classifier_performance(y_true,y_pred)
  N = length(y_true);
  tp = sum(y_true & y_pred);     % true positive
  tn = sum(~y_true & ~y_pred);   % true negative
  fp = sum(~y_true & y_pred);    % false positive (type 1 error)
  fn = sum(y_true & ~y_pred);    % false negative (type 2 error)
  accuracy = (tp + tn) / N;
  precision = tp / (tp + fp);
  recall = tp / (tp + fn);
  if precision==0 && recall==0
    f1_score = 0;
  else
    f1_score = 2 * precision * recall / (precision + recall);
  end
return