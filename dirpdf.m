function pdf = dirpdf(x,a)
pdf = exp(gammaln(sum(a)) - sum(gammaln(a)) + sum((a-1).*log(x)));