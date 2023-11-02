function y = mse(x1, x2)

se = (x1 - x2).*(x1- x2);
y = sum(se)/length(se);