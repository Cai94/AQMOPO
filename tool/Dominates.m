function flag = Dominates(x, y)
if isstruct(x)
    x = x.Cost;
end
if isstruct(y)
    y = y.Cost;
end
flag = all(x <= y) && any(x < y);
end