function Params = PerformScalarizing(z, Params)
n_obj = size(z, 1); 
n_pop = size(z, 2);
if ~isempty(Params.smin)
    zmax = Params.zmax;
    smin = Params.smin;
else
    zmax = zeros(n_obj, n_obj);
    smin = inf(1,n_obj);
end

for j = 1: n_obj
    Vector_scalar = GetScalarizingVector(n_obj, j);
    s = zeros(1, n_pop);
    for i = 1: n_pop
        s(i) = max(z(:, i) ./ Vector_scalar);
    end
    [sminj, ind] = min(s);
    if sminj < smin(j)
        zmax(:, j) = z(:, ind);
        smin(j) = sminj;
    end
end
Params.zmax = zmax;
Params.smin = smin;
end

function Vector_scalar = GetScalarizingVector(nObj, j)
epsilon = 1e-10;
Vector_scalar = epsilon * ones(nObj, 1);
Vector_scalar(j) = 1;
end