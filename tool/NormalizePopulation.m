function [Pop, Params] = NormalizePopulation(Pop, Params)
Params.zmin = UpdateIdealPoint(Pop, Params.zmin);
fp = [Pop.Cost] - repmat(Params.zmin, 1, numel(Pop));
Params = PerformScalarizing(fp, Params);
a = FindHyperplaneIntercepts(Params.zmax);
for i = 1:numel(Pop)
    Pop(i).NormalizedCost = fp(:,i) ./ a;
end
end

function a = FindHyperplaneIntercepts(zmax)  
w = ones(1, size(zmax,2)) / zmax;
a = (1 ./ w)';
end