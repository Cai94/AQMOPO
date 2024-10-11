function [Pop, d, rho] = AssociateToReferencePoint(Pop, Params)
Zr = Params.Zr;
nZr = Params.nZr;
rho = zeros(1,nZr);
d = zeros(numel(Pop), nZr);

for i = 1:numel(Pop)
    for j = 1: nZr
        w = Zr(:, j) / norm(Zr(:,j));
        z = Pop(i).NormalizedCost;
        d(i,j) = norm(z - w' * z * w);
    end
    [dmin, jmin] = min(d(i, :));
    Pop(i).AssociatedRef = jmin;
    Pop(i).DistanceToAssociatedRef = dmin;
    rho(jmin) = rho(jmin) + 1;
end
end