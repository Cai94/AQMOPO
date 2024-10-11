function zmin = UpdateIdealPoint(Pop, zmin_prev)
if ~exist('prev_zmin', 'var') || isempty(zmin_prev)
    zmin_prev = inf(size(Pop(1).Cost));
end
zmin = zmin_prev;
for i = 1: numel(Pop)
    zmin = min(zmin, Pop(i).Cost);
end
end