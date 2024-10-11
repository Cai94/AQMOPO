function [Map_Fast, Ncoved_Fast, Centre_Fast, Link_Fast, Dis_Fast] = covFast(Zone_main, ModCov)
Box = 1 - Zone_main;
Zone_Fast = cell(size(Zone_main));
Map_Fast = {};
Ncoved_Fast = zeros(size(Zone_main));
Centre_Fast = [];
Link_Fast = [];

num_block = 1;
while ~isempty(find(Box > 0))
    Map_Fast{num_block} = zeros(size(Zone_main));
    Point_s = find(Box > 0);
    point = Point_s(randi(length(Point_s)));
    [sx,sy] = ind2sub(size(Box),point);
    Centre_Fast(num_block, :) = [sx, sy];
    Link_Fast = [Link_Fast, zeros([size(Link_Fast, 1), 1])];
    Link_Fast = [Link_Fast; zeros([1, size(Link_Fast, 2)])];
    for i = 1: size(ModCov, 1)
        for j = 1: size(ModCov, 2)
            if ~isempty(ModCov{i, j})
                dx = sx + ModCov{i, j}(1);
                dy = sy + ModCov{i, j}(2);
                if dx > 0 && dy > 0 && dx <= size(Zone_main, 1) && dy <= size(Zone_main, 2)
                    if Zone_main(dx, dy) == 0
                        Map_Fast{num_block}(dx, dy) = 1;
                        Zone_Fast{dx, dy} = [[Zone_Fast{dx, dy}], num_block];
                        Box(dx, dy) = 0;
                        if length(Zone_Fast{dx, dy}) > 1
                            for ic = 1: (length(Zone_Fast{dx, dy}) - 1)
                                Link_Fast(num_block, Zone_Fast{dx, dy}(ic)) = 1;
                                Link_Fast(Zone_Fast{dx, dy}(ic) ,num_block) = 1;
                            end
                        end
                        Ncoved_Fast(dx, dy) = length(Zone_Fast{dx, dy});
                    end
                end
            end
        end
    end
    num_block = num_block + 1;
end
num_block = num_block - 1;

Dis_Fast = inf * ones(size(Link_Fast));
for i = 1: num_block
    for j = i: num_block
        dd = ((Centre_Fast(i, 1) - Centre_Fast(j, 1))^2 + (Centre_Fast(i, 2) - Centre_Fast(j, 2))^2)^0.5;
        Dis_Fast(i, j) = dd;
        Dis_Fast(j, i) = dd;
    end
end
end