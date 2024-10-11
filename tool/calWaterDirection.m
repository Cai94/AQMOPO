function Dir_water = calWaterDirection(Chromo, Block)
n_UUV = size(Chromo.Task ,1);
ModCov = Block.ModCov;
R_mod = length(ModCov) / 2;

V_water = cell(n_UUV, 1);
for iw = 1: size(Block.Water.coordinate, 1)
    for jw = 1: size(Block.Water.coordinate, 2)
        wx = Block.Water.coordinate{iw, jw}(1);
        wy = Block.Water.coordinate{iw, jw}(2);
        if Block.Zone_main(wx, wy) == 0
            for i = 1: n_UUV
                Index_dis = [Chromo.Task{i}];
                for j = 1: length(Index_dis)
                    n_select = Index_dis(j);
                    sx = Chromo.Centre(n_select, 1);
                    sy = Chromo.Centre(n_select, 2);
                    if (wx - sx)^2 + (wy - sy)^2 <= R_mod^2
                        V_water{i} = [[V_water{i}]; [Block.Water.V{iw, jw}(1), Block.Water.V{iw, jw}(2)]];
                        break;
                    end
                end
            end
        end
    end
end

V_mean = ones(n_UUV, 2);
for i = 1: n_UUV
    if ~isempty(V_water{i})
        V_mean(i, 1) = sum(rmoutliers([V_water{i}(:, 1)]));
        V_mean(i, 2) = sum(rmoutliers([V_water{i}(:, 2)]));
        V_mean(i, :) = V_mean(i, :)./size(V_water{i}, 1);
    else
        V_mean(i, :) = [0, 0];
    end
end

Dir_water = V_mean;
for i = 1: n_UUV
    Dir_water(i, :) = Dir_water(i, :)./(Dir_water(i, 1)^2 + Dir_water(i, 2)^2)^0.5;
end
end