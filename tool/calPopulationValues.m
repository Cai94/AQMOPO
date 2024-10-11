function Value = calPopulationValues(Chromo, Plan, Block, Map_Fast)
addpath(genpath('..\tool'));
n_UUV = size(Chromo.Task ,1);
n_cov = size(Chromo.Centre, 1);
Initial = Plan.Initial;
v_UUV = Plan.v_UUV;
p_task = Plan.p_task;
ModCov = Block.ModCov;
R_mod = length(ModCov) / 2;

dis_task = [];
for i = 1: n_UUV
    Index_dis = [Chromo.Task{i}];
    if length(Index_dis) > 1
        Box_Initial = Chromo.DisInitial(Index_dis);
        Box_Dis = Chromo.Dis(Index_dis, :);
        D = Box_Dis(:, Index_dis);
        dis_task(i) = disCal(D) + min(Box_Initial) * 2;
    elseif length(Index_dis) == 1
        dis_task(i) = Chromo.DisInitial(Index_dis) * 2;
    else
        dis_task(i) = 0;
    end
end
dis_all = 2 * (Block.nPoint_main / (pi * R_mod^2)) * length(Block.Zone);
value_dis = sum(dis_task) / dis_all;

for i = 1: n_UUV
    t_move = dis_task(i) ./ v_UUV;
    s_task = 0;
    Index_dis = [Chromo.Task{i}];
    Box = zeros(size(Block.Zone));
    for j = 1: length(Index_dis)
        n_select = Index_dis(j);
        Box = Box + Map_Fast{n_select};
    end
    Box = (Box > 0);
    s_task = sum(sum(Box));
    t_task = s_task ./ p_task;
    T_sum(i) = t_move + t_task;
end
value_equil = std(T_sum) ./ mean(T_sum);
value_equil = value_equil^0.1/2;

Dir_water = Chromo.DirWater;
xmin_mian = max(1, min(Block.Border_main(:, 1)));
ymin_mian = max(1, min(Block.Border_main(:, 2)));
xmax_mian = min(max(Block.Border_main(:, 1)), size(Block.Zone, 1));
ymax_mian = min(max(Block.Border_main(:, 2)), size(Block.Zone, 2));
NTurn = zeros(n_UUV, 1);
for i = 1: n_UUV
    if sign(Dir_water(i, 1)) == 0
        dd = 4 * Plan.R;
        dn = 0;
        x = xmin_mian + dd/2 + dn * dd;
        while x <= xmax_mian
            flag = 0;
            for y = ymin_mian: ymax_mian
                if flag == 0 && Chromo.CovMap{i}(x, y) == 0
                    flag = 1;
                    NTurn(i) = NTurn(i) + 1;
                elseif flag == 1 && Chromo.CovMap{i}(x, y) == 1
                    flag = 0;
                    NTurn(i) = NTurn(i) + 1;
                elseif flag == 1 && Chromo.CovMap{i}(x, y) == 0 && y == ymax_mian
                    flag = 0;
                    NTurn(i) = NTurn(i) + 1;
                end
            end
            dn = dn + 1;
            x = xmin_mian + dd/2 + dn * dd;
        end
    elseif sign(Dir_water(i, 2)) == 0 && sign(Dir_water(i, 1)) ~= 0
        dd = 4 * Plan.R;
        dn = 0;
        y = ymin_mian + dd/2 + dn * dd;
        while y <= ymax_mian
            flag = 0;
            for x = xmin_mian: xmax_mian
                if flag == 0 && Chromo.CovMap{i}(x, y) == 0
                    flag = 1;
                    NTurn(i) = NTurn(i) + 1;
                elseif flag == 1 && Chromo.CovMap{i}(x, y) == 1
                    flag = 0;
                    NTurn(i) = NTurn(i) + 1;
                elseif flag == 1 && Chromo.CovMap{i}(x, y) == 0 && x == xmax_mian
                    flag = 0;
                    NTurn(i) = NTurn(i) + 1;
                end
            end
            dn = dn + 1;
            y = ymin_mian + dd/2 + dn * dd;
        end
    elseif sign(Dir_water(i, 1)) == sign(Dir_water(i, 2))
        dd = 4 * Plan.R * ((Dir_water(i, 1))^2 + (Dir_water(i, 2))^2)^0.5 / abs(Dir_water(i, 2));
        ddx = abs(Dir_water(i, 1)) / ((Dir_water(i, 1))^2 + (Dir_water(i, 2))^2)^0.5;
        ddy = abs(Dir_water(i, 2)) / ((Dir_water(i, 1))^2 + (Dir_water(i, 2))^2)^0.5;
        dn = 0;
        xx = xmin_mian - ymax_mian * abs(Dir_water(i, 1)) / abs(Dir_water(i, 2)) + dd/2 + dn * dd;
        while xx <= xmax_mian - ymin_mian * abs(Dir_water(i, 1)) / abs(Dir_water(i, 2))
            flag = 0;
            ddn = 0;
            x = xx + ddn * ddx;
            y = ddn * ddy;
            x = round(x);
            y = round(y);
            while y <= ymax_mian
                if x >= xmin_mian && x <= xmax_mian && y >= ymin_mian
                    if flag == 0 && Chromo.CovMap{i}(x, y) == 0
                        flag = 1;
                        NTurn(i) = NTurn(i) + 1;
                    elseif flag == 1 && Chromo.CovMap{i}(x, y) == 1
                        flag = 0;
                        NTurn(i) = NTurn(i) + 1;
                    elseif flag == 1 && Chromo.CovMap{i}(x, y) == 0 && y >= ymax_mian
                        flag = 0;
                        NTurn(i) = NTurn(i) + 1;
                    end
                elseif flag == 1 && x >= xmax_mian
                    flag = 0;
                    NTurn(i) = NTurn(i) + 1;
                end
                ddn = ddn + 1;
                x = xx + ddn * ddx;
                y = ddn * ddy;
                x = round(x);
                y = round(y);
            end
            dn = dn + 1;
            xx = xmin_mian - ymax_mian * abs(Dir_water(i, 1)) / abs(Dir_water(i, 2)) + dd/2 + dn * dd;
        end
    elseif sign(Dir_water(i, 1)) ~= sign(Dir_water(i, 2))
        dd = 4 * Plan.R * ((Dir_water(i, 1))^2 + (Dir_water(i, 2))^2)^0.5 / abs(Dir_water(i, 2));
        ddx = abs(Dir_water(i, 1)) / ((Dir_water(i, 1))^2 + (Dir_water(i, 2))^2)^0.5;
        ddy = abs(Dir_water(i, 2)) / ((Dir_water(i, 1))^2 + (Dir_water(i, 2))^2)^0.5;
        dn = 0;
        xx = xmin_mian + ymin_mian * abs(Dir_water(i, 1)) / abs(Dir_water(i, 2)) + dd/2 + dn * dd;
        while xx <= xmax_mian + ymax_mian * abs(Dir_water(i, 1)) / abs(Dir_water(i, 2))
            flag = 0;
            ddn = 0;
            x = xx - ddn * ddx;
            y = ddn * ddy;
            x = round(x);
            y = round(y);
            while y <= ymax_mian
                if x >= xmin_mian && x <= xmax_mian && y >= ymin_mian
                    if flag == 0 && Chromo.CovMap{i}(x, y) == 0
                        flag = 1;
                        NTurn(i) = NTurn(i) + 1;
                    elseif flag == 1 && Chromo.CovMap{i}(x, y) == 1
                        flag = 0;
                        NTurn(i) = NTurn(i) + 1;
                    elseif flag == 1 && Chromo.CovMap{i}(x, y) == 0 && y >= ymax_mian
                        flag = 0;
                        NTurn(i) = NTurn(i) + 1;
                    end
                elseif flag == 1 && x <= xmin_mian
                    flag = 0;
                    NTurn(i) = NTurn(i) + 1;
                end
                ddn = ddn + 1;
                x = xx - ddn * ddx;
                y = ddn * ddy;
                x = round(x);
                y = round(y);
            end
            dn = dn + 1;
            xx = xmin_mian + ymin_mian * abs(Dir_water(i, 1)) / abs(Dir_water(i, 2)) + dd/2 + dn * dd;
        end
    end
end
nTurn_all = (Block.nPoint_main / (pi * (R_mod /2)^2)) * R_mod / Plan.R;
nTurn_all = round(nTurn_all);
value_rule = sum(NTurn) / nTurn_all;

Value = [value_dis; value_equil; value_rule];
end

function dis = disCal(D)
n_cov = length(D);
gn = n_cov - 1;

path = [1];
dis = 0;
for i = 1: gn
    Index_Un = [1: n_cov];
    Index_Un(path) = [];
    
    Box = D(path(i), Index_Un);
    [dis_min, index_min] = min(Box);
    
    path = [path, Index_Un(index_min)];
    dis = dis + dis_min;
end
end