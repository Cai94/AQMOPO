function P = sumProbQ(Block, Centre, Table, n_bit)
n_c = size(Centre, 1);
r = Block.R_cov;
P = zeros(n_c, n_bit);
for i = 1: n_c
    cx = Centre(i, 1);
    cy = Centre(i, 2);
    P_box = zeros(1, n_bit);
    n_sum = 0;
    for j = 1: size(Table, 1)
        tx = Table{j, 1}(1);
        ty = Table{j, 1}(2);

        dL = ((cx - tx)^2 + (cy - ty)^2)^0.5;
        theta = 2 * acos((dL/2) / r);
        area_sector = theta * r^2 /2;
        area_tri = dL^2 * tan(theta/2) /4;
        area_intersect = 2 * (area_sector - area_tri);
                
        if area_intersect > 0
            Q_pop = Table{j, 2};
            P_pop = Q_pop.^2;
            P_box = P_box + P_pop * area_intersect / (2*pi * r)^2;
            n_sum = n_sum + 1;
        end
    end
    
    if n_sum == 0
        P_box = (1/n_bit) * ones(1, n_bit);
    else
        P_box = P_box ./ sum(P_box);
    end
    P(i, :) = P_box;
end
end