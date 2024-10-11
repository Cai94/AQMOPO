function CovMap = calCovMap(Chromo, Block, Map_Fast)
n_UUV = size(Chromo.Task ,1);
ModCov = Block.ModCov;
CovMap = {};
for i = 1: n_UUV
    CovMap{i} = ones(size(Block.Zone));
    Index_dis = [Chromo.Task{i}];
    Box = zeros(size(Block.Zone));
    for j = 1: length(Index_dis)
        n_select = Index_dis(j);
        Box = Box + Map_Fast{n_select};
    end
    Box = (Box > 0);
    CovMap{i} = 1 - Box;
end
end