function Record = AQMOPO(Plan, Block, Params)
addpath(genpath('.\tool'));
addpath(genpath('.\QMOPO'));

n_UUV = Plan.n_UUV;
R = Plan.R;
Initial = Plan.Initial;

Zone = Block.Zone;
Limit = Block.Limit;
Border_main = Block.Border_main;
Zone_main = Block.Zone_main;
ModCov = Block.ModCov;

n_pop = Params.n_pop;
gn_max = Params.gn_max;
n_obj = Params.n_obj;
Params.nPop = n_pop;
Params.Prob = cell(n_UUV, 1);
n_pop = Params.n_pop;
gn_max = Params.gn_max;
n_obj = Params.n_obj;
n_division = Params.n_division;
Params.nPop = n_pop;
Point_refer = GenerateReferencePoints(n_obj, n_division);
Params.Zr = Point_refer;
Params.nZr = size(Point_refer, 2);
Params.TaskTableF1 = [];
Params.TaskTableAll = [];
Params.zmin = [];
Params.zmax = [];
Params.smin = [];

Individual.Ncoved = [];
Individual.Centre = [];
Individual.Link = [];
Individual.Dis = [];
Individual.DisInitial = [];
Individual.Task= {};
Individual.Cov2Task = [];
Individual.DirWater = [];
Individual.CovMap = {};
Individual.Cost = [];
Individual.NormalizedCost = [];
Individual.DominationSet = [];
Individual.DominatedCount = [];
Individual.Rank = [];
Individual.AssociatedRef = [];
Individual.DistanceToAssociatedRef = [];
Individual.Index_old = [];

Record.GN = [];
Record.F1 = {};
Record.BestCost = [];
Record.Time = [];

Pop = repmat(Individual, 1);
for i = 1: n_pop
    [Map_Fast, Ncoved_Fast, Centre_Fast, Link_Fast, Dis_Fast] = covFast(Zone_main, ModCov);
    Pop(i).Index_old = i;
    Pop(i).Ncoved = Ncoved_Fast; 
    Pop(i).Centre = Centre_Fast;
    Pop(i).Link = Link_Fast;
    Pop(i).Dis = Dis_Fast;
    Pop(i).DisInitial = disInitial(Initial, Centre_Fast);
    Pop(i).TaskQ = (1/sqrt(n_UUV)) * ones(size(Centre_Fast, 1), n_UUV);
    [Pop(i).Task, Pop(i).Cov2Task] = measurementQ(Pop(i).TaskQ, Pop(i).Centre);
    Pop(i).DirWater = calWaterDirection(Pop(i), Block);
    Pop(i).CovMap = calCovMap(Pop(i), Block, Map_Fast);
    Pop(i).Cost = calPopulationValues(Pop(i), Plan, Block, Map_Fast);
end
[Pop, Params] = NormalizePopulation(Pop, Params);
[Pop, F] = NonDominatedSorting(Pop);
F1 = Pop(F{1});
CostBest_0 = min([F1.Cost],[], 2);
St_Q = (1/sqrt(4)) * ones(n_pop, 4);

for gn = 1: gn_max
    TableQ_All = makeTableQ(Pop);
    TableQ_F1 = makeTableQ(F1);
    Record_St = [];
    
    Pop_new = repmat(Individual, 1);
    [Map_Fast, Ncoved_Fast, Centre_Fast, Link_Fast, Dis_Fast] = covFast(Zone_main, ModCov);
    for i = 1: n_pop
        Pop_new(i).Ncoved = Ncoved_Fast;
        Pop_new(i).Centre = Centre_Fast;
        Pop_new(i).Link = Link_Fast;
        Pop_new(i).Dis = Dis_Fast;
        Pop_new(i).DisInitial = disInitial(Initial, Centre_Fast);
        TableQ_i = makeTableQ(Pop(i));
        
        St_prob = St_Q.^2;
        sumSt = 0;
        randSt = rand * sum(St_prob(i, :));
        for nm = 1: 4
            sumSt = sumSt + St_prob(i, nm);
            if sumSt >= randSt
                St = nm;
                break;
            end
        end
        Record_St(i, 1) = St;
        
        if St == 1
            for nc = 1: size(Pop_new(i).Centre, 1)
                P_i = sumProbQ(Block, Pop_new(i).Centre(nc, :), TableQ_i, n_UUV);
                P_F1 = sumProbQ(Block, Pop_new(i).Centre(nc, :), TableQ_F1, n_UUV);
                P_All = sumProbQ(Block, Pop_new(i).Centre(nc, :), TableQ_All, n_UUV);
               
                P_selest = P_i + P_F1 .* levy(1, n_UUV, 1.8) + rand(1, n_UUV) .* P_All * (1 - gn/gn_max)^(2 * gn/gn_max);
                P_selest = P_selest .* (P_selest >= 0);
                if sum(P_selest) == 0 
                    P_selest = ones(1, n_UUV);
                end
                P_selest = P_selest ./ sum(P_selest);
                Pop_new(i).TaskQ(nc, :) = sqrt(P_selest);
            end
            [Pop_new(i).Task, Pop_new(i).Cov2Task] = measurementQ(Pop_new(i).TaskQ, Pop_new(i).Centre);
            Pop_new(i).DirWater = calWaterDirection(Pop_new(i), Block);
            Pop_new(i).CovMap = calCovMap(Pop_new(i), Block, Map_Fast);
            Pop_new(i).Cost = calPopulationValues(Pop_new(i), Plan, Block, Map_Fast);
        elseif St == 2
            for nc = 1: size(Pop_new(i).Centre, 1)
                P_i = sumProbQ(Block, Pop_new(i).Centre(nc, :), TableQ_i, n_UUV);
                P_F1 = sumProbQ(Block, Pop_new(i).Centre(nc, :), TableQ_F1, n_UUV);
                P_All = sumProbQ(Block, Pop_new(i).Centre(nc, :), TableQ_All, n_UUV);
                
                P_selest = P_i + P_F1 .* levy(1, n_UUV, 1.8) + rand(1, n_UUV) .* ones(1, n_UUV);
                P_selest = P_selest .* (P_selest >= 0);
                if sum(P_selest) == 0
                    P_selest = ones(1, n_UUV);
                end
                P_selest = P_selest ./ sum(P_selest);
                Pop_new(i).TaskQ(nc, :) = sqrt(P_selest);
            end
            [Pop_new(i).Task, Pop_new(i).Cov2Task] = measurementQ(Pop_new(i).TaskQ, Pop_new(i).Centre);
            Pop_new(i).DirWater = calWaterDirection(Pop_new(i), Block);
            Pop_new(i).CovMap = calCovMap(Pop_new(i), Block, Map_Fast);
            Pop_new(i).Cost = calPopulationValues(Pop_new(i), Plan, Block, Map_Fast);
        elseif St == 3
            alpha = rand(1, n_UUV) / 5;
            H = rand(1);
            if H < 0.5
                for nc = 1: size(Pop_new(i).Centre, 1)
                    P_i = sumProbQ(Block, Pop_new(i).Centre(nc, :), TableQ_i, n_UUV);
                    P_F1 = sumProbQ(Block, Pop_new(i).Centre(nc, :), TableQ_F1, n_UUV);
                    P_All = sumProbQ(Block, Pop_new(i).Centre(nc, :), TableQ_All, n_UUV);
                    
                    P_selest = alpha .* (1 - gn / gn_max) .* (P_i - P_All);
                    P_selest = P_selest .* (P_selest >= 0);
                    if sum(P_selest) == 0
                        P_selest = ones(1, n_UUV);
                    end
                    P_selest = P_selest ./ sum(P_selest);
                    Pop_new(i).TaskQ(nc, :) = sqrt(P_selest);
                end
                [Pop_new(i).Task, Pop_new(i).Cov2Task] = measurementQ(Pop_new(i).TaskQ, Pop_new(i).Centre);
                Pop_new(i).DirWater = calWaterDirection(Pop_new(i), Block);
                Pop_new(i).CovMap = calCovMap(Pop_new(i), Block, Map_Fast);
                Pop_new(i).Cost = calPopulationValues(Pop_new(i), Plan, Block, Map_Fast);
            else
                for nc = 1: size(Pop_new(i).Centre, 1)
                    P_i = sumProbQ(Block, Pop_new(i).Centre(nc, :), TableQ_i, n_UUV);
                    P_F1 = sumProbQ(Block, Pop_new(i).Centre(nc, :), TableQ_F1, n_UUV);
                    P_All = sumProbQ(Block, Pop_new(i).Centre(nc, :), TableQ_All, n_UUV);
                    
                    P_selest = alpha .* exp(-gn ./ (rand(1, n_UUV) .* gn_max)) .* ones(1, n_UUV);
                    P_selest = P_selest .* (P_selest >= 0);
                    if sum(P_selest) == 0
                        P_selest = ones(1, n_UUV);
                    end
                    P_selest = P_selest ./ sum(P_selest);
                    Pop_new(i).TaskQ(nc, :) = sqrt(P_selest);
                end
                [Pop_new(i).Task, Pop_new(i).Cov2Task] = measurementQ(Pop_new(i).TaskQ, Pop_new(i).Centre);
                Pop_new(i).DirWater = calWaterDirection(Pop_new(i), Block);
                Pop_new(i).CovMap = calCovMap(Pop_new(i), Block, Map_Fast);
                Pop_new(i).Cost = calPopulationValues(Pop_new(i), Plan, Block, Map_Fast);
            end
        else
            sita = rand(1) * pi;
            for nc = 1: size(Pop_new(i).Centre, 1)
                P_i = sumProbQ(Block, Pop_new(i).Centre(nc, :), TableQ_i, n_UUV);
                P_F1 = sumProbQ(Block, Pop_new(i).Centre(nc, :), TableQ_F1, n_UUV);
                P_All = sumProbQ(Block, Pop_new(i).Centre(nc, :), TableQ_All, n_UUV);
                
                P_selest = P_i + rand(1, n_UUV) .* cos((pi * gn )/ (2 * gn_max)) .* (P_F1 - P_i) - cos(sita) * (gn / gn_max)^(2 / gn_max) * (P_i - P_F1);
                P_selest = P_selest .* (P_selest >= 0);
                if sum(P_selest) == 0
                    P_selest = ones(1, n_UUV);
                end
                P_selest = P_selest ./ sum(P_selest);
                Pop_new(i).TaskQ(nc, :) = sqrt(P_selest);
            end
            [Pop_new(i).Task, Pop_new(i).Cov2Task] = measurementQ(Pop_new(i).TaskQ, Pop_new(i).Centre);
            Pop_new(i).DirWater = calWaterDirection(Pop_new(i), Block);
            Pop_new(i).CovMap = calCovMap(Pop_new(i), Block, Map_Fast);
            Pop_new(i).Cost = calPopulationValues(Pop_new(i), Plan, Block, Map_Fast);
        end
        Pop_new(i).Index_old = i;
    end
    
    Record.St_P{gn} = St_prob;
    for i = 1: n_pop
        num_W = Record_St(i);
        Index_other = find([1:4]~=num_W);
        Box = ones(1, 4);
        if Dominates(Pop_new(i), Pop(i))
            if St_prob(i, num_W) < 1/2
                dprob = 1/2 - St_prob(i, num_W);
                dprob_up = dprob * (gn/gn_max)^(0.8 * (gn_max - gn)/gn_max) * rand;
                Box(num_W) = St_prob(i, num_W) + dprob_up;
                dprob_other = dprob_up / 3;
                Box(Index_other) = St_prob(i, Index_other) - dprob_other;
            else
                dprob_up = 1/2 - St_prob(i, num_W);
                Box(num_W) = St_prob(i, num_W) + dprob_up;
                dprob_other = dprob_up / 3;
                Box(Index_other) = St_prob(i, Index_other) - dprob_other;
            end
        else
            if St_prob(i, num_W) > 1/6
                dprob = St_prob(i, num_W) - 1/6;
                dprob_up = dprob * (gn/gn_max)^(0.8 * (gn_max - gn)/gn_max) * rand;
                Box(num_W) = St_prob(i, num_W) - dprob_up;
                dprob_other = dprob_up / 3;
                Box(Index_other) = St_prob(i, Index_other) + dprob_other;
            else
                dprob_up = St_prob(i, num_W) - 1/6;
                Box(num_W) = St_prob(i, num_W) - dprob_up;
                dprob_other = dprob_up / 3;
                Box(Index_other) = St_prob(i, Index_other) + dprob_other;
            end
        end
        Box = Box ./ sum(Box);
        St_prob(i, :) = Box;
    end
    St_Q = St_prob.^0.5;
    
    Pop_box = [Pop_new, F1];
    [Pop_box, Params] = NormalizePopulation(Pop_box, Params);
    [Pop_box, F] = NonDominatedSorting(Pop_box);
    [Pop_box, d, rho] = AssociateToReferencePoint(Pop_box, Params);
    Pop_selected = selectPopulation(Pop_box, F, d, rho, Params);
    
    Pop = Pop_selected;
    [Pop, Params] = NormalizePopulation(Pop, Params);
    [Pop, F] = NonDominatedSorting(Pop);
    F1 = Pop(F{1});
    
    Box = St_prob;
    for i = 1: n_pop
        num_old = Pop(i).Index_old;
        St_prob(i, :) = Box(num_old, :);
        Pop(i).Index_old = i;
    end
    
    Record.GN(gn) = gn;
    Record.F1{gn} = {F1.Cov2Task};
    CostBest_F1 = min([F1.Cost],[], 2);
    n_obj = size(CostBest_F1 ,1);
    for nobj = 1: n_obj
        if gn == 1
            bestCost = CostBest_0(nobj);
        else
            bestCost = Record.BestCost(nobj, 1: (gn - 1));
        end
        Record.BestCost(nobj, gn) = min([CostBest_F1(nobj), bestCost]);
    end
    
    fprintf('Generation = %d\n', gn);
    fprintf('Best Cost = ');
    for nobj = 1: n_obj
        fprintf('%f  ', Record.BestCost(nobj, gn));
    end
    fprintf('\n');

end
end