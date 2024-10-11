function newpop = selectPopulation(Pop, F, d, rho, Params)
n_pop = Params.nPop;
if numel(Pop) == n_pop
    return;
end

newpop = [];
for l = 1: numel(F)
    if numel(newpop) + numel(F{l}) >= n_pop 
        LastFront = F{l};
        break;
    end
    newpop = [newpop, Pop(F{l})]; 
end

while true
    [~, j] = min(rho);
    AssocitedFromLastFront = [];
    for i = LastFront
        if Pop(i).AssociatedRef == j
            AssocitedFromLastFront = [AssocitedFromLastFront, i];
        end
    end
    if isempty(AssocitedFromLastFront)
        rho(j) = inf;
        continue;
    end
    if rho(j) == 0
        ddj = d(AssocitedFromLastFront, j);
        [~, new_member_ind] = min(ddj);
    else
        new_member_ind = randi(numel(AssocitedFromLastFront));
    end
    MemberToAdd = AssocitedFromLastFront(new_member_ind);
    LastFront(LastFront == MemberToAdd) = [];
    newpop = [newpop, Pop(MemberToAdd)];
    rho(j) = rho(j) + 1;
    if numel(newpop) >= n_pop
        break;
    end
end
end
