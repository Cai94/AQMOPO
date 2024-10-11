function [Pop, F] = NonDominatedSorting(Pop)
n_pop = numel(Pop);
for i = 1 : n_pop
    Pop(i).DominationSet = [];
    Pop(i).DominatedCount = 0;
end

F{1} = []; 
for i = 1: n_pop
    for j = i + 1 : n_pop
        p = Pop(i);
        q = Pop(j);
        if Dominates(p, q)
            p.DominationSet = [p.DominationSet, j];
            q.DominatedCount = q.DominatedCount + 1;
        end
        if Dominates(q, p)
            q.DominationSet = [q.DominationSet, i];
            p.DominatedCount = p.DominatedCount + 1;
        end
        Pop(i) = p;
        Pop(j) = q;
    end
    if Pop(i).DominatedCount == 0
        F{1} = [F{1}, i];
        Pop(i).Rank = 1;
    end
end

k = 1;
while true  
    Q = [];
    for i = F{k}
        p = Pop(i);
        for j = p.DominationSet
            q = Pop(j);
            q.DominatedCount = q.DominatedCount - 1;
            if q.DominatedCount == 0
                Q = [Q, j];
                q.Rank = k+1;
            end
            Pop(j) = q;
        end
    end
    if isempty(Q)
        break;
    end
    F{k+1} = Q;
    k = k + 1;
end
end