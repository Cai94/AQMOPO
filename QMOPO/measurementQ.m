function [Task, Cov2Task] = measurementQ(Q, Centre)
[N_Q, N_bit] = size(Q);
Task = cell(N_bit ,1);
for i = 1: N_Q
    Box = Q(i, :);
    Box = Box.^2;
    dp = sum(Box);
    n2Bit = [];
    
    sumP = 0;
    randNum = rand * dp;
    for rn = 1: N_bit
        sumP = sumP + Box(rn);
        if sumP >= randNum
            n2Bit = rn;
            break;
        end
    end
    Task{n2Bit} = [Task{n2Bit}, i];
end

Cov2Task = [];
dn = 1;
for j = 1: N_bit
    for nc = 1: length(Task{j})
        num = Task{j}(nc);
        Cov2Task(dn, 1) = Centre(num, 1);
        Cov2Task(dn, 2) = Centre(num, 2);
        Cov2Task(dn, 3) = j;
        dn = dn + 1;
    end
end
end