function Table = makeTable(Pop)
Table = [];
dn = 1;
for i = 1: length(Pop)
    for j = 1: size(Pop(i).Task, 1) % n_UUV
        for nc = 1: length(Pop(i).Task{j})
            num = Pop(i).Task{j}(nc);
            Table(dn, 1) = Pop(i).Centre(num, 1);
            Table(dn, 2) = Pop(i).Centre(num, 2);
            Table(dn, 3) = j;
            dn = dn + 1;
        end
    end
end
end