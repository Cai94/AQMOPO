function Table = makeTableQ(Pop)
Table = {};
dn = 1;
for i = 1: length(Pop)
    for j = 1: size(Pop(i).Task, 1)
        for nc = 1: length(Pop(i).Task{j})
            num = Pop(i).Task{j}(nc);
            Table{dn, 1} = Pop(i).Centre(num, :);
            Table{dn, 2} = Pop(i).TaskQ(num, :);
            dn = dn + 1;
        end
    end
end
end