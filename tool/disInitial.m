function D = disInitial(Initial, Centre_Fast)
n = size(Centre_Fast, 1);
D = [];
for i = 1: n
    D(i) = ((Initial(1) - Centre_Fast(i, 1))^2 + (Initial(2) - Centre_Fast(i, 2))^2)^0.5;
end
end