function Zr = GenerateReferencePoints(M, p)
    Zr = GetFixedRowSumIntegerMatrix(M, p)' / p;
end
function A = GetFixedRowSumIntegerMatrix(M, RowSum)
    if M < 1
        error('error');
    end    
    if floor(M) ~= M
        error('error');
    end
    if M == 1
        A = RowSum;
        return;
    end
    A = [];
    for i = 0: RowSum
        B = GetFixedRowSumIntegerMatrix(M - 1, RowSum - i);
        A = [A; i * ones(size(B, 1), 1), B];
    end
end
