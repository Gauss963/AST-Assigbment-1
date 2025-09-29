function mat = ReadMatrix(IREAD, row, col)
% for INPUT function
    mat = zeros(row,col);
    for j = 1:col
        line = fgets(IREAD);
        num = str2num(line);
        mat(:,j) = num(2:row+1);
    end
end