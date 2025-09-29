function HEADLINE(ID,IREAD)
    while ~feof(IREAD)
        temp = fgets(IREAD);
        if ~isempty(temp) && temp(1) == ID
            return
        end
    end
end