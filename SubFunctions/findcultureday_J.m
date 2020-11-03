function CultureTag = findcultureday_J(date)

DateList = {'dd-mm-yy'};

Ind = find(strcmp(DateList,date)==1);

if isempty(Ind)
    CultureTag = 'NoRef';
else
    
CultureList = {'Example'};

CultureTag = CultureList{Ind};
end

end