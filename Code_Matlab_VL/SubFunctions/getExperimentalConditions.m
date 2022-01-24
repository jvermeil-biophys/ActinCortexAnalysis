function currentExperimentalConditions = getExperimentalConditions(table, date, manip)
    isExperimentFound = any((table.date == date & table.manip == manip));
    if isExperimentFound
        currentExperimentalConditions = table((table.date == date & table.manip == manip),:);
    else
        currentExperimentalConditions = table((table.date == "DEFAULT" & table.manip == "DEFAULT"),:);
    end
    currentExperimentalConditions = table2struct(currentExperimentalConditions);

end