function tableExperimentalConditions = loadTableExperimentalConditions(path)
%     tableExperimentalConditions = readtable(path,'DatetimeType','text','TextType','string','Delimiter',';','HeaderLines',0);
    tableExperimentalConditions = readtable(path,'DatetimeType','text','TextType','string','Format','auto','Delimiter',';','HeaderLines',0);
    if isa(tableExperimentalConditions.magneticFieldCorrection,'string')
        tableExperimentalConditions.magneticFieldCorrection = str2double(replace(tableExperimentalConditions.magneticFieldCorrection,',','.'));
    end
    if isa(tableExperimentalConditions.scalePixelPerUm,'string')
        tableExperimentalConditions.scalePixelPerUm = str2double(replace(tableExperimentalConditions.scalePixelPerUm,',','.'));
    end
    if isa(tableExperimentalConditions.opticalIndexCorrection,'string')
        numbers = split(tableExperimentalConditions.opticalIndexCorrection, '/');
        if isa(numbers,'string')
            numbers = replace(numbers,',','.');
            numbers = str2double(numbers);
        end
        tableExperimentalConditions.opticalIndexCorrection = numbers(:,1) ./ numbers(:,2);
    end
end