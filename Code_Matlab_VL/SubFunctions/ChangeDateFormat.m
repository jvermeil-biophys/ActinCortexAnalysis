%% Change Dates into BigTables

load('D:\Matlab Analysis\Data_Joseph\MatFiles\MecaDataTable.mat')
N = size(BigTable);
ColsOfInterest = ["ExpDay" "CellName"];
for iCol = 1:length(ColsOfInterest)
    currentCol = BigTable.(ColsOfInterest(iCol));
    for i = 1:N
        currentString = char(currentCol(i));
        if strlength(currentString) == 8
            newString = [currentString(7:8) currentString(3:6) currentString(1:2)];
        elseif strlength(currentString) >= 8
            newString = [currentString(7:8) currentString(3:6) currentString(1:2) currentString(9:end)];
        else
            error
        end
        currentCol(i) = string(newString);
    end
    BigTable.(ColsOfInterest(iCol)) = currentCol;
end

% save('D:\Matlab Analysis\Data_Joseph\MatFiles\MecaDataTable.mat','BigTable');

load('D:\Matlab Analysis\Data_Joseph\MatFiles\MecaData3T3Table.mat')
N = size(SubsetTable);
ColsOfInterest = ["ExpDay" "CellName"];
for iCol = 1:length(ColsOfInterest)
    currentCol = BigTable.(ColsOfInterest(iCol));
    for i = 1:N
        currentString = char(currentCol(i));
        if strlength(currentString) == 8
            newString = [currentString(7:8) currentString(3:6) currentString(1:2)];
        elseif strlength(currentString) >= 8
            newString = [currentString(7:8) currentString(3:6) currentString(1:2) currentString(9:end)];
        else
            error
        end
        currentCol(i) = string(newString);
    end
    BigTable.(ColsOfInterest(iCol)) = currentCol;
end

% save('D:\Matlab Analysis\Data_Joseph\MatFiles\MecaData3T3Table.mat','SubsetTable');

%% Change Dates into R2Vs

targetDirectory = 'D:\Matlab Analysis\Data_Joseph\MatFiles\R2V';
listFiles = dir(targetDirectory);
listFiles = listFiles(3:end);
[nFiles, ~] = size(listFiles);
FieldsOfInterest = ["exp" "stack"];
for iFile = 1:nFiles
    fileName = listFiles(iFile).name;
    load([targetDirectory filesep fileName])
    nCells = length(MT);
    for iCell = 1:nCells
        for iField = 1:length(FieldsOfInterest)
            currentString = char(MT{1,iCell}.(FieldsOfInterest(iField)));
            newString = currentString;
%             if strlength(currentString) == 8
%                 newString = [currentString(7:8) currentString(3:6) currentString(1:2)];
%             elseif strlength(currentString) >= 8
%                 newString = [currentString(7:8) currentString(3:6) currentString(1:2) currentString(9:end)];
%             else
%                 error
%             end
            MT{1,iCell}.(FieldsOfInterest(iField)) = char(newString);
        end
    end
    savePath = [targetDirectory filesep fileName];
%     save([targetDirectory filesep fileName],'MT','ListD','ListF')
end


%% Change Dates into V2Ds

targetDirectory = 'D:\Matlab Analysis\Data_Joseph\MatFiles\V2D';
listFiles = dir(targetDirectory);
listFiles = listFiles(3:end);
[nFiles, ~] = size(listFiles);
FieldsOfInterest = ["name" "stack"];
for iFile = 1:nFiles
    fileName = listFiles(iFile).name;
    if contains(fileName, 'V2D')
        load([targetDirectory filesep fileName])
        nCells = length(MR);
        for iCell = 1:nCells
            if ~isempty(MR{1,iCell})
                for iField = 1:length(FieldsOfInterest)
                    if isfield(MR{1,1},FieldsOfInterest(iField))
                        currentString = char(MR{1,iCell}.(FieldsOfInterest(iField)));
                        newString = currentString;
%                         if strlength(currentString) == 8
%                             newString = [currentString(7:8) currentString(3:6) currentString(1:2)];
%                         elseif strlength(currentString) >= 8
%                             newString = [currentString(7:8) currentString(3:6) currentString(1:2) currentString(9:end)];
%                         else
%                             error
%                         end
%                         if ~strcmp(newString(1:2),'21') && ~strcmp(newString(1:2),'20')
%                             cprintf('shit\n');
%                         else
%                             cprintf('text', [newString ' (' char(FieldsOfInterest(iField)) ')\n']);
%                         end
                        MR{1,iCell}.(FieldsOfInterest(iField)) = char(newString);
                    end
                end
            end
        end
        savePath = [targetDirectory filesep fileName];
        test1 = exist('AllRamp', 'var');
        if exist('AllRamp', 'var')
%             save([targetDirectory filesep fileName],'MR','AllRamp')
            clear('AllRamp', 'MR')
        else
%             save([targetDirectory filesep fileName],'MR')
            clear('AllRamp', 'MR')
        end
    end
end