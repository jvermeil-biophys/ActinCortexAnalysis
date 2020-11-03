function [VarName1,Area,Mean,StdDev,Min,Max,X,Y,XM,YM,Major,Minor,Angle1,Slice] = importfileOld(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [VARNAME1,AREA,MEAN,STDDEV,MIN,MAX,X,Y,XM,YM,MAJOR,MINOR,ANGLE1,SLICE]
%   = IMPORTFILE(FILENAME) Reads data from text file FILENAME for the
%   default selection.
%
%   [VARNAME1,AREA,MEAN,STDDEV,MIN,MAX,X,Y,XM,YM,MAJOR,MINOR,ANGLE1,SLICE]
%   = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows STARTROW
%   through ENDROW of text file FILENAME.
%
% Example:
%   [VarName1,Area,Mean,StdDev,Min,Max,X,Y,XM,YM,Major,Minor,Angle1,Slice] = importfile('K1_20n_Results.txt',2, 237);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2019/08/05 11:07:09

%% Initialize variables.
delimiter = '\t';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Format for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
%   column11: double (%f)
%	column12: double (%f)
%   column13: double (%f)
%	column14: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
VarName1 = dataArray{:, 1};
Area = dataArray{:, 2};
Mean = dataArray{:, 3};
StdDev = dataArray{:, 4};
Min = dataArray{:, 5};
Max = dataArray{:, 6};
X = dataArray{:, 7};
Y = dataArray{:, 8};
XM = dataArray{:, 9};
YM = dataArray{:, 10};
Major = dataArray{:, 11};
Minor = dataArray{:, 12};
Angle1 = dataArray{:, 13};
Slice = dataArray{:, 14};


