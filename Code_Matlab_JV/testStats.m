%% Table 1

df = readtable("C:/Users/JosephVermeil/Desktop/ActinCortexAnalysis/DataAnalysis/TestDataSet/testStats01_avgPerCell.csv",'DatetimeType','text','TextType','string','Delimiter',',');
df_1 = df(df.substrate_Drug == "20um fibronectin discs & doxycyclin",:);
[~,b]=unique(df.substrate_Drug,'stable');
categories=df.substrate_Drug(b,:);
Ncat = length(categories);
for i=1:Ncat
    for j=i:Ncat
        if i ~= j
            fprintf(categories(i));
            fprintf(' vs ');
            fprintf(categories(j));
            fprintf('\n');
            % df(df.substrate_Drug == categories(i),:).EChadwick
            % fprintf('RS: ');
            % fprintf(RS);
            [~,TT] = ttest2(df(df.substrate_Drug == categories(i),:).EChadwick, df(df.substrate_Drug == categories(j),:).EChadwick)
            RS = ranksum(df(df.substrate_Drug == categories(i),:).EChadwick, df(df.substrate_Drug == categories(j),:).EChadwick)
            fprintf('\n');
            fprintf('\n');
        end
    end
end

%% Table 2

%% Table 3

df = readtable("C:/Users/JosephVermeil/Desktop/ActinCortexAnalysis/DataAnalysis/TestDataSet/testStats03_FakeData.csv",'DatetimeType','text','TextType','string','Delimiter',',');

categories=string(df.Properties.VariableNames);
Ncat = length(categories);
for i=3:Ncat
%     char(A)
    fprintf(categories(2));
    fprintf(' vs ');
    fprintf(categories(i));
    fprintf('\n');
    [~,TT] = ttest2(df(:,categories(2)).Variables, df(:,categories(i)).Variables)
    RS = ranksum(df(:,categories(2)).Variables, df(:,categories(i)).Variables)
    fprintf('\n');
    fprintf('\n');

end

%% Table 4

fprintf('Large fake data');
fprintf('\n');
fprintf('\n');

df = readtable("C:/Users/JosephVermeil/Desktop/ActinCortexAnalysis/DataAnalysis/TestDataSet/testStats04_FakeDataLarge.csv",'DatetimeType','text','TextType','string','Delimiter',',');

categories=string(df.Properties.VariableNames);
Ncat = length(categories);
for i=3:Ncat
%     char(A)
    fprintf(categories(2));
    fprintf(' vs ');
    fprintf(categories(i));
    fprintf('\n');
    [~,TT] = ttest2(df(:,categories(2)).Variables, df(:,categories(i)).Variables)
    RS = ranksum(df(:,categories(2)).Variables, df(:,categories(i)).Variables)
    fprintf('\n');
    fprintf('\n');

end