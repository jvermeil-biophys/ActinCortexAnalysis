
function CutCurves(CurveList)

for kc = 1:length(CurveList)
    ptrs = find(CurveList{kc} == '#');
    NameList{kc} = CurveList{kc}(1:ptrs(1)-1);
    Pos{kc} = CurveList{kc}(ptrs(1)+1:ptrs(2)-1);
    Time{kc} = CurveList{kc}(ptrs(2)+1:end);
    
end

listV2D = ls('D:\Data\MatFile\V2D');
listV2D = listV2D(3:end,:);

BOOL = 0;

for k = 1:size(listV2D,1)
    loadname = ['D:\Data\MatFile\V2D\' listV2D(k,:)];
    load(loadname);
    for kmr = 1:length(MR)
        if not(isempty(MR{kmr}))
            if ismember(MR{kmr}.name,NameList)
                ind = find(ismember(NameList,MR{kmr}.name));
                for kind = 1:length(ind)
                    
                    
                    if not(isempty(strfind(MR{kmr}.name,'R80')))
                        try
                            Nrmp = length(MR{kmr}.RampData{1});
                            
                            for krmp = 1:Nrmp
                                if not(ischar(MR{kmr}.RampData{1}{krmp}))
                                    RmpDatTmp = MR{kmr}.RampData{1}{krmp};
                                    
                                    if strcmp(Pos{ind(kind)},'ap')
                                        if RmpDatTmp(end,2)>str2double(Time{ind(kind)});
                                            MR{kmr}.RampData{1}{krmp} = 'Removed';
                                        end
                                    elseif strcmp(Pos{ind(kind)},'av')
                                        if RmpDatTmp(1,2)<str2double(Time{ind(kind)});
                                            MR{kmr}.RampData{1}{krmp} = 'Removed';
                                        end
                                    elseif strcmp(Pos{ind(kind)},'ex')
                                        ptrt = find(Time{ind(kind)} == '-');
                                        MinPos = str2double(Time{ind(kind)}(1:ptrt-1));
                                        MaxPos = str2double(Time{ind(kind)}(ptrt+1:end));
                                        if (RmpDatTmp(1,2)>MinPos)&&(RmpDatTmp(1,2)<MaxPos)
                                            MR{kmr}.RampData{1}{krmp} = 'Removed';
                                        end
                                    end
                                    
                                    
                                    clear MinPos MaxPos RmpDatTmp
                                    
                                    
                                end
                            end
                            
                            if strcmp(Pos{ind(kind)},'ap')
                                CutMask = MR{kmr}.time<str2double(Time{ind(kind)});
                            elseif strcmp(Pos{ind(kind)},'av')
                                CutMask = MR{kmr}.time>str2double(Time{ind(kind)});
                            elseif strcmp(Pos{ind(kind)},'ex')
                                ptrt = find(Time{ind(kind)} == '-');
                                MinPos = str2double(Time{ind(kind)}(1:ptrt-1));
                                MaxPos = str2double(Time{ind(kind)}(ptrt+1:end));
                                CutMask = (MR{kmr}.time<MinPos)|(MR{kmr}.time>MaxPos);
                                
                            end
                            
                            MR{kmr}.F      = MR{kmr}.F(CutMask);
                            MR{kmr}.D3     = MR{kmr}.D3(CutMask);
                            MR{kmr}.dz     = MR{kmr}.dz(CutMask);
                            MR{kmr}.time   = MR{kmr}.time(CutMask);
                            MR{kmr}.S      = MR{kmr}.S(CutMask);
                            
                            
                            clear CutMask*
                            
                            
                            if strcmp(Pos{ind(kind)},'ap')
                                CutMask = MR{kmr}.CstData(:,2)<str2double(Time{ind(kind)});
                                CutMask2 = MR{kmr}.RampData{2}(:,2)<str2double(Time{ind(kind)});
                            elseif strcmp(Pos{ind(kind)},'av')
                                CutMask = MR{kmr}.CstData(:,2)>str2double(Time{ind(kind)});
                                CutMask2 = MR{kmr}.RampData{2}(:,2)>str2double(Time{ind(kind)});
                            elseif strcmp(Pos{ind(kind)},'ex')
                                ptrt = find(Time{ind(kind)} == '-');
                                MinPos = str2double(Time{ind(kind)}(1:ptrt-1));
                                MaxPos = str2double(Time{ind(kind)}(ptrt+1:end));
                                CutMask = (MR{kmr}.CstData(:,2)<MinPos)|(MR{kmr}.CstData(:,2)>MaxPos);
                                CutMask2 = (MR{kmr}.RampData{2}(:,2)<MinPos)|(MR{kmr}.RampData{2}(:,2)>MaxPos);
                                
                            end
                            
                            MR{kmr}.CstData = MR{kmr}.CstData(CutMask,:);
                            MR{kmr}.RampData{2} = MR{kmr}.RampData{2}(CutMask2,:);
                            
                            clear CutMask*
                        catch
                            BOOL = 1;
                        end
                    end
                    
                    
                    if not(isempty(strfind(MR{kmr}.name,'5mT'))) || (BOOL==1)
                        
                        if strcmp(Pos{ind(kind)},'ap')
                            CutMask = MR{kmr}.time<str2double(Time{ind(kind)});
                        elseif strcmp(Pos{ind(kind)},'av')
                            CutMask = MR{kmr}.time>str2double(Time{ind(kind)});
                        elseif strcmp(Pos{ind(kind)},'ex')
                            ptrt = find(Time{ind(kind)} == '-');
                            MinPos = str2double(Time{ind(kind)}(1:ptrt-1));
                            MaxPos = str2double(Time{ind(kind)}(ptrt+1:end));
                            CutMask = (MR{kmr}.time<MinPos)|(MR{kmr}.time>MaxPos);
                            
                        end
                        
                        MR{kmr}.F      = MR{kmr}.F(CutMask);
                        MR{kmr}.D3     = MR{kmr}.D3(CutMask);
                        MR{kmr}.dz     = MR{kmr}.dz(CutMask);
                        MR{kmr}.time   = MR{kmr}.time(CutMask);
                        MR{kmr}.S      = MR{kmr}.S(CutMask);
                    end
                    
                    clear CutMask
                end
            end
        end
    end
    save(['D:\Data\MatFile\V2D\' listV2D(k,:)],'MR');
    
end











end