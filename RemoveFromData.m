function RemoveFromData(Folder,NameList,specif)

RemovedList = {};

listV2D = ls([Folder filesep 'V2D']);
listV2D = listV2D(3:end,:);

for k = 1:size(listV2D,1)
    if sum(strfind(listV2D(k,:),specif))>0
    load(['D:\Data\MatFile\V2D\' listV2D(k,:)])
    for kmr = 1:length(MR)
        if not(isempty(MR{kmr}))
            if ismember(MR{kmr}.name,NameList)
                RemovedList = {RemovedList{:},MR{kmr}.name};
                MR{kmr} = [];
            end
        end
    end
    save(['D:\Data\MatFile\V2D\' listV2D(k,:)],'MR');
    end
end

fprintf('Following data could not be removed :\n\n')

NotRemovedList = char(setdiff(NameList,RemovedList))











end