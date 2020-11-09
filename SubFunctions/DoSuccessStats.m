function DoSuccessStats(ExpTypes,exclusionvect)

load('D:\Data\MatFile\MecaDataTable.mat')

if ~exist('exclusionvect','var')
    
    exclusionvect = zeros(1,heigth(BigTable));
    
end

for ke = 1:length(ExpTypes)
    
    cprintf('_text',['\n\n\nFor experiments : ' ExpTypes{ke} '\n\n']) 
    
   
    ptrnotvalid = BigTable(:,'Validated').Variables == 0;  
    
    ptrexp = (BigTable(:,'ExpType').Variables == ExpTypes{ke}) & ~exclusionvect;
    ptrexpnotvalid = ptrexp & ptrnotvalid;    
     
    
    pctgtot = sum(ptrexpnotvalid)/sum(ptrexp)*100;
    
    cprintf('*red',['\n' num2str(round(pctgtot*10)/10) '%%'])
    cprintf('*blue',[' (' num2str(sum(ptrexpnotvalid)) '/' num2str(sum(ptrexp)) ')'])
    cprintf('text', ' of compressions not kept.\n\n') 
    
    reasons = BigTable(ptrexpnotvalid,'Comments').Variables;
    
    possiblereasons = unique(reasons);
    
    
    for kr = 1:length(possiblereasons)
        
        pctg = sum(reasons == possiblereasons{kr})/sum(ptrexp)*100;
        
        cprintf('*Red',['\n' num2str(round(pctg*10)/10) '%%'])
    cprintf('*blue',[' (' num2str(sum(reasons == possiblereasons{kr})) ')'])
        cprintf('text',[' ' possiblereasons{kr} '\n'])
        
    end
    
    
end


end
