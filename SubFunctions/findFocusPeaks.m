function [FPks,Sf] = findFocusPeaks(S,STD)

% STD = smooth(STD,4);

ptr_seuil = find(STD > 3000);

STD_seuil = STD(ptr_seuil);

S_seuil = S(ptr_seuil);

DS_seuil = diff(S_seuil);

i = 1;
k = 1;
S_tmp = [];
STD_tmp = [];
Sf = [];
STDf = [];

% plot(S,STD,'b-*')
% hold on

while k < length(DS_seuil)+1
    if ((DS_seuil(k) < 2) && (k ~= length(DS_seuil)))
        
       
            
     
        
        S_tmp   = [S_tmp S_seuil(k)];
        STD_tmp = [STD_tmp STD_seuil(k)];
        
        k       = k +1;
      
    else

        S_tmp = [S_tmp S_seuil(k)];
        STD_tmp = [STD_tmp STD_seuil(k)];
            
        
        Dstd = diff(STD_tmp)./abs(diff(STD_tmp));
        DDstd = abs(diff(Dstd));
        Mptr = find(DDstd~=0);
        
        SP = isequal(length(Mptr),1); % Booleen =1
        % pour un pic simple 
        DP = (isequal(length(Mptr),3)) && ...
            (all(diff(Mptr)>3)); % Booleen =1
        % pour un pic double
        
                
        if (length(STD_tmp) > 2) && SP
            [STDf_tmp,I] = max(STD_tmp);
            Sf_tmp = S_tmp(I);
            Sf = [Sf Sf_tmp];   
            
            l = length(S_tmp)-1;
            
            wi = 1;
            while (ptr_seuil(k-l)-wi > 0) && ...
                    (STD(ptr_seuil(k-l)-wi)<STD_tmp(1)) && ...
                    (S(ptr_seuil(k-l)-wi)==S_tmp(1)-1)
                S_tmp = [S(ptr_seuil(k-l)-wi) S_tmp];
                STD_tmp = [STD(ptr_seuil(k-l)-wi) STD_tmp];
                wi = wi+1;
            end  
                     
            wi = 1;
            while (ptr_seuil(k)+wi < length(STD)) && ...
                   (STD(ptr_seuil(k)+wi)<STD_tmp(end)) && ...
                   (S(ptr_seuil(k)+wi)==S_tmp(end)+1)
                S_tmp = [S_tmp S(ptr_seuil(k)+wi)];
                STD_tmp = [STD_tmp STD(ptr_seuil(k)+wi)];
                wi = wi+1;
            end    
            
            
            
            FPks{i}.STDc = STD_tmp; % data curve
            FPks{i}.Sc = S_tmp; % sampling curve
            FPks{i}.Sf = Sf_tmp; % approximate focus points
            i = i + 1;
            
% %             figure 
% %             hold on
%             plot(S_tmp,STD_tmp,'--o')
%             plot(FPks{i-1}.Sf,max(STD_tmp),'+')
% %             hold off
            
        elseif (length(STD_tmp) > 2) && DP
            
             l = length(S_tmp)-1;
            
            Smin = Mptr(2) +1;
            
            STD1_tmp = STD_tmp(1:Smin);
            STD2_tmp = STD_tmp(Smin:end);
            S1_tmp = S_tmp(1:Smin);
            S2_tmp = S_tmp(Smin:end);
            
            wi = 1;
            while (ptr_seuil(k-l)-wi > 0) && ...
                    (STD(ptr_seuil(k-l)-wi)<STD1_tmp(1)) && ...
                    (S(ptr_seuil(k-l)-wi)== S1_tmp(1)-1)
                S1_tmp = [S(ptr_seuil(k-l)-wi) S1_tmp];
                STD1_tmp = [STD(ptr_seuil(k-l)-wi) STD1_tmp];
                wi = wi+1;
            end  
                
            wi = 1;
            while (ptr_seuil(k)+wi < length(STD)) && ...
                   (STD(ptr_seuil(k)+wi)<STD2_tmp(end)) && ...
                   (S(ptr_seuil(k)+wi)==S2_tmp(end)+1)
                S2_tmp = [S2_tmp S(ptr_seuil(k)+wi)];
                STD2_tmp = [STD2_tmp STD(ptr_seuil(k)+wi)];
                wi = wi+1;
            end  
%             
            
            
            [STDf_tmp,I] = max(STD1_tmp);
            Sf = [Sf S1_tmp(I)];
            FPks{i}.STDc = STD1_tmp; % data curve
            FPks{i}.Sc = S1_tmp; % sampling curve
            FPks{i}.Sf = S1_tmp(I); % approximate focus points
            i = i + 1;
            
            [STDf_tmp,I] = max(STD2_tmp);
            Sf = [Sf S2_tmp(I)];
            FPks{i}.STDc = STD2_tmp; % data curve
            FPks{i}.Sc = S2_tmp; % sampling curve
            FPks{i}.Sf = S2_tmp(I); % approximate focus points
            i = i + 1;
            
% %             figure 
% %             hold on
%             plot(S1_tmp,STD1_tmp,'--o')
%             plot(S2_tmp,STD2_tmp,'--o')
% %             plot(FPks{i-1}.Sf,max(STD_tmp),'+')
% %             plot(FPks{i-2}.Sf,max(STD_tmp),'+')
% %             hold off
         
        end
        S_tmp = [];
        STD_tmp = [];
        k = k +1;
    end
end
Sf = Sf';


end
