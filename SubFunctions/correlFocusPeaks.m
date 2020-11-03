function [dt,Sxy] = correlFocusPeaks(FPks1,FPks2,ptr);

l = min(length(FPks1),length(FPks2));
Sxy = [];
dt = [];

for i=1:l
    
    STD1 = FPks1{ptr(1,i)}.STDc;
    STD2 = FPks2{ptr(2,i)}.STDc;
    S1 = FPks1{ptr(1,i)}.Sc;
    S2 = FPks2{ptr(2,i)}.Sc;
    
    
    [S,p1,p2] = intersect(S1,S2);
    
    if length(S) > 5
        
        STD1 = STD1(p1);
        STD2 = STD2(p2);
        
        STD1 = (STD1-min(STD1))/max(STD1-min(STD1));
        STD2 = (STD2-min(STD2))/max(STD2-min(STD2));
        
%                 figure
%                 hold on
%                 plot(S,STD1,'-*')
%                 plot(S,STD2,'-*')
%                 xlim([min(S)-1 max(S)+1])
%                 hold off
%         
        
    
       
    [val,lag] = xcorr(STD1,STD2);
    lagint = lag(1):0.01:lag(end); % x for interp
    valint = interp1(lag,val,lagint,'spline');
    [valm,I] = max(valint);
    dt_tmp = lagint(I);
    
    
    
    
%         figure
%         hold on
%         plot(lagint,valint,'-o')
%         plot(lag,val,'-*')
%         hold off
%     
%         figure
%         hold on
%         plot(S-dt_tmp,STD1,'-*')
%         plot(S,STD2,'-*')
%         xlim([min(S)-1 max(S)+1])
%         hold off
    
    dt = [dt dt_tmp];
    
    Sxy_tmp = mean(FPks2{ptr(2,i)}.Sf,FPks1{ptr(1,i)}.Sf);
    Sxy = [Sxy Sxy_tmp];
    
    
    end
    
    [Sxy, ptrS] = sort(Sxy);
    dt = dt(ptrS);
    
end

end