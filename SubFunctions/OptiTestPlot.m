pixelconv=1/15.8;

PC = [0 3 5 6 7 8 10]; % pourcentage optiprep

Dist = cell(length(PC),1);
Group = [];

Nb = [1 2 3 4 5]; % n° d'exp




for i = 1:length(PC)
    for ii = 1:length(Nb)
        LOOP = 1;
        stackname = ['D:\Data\Raw\18.05.24_OptiTest-DC\DC-OptiTest_' num2str(PC(i)) 'pc_' num2str(Nb(ii))];
        kimoname = ['Kimo_DCmed-' num2str(PC(i)) 'pc_170518'];
        
        filename = [stackname '_Results.txt'];
        if not(exist(filename))
            filename = [stackname '-1' '_Results.txt'];
            LOOP = LOOP + 1;
            if not(exist(filename))
                error('ERROR : No file?????')
            end
        end
        
       
        while LOOP > 0
            
         filename
        [VarName1,Area,Mean,StdDev,Min,Max,X,Y,XM,YM,Major,Minor,Angle,Slice] = importfileDC_OptiTest(filename);
        
            
            [X1,Y1,S1,ptr1,X2,Y2,S2,ptr2] = splitBeads_Fluo([XM Slice],stackname,YM);
            
            STD1 = StdDev(ptr1);
            STD2 = StdDev(ptr2);
            
            Z1 = Zcorr(STD1,X1,Y1,S1,stackname,kimoname);
            Z2 = Zcorr(STD2,X2,Y2,S2,stackname,kimoname);
            
            X1 = X1*pixelconv;
            X2 = X2*pixelconv;
            Y1 = Y1*pixelconv;
            Y2 = Y2*pixelconv;
            
            D = sqrt((X1-X2).^2+(Y1-Y2).^2+(Z1-Z2).^2);
            D = D*1000 - 4432;
            
            
            figure
            plot(S1,D)
            
            Dist{i} = [Dist{i} mean(D)];
            Group = [Group PC(i)];
            
            clear D X1 X2 Y1 Y2 Z1 Z2
            
            if LOOP > 1
                filename = [stackname '-' num2str(LOOP) '_Results.txt'];
                if not(exist(filename))
                    LOOP = 0;
                else
                    LOOP = LOOP + 1;
                end
                
            else
                LOOP = 0;
            end
            
            
        end
        
    end
end


boxplot([Dist{:}],Group)

