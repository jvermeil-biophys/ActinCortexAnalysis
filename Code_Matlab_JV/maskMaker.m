%% create a mask in CleWin5

%% Useful commands
% rectangle(x1, y1, x2, y2)
% box(x, y, width, height, angle)
% polygon(nodes)
% wire(style, width, nodes)
% circle(x, y, radius)
% ring(x, y, radius, width)
% text(text, matrix)
% symbol(symbolname, matrix)

%% Zone Rectangle Drawer

S1 = 101600-2000;
O1 = -S1/2;
s1_3 = 33860;

S2 = 24000;
margin_2 = (s1_3 - S2)/2;
s2_2 = S2/2;
s2_3 = S2/3;

margin_3 = 200;

for x1=0:2
    for y1=0:2
        O2x = O1 + x1*s1_3;
        O2y = O1 + y1*s1_3;
        
        
        O3x_5 = O2x + margin_2;
        O3y_5 = O2y + margin_2 + s2_3*2 + margin_3;
        rectangle(O3x_5, O3y_5, O3x_5 + 7800, O3y_5 + 7800);
        
        O3x_10 = O2x + margin_2;
        O3y_10 = O2y + margin_2 + s2_3*1 + margin_3/2;
        rectangle(O3x_10, O3y_10, O3x_10 + 7800, O3y_10 + 7800);
        
        O3x_15 = O2x + margin_2;
        O3y_15 = O2y + margin_2;
        rectangle(O3x_15, O3y_15, O3x_15 + 7800, O3y_15 + 7800);
        
        O3x_20 = O2x + margin_2 + s2_3*2 + margin_3;
        O3y_20 = O2y + margin_2 + s2_3*1 + margin_3/2;
        rectangle(O3x_20, O3y_20, O3x_20 + 7800, O3y_20 + 7800);
        
        O3x_25 = O2x + margin_2 + s2_3*2 + margin_3;
        O3y_25 = O2y + margin_2 + s2_3*2 + margin_3;
        rectangle(O3x_25, O3y_25, O3x_25 + 7800, O3y_25 + 7800);
        
        O3x_30 = O2x + margin_2 + s2_3*2 + margin_3;
        O3y_30 = O2y + margin_2;
        rectangle(O3x_30, O3y_30, O3x_30 + 7800, O3y_30 + 7800);
        
        O3x_mixed1 = O2x + margin_2 + s2_3*1 + margin_3/2;
        O3y_mixed1 = O2y + margin_2 + s2_2*1 + margin_3*1.5;
        rectangle(O3x_mixed1, O3y_mixed1, O3x_mixed1 + 7800, O3y_mixed1 + 11850);
        
        O3x_mixed2 = O2x + margin_2 + s2_3*1 + margin_3/2;
        O3y_mixed2 = O2y + margin_2;
        rectangle(O3x_mixed2, O3y_mixed2, O3x_mixed2 + 7800, O3y_mixed2 + 11850);
        
    end
end

%% limits (1/2)

S1 = 101600-2000;
O1 = -S1/2;
s1_3 = 33860;

S2 = 24000;
margin_2 = (s1_3 - S2)/2;
s2_2 = S2/2;
s2_3 = S2/3;

margin_3 = 200;

for x1=0:2
    for y1=0:2
        O2x = O1 + x1*s1_3;
        O2y = O1 + y1*s1_3;
        
        O3x_5 = O2x + margin_2;
        O3y_5 = O2y + margin_2 + s2_3*2 + margin_3;
        
        O3x_10 = O2x + margin_2;
        O3y_10 = O2y + margin_2 + s2_3*1 + margin_3/2;
        
        O3x_15 = O2x + margin_2;
        O3y_15 = O2y + margin_2;
        
        O3x_20 = O2x + margin_2 + s2_3*2 + margin_3;
        O3y_20 = O2y + margin_2 + s2_3*1 + margin_3/2;
        
        O3x_25 = O2x + margin_2 + s2_3*2 + margin_3;
        O3y_25 = O2y + margin_2 + s2_3*2 + margin_3;
        
        O3x_30 = O2x + margin_2 + s2_3*2 + margin_3;
        O3y_30 = O2y + margin_2;
        
        O3x_mixed1 = O2x + margin_2 + s2_3*1 + margin_3/2;
        O3y_mixed1 = O2y + margin_2 + s2_2*1 + margin_3*1.5;
        
        O3x_mixed2 = O2x + margin_2 + s2_3*1 + margin_3/2;
        O3y_mixed2 = O2y + margin_2;
        
        O3x_list = [O3x_5 O3x_10 O3x_15 O3x_20 O3x_25 O3x_30];
        O3y_list = [O3y_5 O3y_10 O3y_15 O3y_20 O3y_25 O3y_30];
%         O3x_list = [O3x_5 O3x_10 O3x_15];
%         O3y_list = [O3y_5 O3y_10 O3y_15];
        
        nL = 7800/120;
        wid = 10;
        len = 60;
        
        for kO = 1:length(O3x_list)
            O3x = O3x_list(kO);
            O3y = O3y_list(kO);
            for iLV = 0:5
                for jLV = 0:nL
                    xLV = O3x + iLV * 7800/5;
                    yLV = O3y + jLV * 7800/nL;
                    rectangle(xLV - wid/2, yLV - len/2, xLV + wid/2, yLV + len/2);
                end
            end
            
            for iLH = 0:nL
                for jLH = 0:10
                    xLH = O3x + iLH * 7800/nL;
                    yLH = O3y + jLH * 7800/10;
                    rectangle(xLH - len/2, yLH - wid/2, xLH + len/2, yLH + wid/2);
                end
            end
            
        end
        
    end
end
