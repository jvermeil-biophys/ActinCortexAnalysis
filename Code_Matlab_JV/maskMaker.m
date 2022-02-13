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

S2 = 24600;
margin_2 = (s1_3 - S2)/2;
s2_2 = S2/2;
s2_3 = S2/3;

margin_3 = 200;
S3 = 8000;

for x1=0:2
    for y1=0:2
        O2x = O1 + x1*s1_3;
        O2y = O1 + y1*s1_3;
        
        
        O3x_5 = O2x + margin_2;
        O3y_5 = O2y + margin_2 + s2_3*2 + margin_3;
        rectangle(O3x_5, O3y_5, O3x_5 + S3, O3y_5 + S3);
        
        O3x_10 = O2x + margin_2;
        O3y_10 = O2y + margin_2 + s2_3*1 + margin_3/2;
        rectangle(O3x_10, O3y_10, O3x_10 + S3, O3y_10 + S3);
        
        O3x_15 = O2x + margin_2;
        O3y_15 = O2y + margin_2;
        rectangle(O3x_15, O3y_15, O3x_15 + S3, O3y_15 + S3);
        
        O3x_20 = O2x + margin_2 + s2_3*2 + margin_3;
        O3y_20 = O2y + margin_2 + s2_3*1 + margin_3/2;
        rectangle(O3x_20, O3y_20, O3x_20 + S3, O3y_20 + S3);
        
        O3x_25 = O2x + margin_2 + s2_3*2 + margin_3;
        O3y_25 = O2y + margin_2 + s2_3*2 + margin_3;
        rectangle(O3x_25, O3y_25, O3x_25 + S3, O3y_25 + S3);
        
        O3x_30 = O2x + margin_2 + s2_3*2 + margin_3;
        O3y_30 = O2y + margin_2;
        rectangle(O3x_30, O3y_30, O3x_30 + S3, O3y_30 + S3);
        
        O3x_mixed1 = O2x + margin_2 + s2_3*1 + margin_3/2;
        O3y_mixed1 = O2y + margin_2 + s2_2*1 + 800 + margin_3*0.75;
        rectangle(O3x_mixed1, O3y_mixed1, O3x_mixed1 + S3, O3y_mixed1 + ((S2/2)-margin_3*1.5));
        
        O3x_mixed2 = O2x + margin_2 + s2_3*1 + margin_3/2;
        O3y_mixed2 = O2y + margin_2 + margin_3*0.75;
        rectangle(O3x_mixed2, O3y_mixed2, O3x_mixed2 + S3, O3y_mixed2 + ((S2/2)-margin_3*1.5) + 800);
        
    end
end

%% Text
t = '20 µm';
xt = 1000;
yt = 1000;
matrix = [0.1 0 xt; 0 0.1 yt; 0 0 1];
text(t, matrix)


%% Limits

S1 = 101600-2000;
O1 = -S1/2;
s1_3 = 33860;

S2 = 24600;
margin_2 = (s1_3 - S2)/2;
s2_2 = S2/2;
s2_3 = S2/3;

margin_3 = 200;
S3 = 8000;

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

        O3x_list = [O3x_5 O3x_10 O3x_15 O3x_20 O3x_25 O3x_30];
        O3y_list = [O3y_5 O3y_10 O3y_15 O3y_20 O3y_25 O3y_30];
        
        space = 100;
        nL = S3/space;
        wid = 10;
        len = 50;
        nCol3 = 5;
        nRow3 = 10;
        
        for kO = 1:length(O3x_list)
            O3x = O3x_list(kO);
            O3y = O3y_list(kO);
            for iLV = 0:nCol3
                for jLV = 0:nL
                    xLV = O3x + iLV * S3/nCol3;
                    yLV = O3y + jLV * S3/nL;
                    rectangle(xLV - wid/2, yLV - len/2, xLV + wid/2, yLV + len/2);
                end
            end
            
            for iLH = 0:nL
                for jLH = 0:nRow3
                    xLH = O3x + iLH * S3/nL;
                    yLH = O3y + jLH * S3/nRow3;
                    rectangle(xLH - len/2, yLH - wid/2, xLH + len/2, yLH + wid/2);
                end
            end
        end
    end
end
  


S1 = 101600-2000;
O1 = -S1/2;
s1_3 = 33860;

S2 = 24600;
margin_2 = (s1_3 - S2)/2;
s2_2 = S2/2;
s2_3 = S2/3;

margin_3 = 200;
S3 = 8000;

for x1=0:2
    for y1=0:2
        O2x = O1 + x1*s1_3;
        O2y = O1 + y1*s1_3;
        
        O3x_mixed1 = O2x + margin_2 + s2_3*1 + margin_3/2;
        O3y_mixed1 = O2y + margin_2 + s2_2*1 + 800 + margin_3*0.75;

        O3x_mixed2 = O2x + margin_2 + s2_3*1 + margin_3/2;
        O3y_mixed2 = O2y + margin_2 + margin_3*0.75;

        O3x_list = [O3x_mixed1 O3x_mixed2];
        O3y_list = [O3y_mixed1 O3y_mixed2];
        
        space = 100;
        W = S3;
        H = [((S2/2)-margin_3*1.5)-800  ((S2/2)-margin_3*1.5)+800];
        nLx = W/space;
        nLy = H./space;
        wid = 10;
        len = 50;
        nRow3 = [14 16];
        nCol3 = 5;
        
        for kO = 1:length(O3x_list)
            
            O3x = O3x_list(kO);
            O3y = O3y_list(kO);
            for iLV = 0:nCol3
                for jLV = 0:nLy(kO)
                    xLV = O3x + iLV * W/nCol3;
                    yLV = O3y + jLV * H(kO)/nLy(kO);
                    rectangle(xLV - wid/2, yLV - len/2, xLV + wid/2, yLV + len/2);
                end
            end
            
            for iLH = 0:nLx
                for jLH = 0:nRow3(kO)
                    xLH = O3x + iLH * W/nLx;
                    yLH = O3y + jLH * H(kO)/nRow3(kO);
                    rectangle(xLH - len/2, yLH - wid/2, xLH + len/2, yLH + wid/2);
                end
            end
        end
    end
end

%% Test Limits

S1 = 101600-2000;
O1 = -S1/2;
s1_3 = 33860;

S2 = 24600;
margin_2 = (s1_3 - S2)/2;
s2_2 = S2/2;
s2_3 = S2/3;

margin_3 = 200;
S3 = 8000;

for x1=0:0
    for y1=0:0
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
        O3y_mixed1 = O2y + margin_2 + s2_2*1 + 800 + margin_3*0.75;
        
        O3x_mixed2 = O2x + margin_2 + s2_3*1 + margin_3/2;
        O3y_mixed2 = O2y + margin_2 + margin_3*0.75;
        
        
        
        O3x_list = [O3x_5 O3x_10 O3x_15 O3x_20 O3x_25 O3x_30];
        O3y_list = [O3y_5 O3y_10 O3y_15 O3y_20 O3y_25 O3y_30];
        
        space = 100;
        nL = S3/space;
        wid = 10;
        len = 50;
        nCol3 = 5;
        nRow3 = 10;
        
        for kO = 1:length(O3x_list)
            O3x = O3x_list(kO);
            O3y = O3y_list(kO);
            for iLV = 0:nCol3
                for jLV = 0:nL
                    xLV = O3x + iLV * S3/nCol3;
                    yLV = O3y + jLV * S3/nL;
%                     cprintf(xLV - wid/2, yLV - len/2, xLV + wid/2, yLV + len/2);
                end
            end
            
            for iLH = 0:nL
                for jLH = 0:nRow3
                    xLH = O3x + iLH * S3/nL;
                    yLH = O3y + jLH * S3/nRow3;
%                     cprintf(xLH - len/2, yLH - wid/2, xLH + len/2, yLH + wid/2);
                end
            end
        end
        
        
        
        O3x_list = [O3x_mixed1 O3x_mixed2];
        O3y_list = [O3y_mixed1 O3y_mixed2];
        
        space = 100;
        W = S3;
        H = [((S2/2)-margin_3*1.5)-800  ((S2/2)-margin_3*1.5)+800];
        nLx = W/space;
        nLy = H./space;
        wid = 10;
        len = 50;
        nRow3 = [14 16];
        nCol3 = 5;
        
        for kO = 1:length(O3x_list)
            
            O3x = O3x_list(kO);
            O3y = O3y_list(kO);
            for iLV = 0:nCol3
                for jLV = 0:nLy(kO)
                    xLV = O3x + iLV * W/nCol3;
                    yLV = O3y + jLV * H(kO)/nLy(kO);
%                     rectangle(xLV - wid/2, yLV - len/2, xLV + wid/2, yLV + len/2);
                end
            end
            
            for iLH = 0:nLx
                for jLH = 0:nRow3(kO)
                    xLH = O3x + iLH * W/nLx;
                    yLH = O3y + jLH * H(kO)/nRow3(kO);
%                     rectangle(xLH - len/2, yLH - wid/2, xLH + len/2, yLH + wid/2);
                end
            end
        end
        
        
    end
end

%% 5um

S1 = 101600-2000;
O1 = -S1/2;
s1_3 = 33860;

S2 = 24600;
margin_2 = (s1_3 - S2)/2;
s2_2 = S2/2;
s2_3 = S2/3;

margin_3 = 200;
S3 = 8000;

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
        O3y_mixed1 = O2y + margin_2 + s2_2*1 + margin_3*0.75;
        
        O3x_mixed2 = O2x + margin_2 + s2_3*1 + margin_3/2;
        O3y_mixed2 = O2y + margin_2 + margin_3*0.75;
        
        nCol3 = 5;
        nRow3 = 10;
        
        W3 = S3 / nCol3;
        H3 = S3 / nRow3;
        diam = 5;
        space = 60;
        
        nCol4 = floor(W3/space) - 1;
        margin_x4 = (W3 - ((nCol4-1)*space)) / 2;
        nRow4 = floor(H3/space) - 1;
        margin_y4 = (H3 - ((nRow4-1)*space)) / 2;
        
        for x3=1:nCol3
            for y3 = 1:nRow3
                Ox4_5 = O3x_5 + (x3-1) * W3;
                Oy4_5 = O3y_5 + (y3-1) * H3;
                for x4 = 1:nCol4
                    for y4 = 1:nRow4
                        xC = Ox4_5 + margin_x4 + (x4-1)*space;
                        yC = Oy4_5 + margin_y4 + (y4-1)*space;
                        circle(xC, yC, diam/2);
                    end
                end
            end
        end
        
        
    end
end

%% 10um

S1 = 101600-2000;
O1 = -S1/2;
s1_3 = 33860;

S2 = 24600;
margin_2 = (s1_3 - S2)/2;
s2_2 = S2/2;
s2_3 = S2/3;

margin_3 = 200;
S3 = 8000;

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
        O3y_mixed1 = O2y + margin_2 + s2_2*1 + margin_3*0.75;
        
        O3x_mixed2 = O2x + margin_2 + s2_3*1 + margin_3/2;
        O3y_mixed2 = O2y + margin_2 + margin_3*0.75;
        
        nCol3 = 5;
        nRow3 = 10;
        
        W3 = S3 / nCol3;
        H3 = S3 / nRow3;
        diam = 10;
        space = 60;
        
        nCol4 = floor(W3/space) - 1;
        margin_x4 = (W3 - ((nCol4-1)*space)) / 2;
        nRow4 = floor(H3/space) - 1;
        margin_y4 = (H3 - ((nRow4-1)*space)) / 2;
        
        for x3=1:nCol3
            for y3 = 1:nRow3
                Ox4_10 = O3x_10 + (x3-1) * W3;
                Oy4_10 = O3y_10 + (y3-1) * H3;
                for x4 = 1:nCol4
                    for y4 = 1:nRow4
                        xC = Ox4_10 + margin_x4 + (x4-1)*space;
                        yC = Oy4_10 + margin_y4 + (y4-1)*space;
                        circle(xC, yC, diam/2);
                    end
                end
            end
        end
        
        
    end
end


%% 15um

S1 = 101600-2000;
O1 = -S1/2;
s1_3 = 33860;

S2 = 24600;
margin_2 = (s1_3 - S2)/2;
s2_2 = S2/2;
s2_3 = S2/3;

margin_3 = 200;
S3 = 8000;

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
        O3y_mixed1 = O2y + margin_2 + s2_2*1 + margin_3*0.75;
        
        O3x_mixed2 = O2x + margin_2 + s2_3*1 + margin_3/2;
        O3y_mixed2 = O2y + margin_2 + margin_3*0.75;
        
        nCol3 = 5;
        nRow3 = 10;
        
        W3 = S3 / nCol3;
        H3 = S3 / nRow3;
        diam = 15;
        space = 70;
        
        nCol4 = floor(W3/space) - 1;
        margin_x4 = (W3 - ((nCol4-1)*space)) / 2;
        nRow4 = floor(H3/space) - 1;
        margin_y4 = (H3 - ((nRow4-1)*space)) / 2;
        
        for x3=1:nCol3
            for y3 = 1:nRow3
                Ox4_15 = O3x_15 + (x3-1) * W3;
                Oy4_15 = O3y_15 + (y3-1) * H3;
                for x4 = 1:nCol4
                    for y4 = 1:nRow4
                        xC = Ox4_15 + margin_x4 + (x4-1)*space;
                        yC = Oy4_15 + margin_y4 + (y4-1)*space;
                        circle(xC, yC, diam/2);
                    end
                end
            end
        end
        
        
    end
end

%% 20um

S1 = 101600-2000;
O1 = -S1/2;
s1_3 = 33860;

S2 = 24600;
margin_2 = (s1_3 - S2)/2;
s2_2 = S2/2;
s2_3 = S2/3;

margin_3 = 200;
S3 = 8000;

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
        O3y_mixed1 = O2y + margin_2 + s2_2*1 + margin_3*0.75;
        
        O3x_mixed2 = O2x + margin_2 + s2_3*1 + margin_3/2;
        O3y_mixed2 = O2y + margin_2 + margin_3*0.75;
        
        nCol3 = 5;
        nRow3 = 10;
        
        W3 = S3 / nCol3;
        H3 = S3 / nRow3;
        diam = 20;
        space = 70;
        
        nCol4 = floor(W3/space) - 1;
        margin_x4 = (W3 - ((nCol4-1)*space)) / 2;
        nRow4 = floor(H3/space) - 1;
        margin_y4 = (H3 - ((nRow4-1)*space)) / 2;
        
        for x3=1:nCol3
            for y3 = 1:nRow3
                Ox4_20 = O3x_20 + (x3-1) * W3;
                Oy4_20 = O3y_20 + (y3-1) * H3;
                for x4 = 1:nCol4
                    for y4 = 1:nRow4
                        xC = Ox4_20 + margin_x4 + (x4-1)*space;
                        yC = Oy4_20 + margin_y4 + (y4-1)*space;
                        circle(xC, yC, diam/2);
                    end
                end
            end
        end
        
        
    end
end

%% 25um

S1 = 101600-2000;
O1 = -S1/2;
s1_3 = 33860;

S2 = 24600;
margin_2 = (s1_3 - S2)/2;
s2_2 = S2/2;
s2_3 = S2/3;

margin_3 = 200;
S3 = 8000;

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
        O3y_mixed1 = O2y + margin_2 + s2_2*1 + margin_3*0.75;
        
        O3x_mixed2 = O2x + margin_2 + s2_3*1 + margin_3/2;
        O3y_mixed2 = O2y + margin_2 + margin_3*0.75;
        
        nCol3 = 5;
        nRow3 = 10;
        
        W3 = S3 / nCol3;
        H3 = S3 / nRow3;
        diam = 25;
        space = 80;
        
        nCol4 = floor(W3/space) - 1;
        margin_x4 = (W3 - ((nCol4-1)*space)) / 2;
        nRow4 = floor(H3/space) - 1;
        margin_y4 = (H3 - ((nRow4-1)*space)) / 2;
        
        for x3=1:nCol3
            for y3 = 1:nRow3
                Ox4_25 = O3x_25 + (x3-1) * W3;
                Oy4_25 = O3y_25 + (y3-1) * H3;
                for x4 = 1:nCol4
                    for y4 = 1:nRow4
                        xC = Ox4_25 + margin_x4 + (x4-1)*space;
                        yC = Oy4_25 + margin_y4 + (y4-1)*space;
                        circle(xC, yC, diam/2);
                    end
                end
            end
        end
        
        
    end
end


%% 30um

S1 = 101600-2000;
O1 = -S1/2;
s1_3 = 33860;

S2 = 24600;
margin_2 = (s1_3 - S2)/2;
s2_2 = S2/2;
s2_3 = S2/3;

margin_3 = 200;
S3 = 8000;

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
        O3y_mixed1 = O2y + margin_2 + s2_2*1 + margin_3*0.75;
        
        O3x_mixed2 = O2x + margin_2 + s2_3*1 + margin_3/2;
        O3y_mixed2 = O2y + margin_2 + margin_3*0.75;
        
        nCol3 = 5;
        nRow3 = 10;
        
        W3 = S3 / nCol3;
        H3 = S3 / nRow3;
        diam = 30;
        space = 80;
        
        nCol4 = floor(W3/space) - 1;
        margin_x4 = (W3 - ((nCol4-1)*space)) / 2;
        nRow4 = floor(H3/space) - 1;
        margin_y4 = (H3 - ((nRow4-1)*space)) / 2;
        
        for x3=1:nCol3
            for y3 = 1:nRow3
                Ox4_30 = O3x_30 + (x3-1) * W3;
                Oy4_30 = O3y_30 + (y3-1) * H3;
                for x4 = 1:nCol4
                    for y4 = 1:nRow4
                        xC = Ox4_30 + margin_x4 + (x4-1)*space;
                        yC = Oy4_30 + margin_y4 + (y4-1)*space;
                        circle(xC, yC, diam/2);
                    end
                end
            end
        end
        
        
    end
end


%% mixed1

S1 = 101600-2000;
O1 = -S1/2;
s1_3 = 33860;

S2 = 24600;
margin_2 = (s1_3 - S2)/2;
s2_2 = S2/2;
s2_3 = S2/3;

margin_3 = 200;
S3 = 8000;

for x1=0:2
    for y1=0:2
        O2x = O1 + x1*s1_3;
        O2y = O1 + y1*s1_3;
        
        O3x_mixed1 = O2x + margin_2 + s2_3*1 + margin_3/2;
        O3y_mixed1 = O2y + margin_2 + s2_2*1 + 800 + margin_3*0.75;

        O3x_mixed2 = O2x + margin_2 + s2_3*1 + margin_3/2;
        O3y_mixed2 = O2y + margin_2 + margin_3*0.75;
        
        W = S3;
        H = ((S2/2)-margin_3*1.5)-800;
        nRow3 = 14;
        nCol3 = 5;
        
        W3 = W / nCol3;
        H3 = H / nRow3;
        
        for x3=1:nCol3
            for y3 = 1:nRow3
                O4x_mixed1 = O3x_mixed1 + (x3-1) * W3;
                O4y_mixed1 = O3y_mixed1 + (y3-1) * H3;
                dX4 = [90 80 80 80 70 70 70 70 60 60];
                dY4 = [90 90 80 80 70 70 70 70 60 60];
                current_y4 = O4y_mixed1;
                D4 =  [40 35 35 30 25 25 20 15 10 5];
                for i_Size = 1:length(dY4)
                    dX = dX4(i_Size);
                    nX = floor(W3/dX) - 1;
                    margin_x = (W3 - ((nX-1)*dX)) / 2;
                    current_y4 = current_y4 + dY4(i_Size);
                    for j = 1:nX
                        xC = O4x_mixed1 + margin_x + (j-1)*dX;
                        yC = current_y4;
                        circle(xC, yC, D4(i_Size)/2);
                    end
                end
            end
        end
    end
end

%% mixed2

S1 = 101600-2000;
O1 = -S1/2;
s1_3 = 33860;

S2 = 24600;
margin_2 = (s1_3 - S2)/2;
s2_2 = S2/2;
s2_3 = S2/3;

margin_3 = 200;
S3 = 8000;

for x1=0:2
    for y1=0:2
        O2x = O1 + x1*s1_3;
        O2y = O1 + y1*s1_3;
        
        O3x_mixed1 = O2x + margin_2 + s2_3*1 + margin_3/2;
        O3y_mixed1 = O2y + margin_2 + s2_2*1 + 800 + margin_3*0.75;

        O3x_mixed2 = O2x + margin_2 + s2_3*1 + margin_3/2;
        O3y_mixed2 = O2y + margin_2 + margin_3*0.75;
        
        W = S3;
        H = ((S2/2)-margin_3*1.5)+800;
        nRow3 = 16;
        nCol3 = 5;
        
        W3 = W / nCol3;
        H3 = H / nRow3;
        
        D3 =     [5  10 15 20 25 30 35 40 40 35 30 25 20 15 10 5];
        space3 = [60 60 70 70 80 80 90 90 90 90 80 80 70 70 60 60];
        for x3=1:nCol3
            for y3 = 1:nRow3
                O4x_mixed2 = O3x_mixed2 + (x3-1) * W3;
                O4y_mixed2 = O3y_mixed2 + (y3-1) * H3;
                
                diam = D3(y3);
                space = space3(y3);
                
                nCol4 = floor(W3/space) - 1;
                margin_x4 = (W3 - ((nCol4-1)*space)) / 2;
                nRow4 = floor(H3/space) - 1;
                margin_y4 = (H3 - ((nRow4-1)*space)) / 2;
                
                for x4 = 1:nCol4
                    for y4 = 1:nRow4
                        xC = O4x_mixed2 + margin_x4 + (x4-1)*space;
                        yC = O4y_mixed2 + margin_y4 + (y4-1)*space;
                        circle(xC, yC, diam/2);
                    end
                end
            end
        end
    end
end