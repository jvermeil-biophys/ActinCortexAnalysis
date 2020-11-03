clear all
close all

stepnm = 20*1.33/1.52;

Bigradius = 4503/2;

Smallradius = 1093;

IMSHOW = 0;

%% grosse bille collée
% 
% names = {'1' '2-1' '2-2'}
% imgnames = {'1' '2' '2'}
% 
% 
% for k = 1:length(names)
% [VarName, Area, StdDev, XM, YM, Slice] = importfile100620(['D:\Data\Raw\20.06.09\DepthoCouplé_GrosseCollée_' names{k} '_Results.txt']);
% 
% 
% 
% ptrsmall = find(Area < 50);
% 
% ptrbig = find(Area > 200);
% 
% 
% [~, indfocussmall] = max(StdDev(ptrsmall));
% 
% focussmall = Slice(ptrsmall(indfocussmall));
% 
% 
% [~, indfocusbig] = max(StdDev(ptrbig));
% 
% focusbig = Slice(ptrbig(indfocusbig));
% 
% 
% Xb = XM(ptrbig(indfocusbig));
% Yb = YM(ptrbig(indfocusbig));
% 
% Xs = XM(ptrsmall(indfocussmall));
% Ys = YM(ptrsmall(indfocussmall));
% 
% 
% Ib = imadjust(imread(['D:\Data\Raw\20.06.09\DepthoCouplé_GrosseCollée_' imgnames{k} '.tif'],focusbig));
% Is = imadjust(imread(['D:\Data\Raw\20.06.09\DepthoCouplé_GrosseCollée_' imgnames{k} '.tif'],focussmall));
% 
% figure
% subplot(121)
% hold on
% imshow(Ib)
% plot(Xb,Yb,'r*')
% subplot(122)
% hold on
% imshow(Is)
% plot(Xs,Ys,'r*')
% 
% focusdist(k) = abs(focusbig-focussmall)*stepnm/1000;
% 
% end
% 
% clear names ptrsmall ptrbig focussmall focusbig indfocussmall indfocusbig

%% tout collé

% 09-06
names = {'1-1' '1-2' '1-3' '1-4' '2-1' '2-2' '2-3' '2-4'};
imgnames = {'1' '1' '1bis' '1bis' '2-1' '2-2' '2-3' '2-4'};


for k = 1:length(names)
[VarName, Area, StdDev, XM, YM, Slice] = importfile100620(['D:\Data\Raw\20.06.09\DepthoCouplé_ToutCollé_' names{k} '_Results.txt']);



ptrsmall = find(Area < 50);

ptrbig = find(Area > 200);


[~, indfocussmall] = max(StdDev(ptrsmall));

focussmall = Slice(ptrsmall(indfocussmall));


[~, indfocusbig] = max(StdDev(ptrbig));

focusbig = Slice(ptrbig(indfocusbig));


if IMSHOW
Xb = XM(ptrbig(indfocusbig));
Yb = YM(ptrbig(indfocusbig));

Xs = XM(ptrsmall(indfocussmall));
Ys = YM(ptrsmall(indfocussmall));



Ib = imadjust(imread(['D:\Data\Raw\20.06.09\DepthoCouplé_ToutCollé_' imgnames{k} '.tif'],focusbig));
Is = imadjust(imread(['D:\Data\Raw\20.06.09\DepthoCouplé_ToutCollé_' imgnames{k} '.tif'],focussmall));

figure
subplot(121)
hold on
imshow(Ib)
plot(Xb,Yb,'r*')
subplot(122)
hold on
imshow(Is)
plot(Xs,Ys,'r*')

end


focusdist2(k) = (abs(focusbig-focussmall)*stepnm + (Bigradius-Smallradius))/1000;

end



histogram(focusdist2,'binwidth', 0.1)


% 11-06

names = {'11-1' '11-2' '12-1' '21-1' '21-2' '22-1' '23-1' '23-2' '31-1' '31-2'...
    '33-1' '33-2' '33-3' '34-1' '34-2' '41-1' '41-2' '41-3' '42-1' '42-2' '43-1'...
    '43-2' '43-3' '44-1' '44-2' '44-3' '51-1' '51-2' '51-3' '51-4' '52-1' '61-1'...
    '61-2' '61-3' '62-1' '63-1' '63-2' '63-3' '71-1' '71-2' '72-1' '81-1' '81-2'...
    '82-2' '82-1' '83-1' '83-2' '91-1' '91-2' '92-1' '92-2' '92-3'};

imgnames = {'11' '11' '12' '21' '21' '22' '23' '23' '31' '31'...
    '33' '33' '33' '34' '34' '41' '41' '41' '42' '42' '43'...
    '43' '43' '44' '44' '44' '51' '51' '51' '51' '52' '61'...
    '61' '61' '62' '63' '63' '63' '71' '71' '72' '81' '81'...
    '82' '82' '83' '83' '91' '91' '92' '92' '92'};



for k = 1:length(names)
[VarName, Area, StdDev, XM, YM, Slice] = importfile100620(['D:\Data\Raw\20.06.11\11-06-20_DepthoCouplé_ToutCollé_' names{k} '.txt']);



ptrsmall = find(Area < 50);

ptrbig = find(Area > 200);


[~, indfocussmall] = max(StdDev(ptrsmall));

focussmall = Slice(ptrsmall(indfocussmall));


[~, indfocusbig] = max(StdDev(ptrbig));

focusbig = Slice(ptrbig(indfocusbig));


if IMSHOW
Xb = XM(ptrbig(indfocusbig));
Yb = YM(ptrbig(indfocusbig));

Xs = XM(ptrsmall(indfocussmall));
Ys = YM(ptrsmall(indfocussmall));


Ib = imadjust(imread(['D:\Data\Raw\20.06.11\11-06-20_DepthoCouplé_ToutCollé_' imgnames{k} '.tif'],focusbig));
Is = imadjust(imread(['D:\Data\Raw\20.06.11\11-06-20_DepthoCouplé_ToutCollé_' imgnames{k} '.tif'],focussmall));

figure
subplot(121)
hold on
imshow(Ib)
plot(Xb,Yb,'r*')
subplot(122)
hold on
imshow(Is)
plot(Xs,Ys,'r*')

end


focusdist3(k) = (abs(focusbig-focussmall)*stepnm + (Bigradius-Smallradius))/1000;

end



focusdist3;

figure
histogram(focusdist3,'binwidth', 0.1)


figure
histogram([focusdist2 focusdist3],'binwidth', 0.15)
title(['Mean = ' num2str(mean([focusdist2 focusdist3])) 'µm - Std = ' num2str(std([focusdist2 focusdist3]))])
    
