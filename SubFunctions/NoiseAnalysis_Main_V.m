datafolder = 'D:\Data\Raw\20.07.03_BeadsTest';
depthoDirPath = 'D:\Data\Raw\';
depthoname = '03-07-20_DepthographM270';

scale=1/15.8; % 100X
% scale=1/9.7; % 63X
% En µm par pixel


%% M270
Diam = 0;

isImageTriplet = 0;
names = {'M270Movies_30ms_Exp_1-1' 'M270Movies_30ms_Exp_1-2' 'M270Movies_30ms_Exp_2-1' ...
    'M270Movies_30ms_Exp_2-2' 'M270Movies_30ms_Exp_4' 'M270Movies_30ms_Exp_5-1' 'M270Movies_30ms_Exp_5-2'}



   PLOTVID = 1;
NoiseAnalysis_NoiseFromBeadsChains_V(datafolder,names,depthoDirPath,depthoname,isImageTriplet,scale,Diam,PLOTVID)



%% MyOne

Diam = 1093;

isImageTriplet = 0;
names = {'MyOneMovies_30msExp_1-1' 'MyOneMovies_30msExp_2-2' 'MyOneMovies_30msExp_3' 'MyOneMovies_30msExp_4' 'MyOneMovies_30msExp_5'}

   PLOTVID = 1;
NoiseAnalysis_NoiseFromBeadsChains_V(datafolder,names,depthoDirPath,depthoname,isImageTriplet,scale,Diam,PLOTVID)

