function SavePaperMyo

date = datestr(now,'yy-mm-dd');

name = [date '-MyoPaperSave'];
cd('C:\Users\Valentin\Dropbox\TheseValentin')
mkdir(['.\SauvegardeFigPaperMyosin\' name])

cd('C:\Users\Valentin\Dropbox\TheseValentin\Papiers\MatlabFig')

copyfile('Make_Fig1.m',['C:\Users\Valentin\Dropbox\TheseValentin\SauvegardeFigPaperMyosin\' name '\'])
copyfile('Make_Fig2.m',['C:\Users\Valentin\Dropbox\TheseValentin\SauvegardeFigPaperMyosin\' name '\'])
copyfile('Make_Fig3.m',['C:\Users\Valentin\Dropbox\TheseValentin\SauvegardeFigPaperMyosin\' name '\'])
copyfile('Make_Fig4.m',['C:\Users\Valentin\Dropbox\TheseValentin\SauvegardeFigPaperMyosin\' name '\'])

pause(0.5)

fprintf('Done\n')
end