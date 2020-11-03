function SaveMatScript(str,EXIT)

cd('C:\Users\Valentin\Dropbox\TheseValentin')

date = datestr(now,'yy-mm-dd');

name = [date '_' str];

% date = '17-09-22'

mkdir(['.\SauvegardeMatlab\' name])

copyfile('Matlab',['.\SauvegardeMatlab\' name '\'])

cd('C:\Users\Valentin\Dropbox\TheseValentin\Matlab')

slocDir('.')

pause(1)


if EXIT
exit
end
end