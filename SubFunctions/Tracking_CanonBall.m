Questions={'Chemin Dossier Parent :','Chemin Dossier Kymographes :','Taille Objet [microns] :','Champ Magnétique [mT]','Gradient [mT/d]'};
definput = {'X','Y','2','35','40'};
Infos=inputdlg(Questions, 'Données Analyse',[1 30; 1 30; 1 30; 1 30; 1 30] ,definput);
j=1;
Dossier=[Infos{1} '\' num2str(j)];
while exist([Infos{1} '\' num2str(j)])
Dossier=[Infos{1} '\' num2str(j)];
Chemin=dir(Dossier);
Chemin = Chemin(~ismember({Chemin.name},{'.','..'}));
Kymograph=zeros(length(Chemin),2048);
for i=1:length(Chemin)
    Path=fullfile(Dossier,Chemin(i).name);
    Images=imread(Path,'tif');
    Tiret=max(strfind(Chemin(i).name,'_'));
    Ordre=str2double(Chemin(i).name(Tiret+1:end-4));
    Kymograph(Ordre,:)=min(Images);
    clear Images Path
end
figure(j)
imshow(Kymograph/max(max(Kymograph)))
xlabel('Position de la bille [micron]', 'FontSize',15)
ylabel('Axe temporel (20fps)','FontSize',15)
title(['Kymograph' num2str(j) ' bille ' Infos{3} 'um ; B=' Infos{4} 'mT ; Grad=' Infos{5} 'mT/d'],'FontSize',20)
saveas(figure(j),fullfile(Infos{2},['Kymograph' num2str(j) '.tif']),'tiffn')
j=j+1;
end
