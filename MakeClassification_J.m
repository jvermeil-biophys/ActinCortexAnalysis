function MakeClassification_J(GeneralFolder)
% dates = {'dd-mm-yy'}; % Example dates
% dates = {'20-02-20'}; 
dates = {'20-02-20','10-03-20'}; 

%% Etat de la cellules

Migrante = {{} {}}; 
%  Migrante = {{'112'}} <=> 'La cellule du fichier M1_P1_C2 est migrante'

Etalee   = {{} {}};

Genee    = {{} {}}; 


    AllMigr = {};
    AllEtal = {};
    AllGene = {};
for kd = 1:length(dates)
    if ~isempty(Migrante{kd})
        for km = 1:length(Migrante{kd})            
            AllMigr = {AllMigr{:} [dates{kd} '_M' Migrante{kd}{km}(1) '_P' Migrante{kd}{km}(2) '_C' Migrante{kd}{km}(3:end)]};
        end
    end
     if ~isempty(Etalee{kd})
        for km = 1:length(Etalee{kd})            
            AllEtal = {AllEtal{:} [dates{kd} '_M' Etalee{kd}{km}(1) '_P' Etalee{kd}{km}(2) '_C' Etalee{kd}{km}(3:end)]};
        end
     end
     if ~isempty(Genee{kd})
        for km = 1:length(Genee{kd})            
            AllGene = {AllGene{:} [dates{kd} '_M' Genee{kd}{km}(1) '_P' Genee{kd}{km}(2) '_C' Genee{kd}{km}(3:end)]};
        end
    end
end


%% Nombre de billes dans la cellule

uneBille = {{} {'111', '121', '122', '123', '124'}};

deuxBille = {{} {'125'}};

troisBille = {{} {'122', '124'}};

quatreBille = {{} {}};

cinqBille = {{} {'123'}};

sixBille = {{} {}};

    AllUne = {};
    AllDeux = {};
    AllTrois = {};
    AllQuatre = {};
    AllCinq = {};
    AllSix = {};
for kd = 1:length(dates)
    if ~isempty(uneBille{kd})
        for km = 1:length(uneBille{kd})            
            AllUne = {AllUne{:} [dates{kd} '_M' uneBille{kd}{km}(1) '_P' uneBille{kd}{km}(2) '_C' uneBille{kd}{km}(3:end)]};
        end
    end

    if ~isempty(deuxBille{kd})
        for km = 1:length(deuxBille{kd})            
            AllDeux = {AllDeux{:} [dates{kd} '_M' deuxBille{kd}{km}(1) '_P' deuxBille{kd}{km}(2) '_C' deuxBille{kd}{km}(3:end)]};
        end
    end

    if ~isempty(troisBille{kd})
        for km = 1:length(troisBille{kd})            
            AllTrois = {AllTrois{:} [dates{kd} '_M' troisBille{kd}{km}(1) '_P' troisBille{kd}{km}(2) '_C' troisBille{kd}{km}(3:end)]};
        end
    end

    if ~isempty(quatreBille{kd})
        for km = 1:length(quatreBille{kd})            
            AllQuatre = {AllQuatre{:} [dates{kd} '_M' quatreBille{kd}{km}(1) '_P' quatreBille{kd}{km}(2) '_C' quatreBille{kd}{km}(3:end)]};
        end
    end

    if ~isempty(cinqBille{kd})
        for km = 1:length(cinqBille{kd})            
            AllCinq = {AllCinq{:} [dates{kd} '_M' cinqBille{kd}{km}(1) '_P' cinqBille{kd}{km}(2) '_C' cinqBille{kd}{km}(3:end)]};
        end
    end

    if ~isempty(sixBille{kd})
        for km = 1:length(sixBille{kd})            
            AllSix = {AllSix{:} [dates{kd} '_M' sixBille{kd}{km}(1) '_P' sixBille{kd}{km}(2) '_C' sixBille{kd}{km}(3:end)]};
        end
    end

end

%% Sauvegarde
mkdir([GeneralFolder 'Data_Joseph\MatFiles\V2D'])
save([GeneralFolder 'Data_Joseph\MatFiles\V2D\Classification.mat'],'AllMigr', 'AllEtal', 'AllGene', 'AllUne', 'AllDeux' ,'AllTrois' ,'AllQuatre' ,'AllCinq' ,'AllSix')

end