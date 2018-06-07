(* ::Package:: *)

function[MaE,MaN,MaH,lieu, nom]=lecture(name)

%% permet de lire un fichier texte de mesures et de stocker les donnees dans des matrices correspondantes aux composantes E, N, H
%% name=nom du fichier sous forme de chaine de caracteres
%% lieu = code de location de la station
%% nom = identifiant de la station

[F]=fopen(name,'rt') %% ouverture du fichier
MaE=[]; %% initialisation de la matrice contenant les informations du fichier pour la composante E
%% chaque ligne de la matrice correspond \[AGrave] l'ann\[EAcute]e, au r\[EAcute]sidu et \[AGrave] l'\[EAcute]cart-type
MaN=[];
MaH=[];

    while ~feof(F) %% tant qu'on n'a pas atteint la fin du fichier
	
	
		
		
		line=fgetl(F); %% lecture de la premiere ligne
		splitline=strsplit(line);
		lieu = cell2mat(splitline(1));
		nom = cell2mat(splitline(2));
        	t = str2num(cell2mat(splitline(3)));
		resE=str2num(cell2mat(splitline(4))); %% r\[CapitalATilde]\[Copyright]sidu
		resN=str2num(cell2mat(splitline(5))); %% r\[CapitalATilde]\[Copyright]sidu
		resH=str2num(cell2mat(splitline(6))); %% r\[CapitalATilde]\[Copyright]sidu
		sigE=str2num(cell2mat(splitline(7))); %% sigma
		sigN=str2num(cell2mat(splitline(8))); %% sigma
		sigH=str2num(cell2mat(splitline(9))); %% sigma
		
		
		vE=[t,resE,sigE];
		MaE=[MaE ; vE];

		vH=[t,resH,sigH];
		MaH=[MaH ; vH];

		vN=[t,resN,sigN];
		MaN=[MaN ; vN];
	
    end


%% on ne consid\[EGrave]re pas les r\[EAcute]sidus ayant un \[EAcute]cart-type nul
ind = find(MaE (:,3) ~= 0);
MaE = MaE(ind,:);

ind = find(MaN (:,3) ~= 0);
MaN = MaN(ind,:);

ind = find(MaH (:,3) ~= 0);
MaH = MaH(ind,:);


status=fclose(F)

end
