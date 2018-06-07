function[l_gnss,c, existe]=trouve(fdoris)

%%retourne la liste des fichiers gnss ayant le même code de location que le fichier doris en entrée, c est la taille de cette liste (vaut 5 maximum)
%%existe est un booléen qui vaut true s'il existe un tel fichier gnss

%extraction du code de location commun gnss et doris
indice=strfind(fdoris,'_');
code_locd = fdoris(indice+1:indice+5);

%recherche du dossier gnss
cd ..
cd gnss

%recherche des fichiers gnss
ls *res;
!ls *res > l_fichiersg;
!cat l_fichiersg;
fichiersg = textread('l_fichiersg','%s');

taille = size(fichiersg);
 
%creation de la liste des fichiers gnss ayant le même code de location
l_gnss = [];
c = 0;
existe = false;
for i = 1:taille
	fichierg = cell2mat(fichiersg(i));
	indice=strfind(fichierg,'_');
	code_locg = fichierg(indice+1:indice+5);
	
	%verification si mÃªme code_loc
	if code_locd == code_locg
		[MaE,MaN,MaH,lieu, nom]=lecture(fichierg)
		if size(MaE)>1000
		%on ne prend que les fichiers avec un temps d'observation assez long
			fg = [fichierg];
			l_gnss = [l_gnss;fichiersg(i)];
			c = c+1;
			existe=true;
		end

		if c==5
		%on ne prend que 5 fichiers gnss maximum
			break
		end
		
	end
end


cd .. 
cd doris

end

