
clc
format long
clear all
more off


%ouverture de fichiers dans lesquels sont rentrés les resultats
resultat_e = fopen('resultat_e.txt','a');
resultat_n= fopen('resultat_n.txt','a');
resultat_h = fopen('resultat_h.txt','a');

%recherche de la liste des fichiers doris du dossier
ls *res
%affichage des fichiers dans une liste l_fichiers
!ls *res > l_fichiers
!cat l_fichiers
%matrice de la liste des fichier
fichiers = textread('l_fichiers','%s');

taille = size(fichiers);

%traitement de tous les fichiers doris du dossier et remplissage des fichiers resultats
for i = 1:taille

[X_DE,X_DN,X_DH]=calcul_coeft(cell2mat(fichiers(i)),resultat_e, resultat_n, resultat_h)



end


%fermeture des fichiers
statut = fclose(resultat_e);
statut = fclose(resultat_n);
statut = fclose(resultat_h);

