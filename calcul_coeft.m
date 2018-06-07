function [X_DE,X_DN,X_DH] = calcul_coeft(fichierd, resultat_e, resultat_n, resultat_h)

%%permet de calculer les coefficients correspondants aux annuel, semi-annuel, 1ère draco jusqu'à la 7ème draco de mesures gnss et doris, ensemble
%%pour les trois composantes E,N,H par moindres carrés
%%entrées : 1) fichierd : fichier texte de mesures DORIS contenant l'identifiant de la station, les années, les résidus E,N et H en et les écarts-types respectifs en mm
%% 2) resultat_e : fichier texte dans lequel s'écrit l'identifiant de la station, le facteur unitaire de variance, X_DE, la matrice de covariance et la matrice de corrélation
%% De même pour resultat_n (avec X_DN) et resultat_h (avec X_DH)
%% sorties : X_DE, X_DN et X_DH sont des vecteurs contenant dans l'ordre les cofficients (cos puis sin) correspondants aux annuel, semi-annuel, 1ère draco à la 7ème draco

[l_gnss,c, existe] = trouve(fichierd) %%on cherche les fichiers gnss ayant le même code de location que le fichier doris en entrée
if existe == true %%s'il y a au moins un fichier gnss correspondant
[M_DEd,M_DNd,M_DHd, lieud, nomd] = lecture(fichierd);  %%on crée les matrices de données doris des composantes E,N,H

%%définition des pulsations de l'annuel, du semi-annuel et des 7 premières draco pour le gnss et le doris

Wa=2*pi;
Wsa=4*pi;


Wd1=2*pi*365.25/351.4;
Wd2=2*Wd1;
Wd3=3*Wd1;
Wd4=4*Wd1;
Wd5=5*Wd1;
Wd6=6*Wd1;
Wd7=7*Wd1;
W_GPS=[Wd1;Wd2;Wd3;Wd4;Wd5;Wd6;Wd7];



Wd8=2*pi*365.25/108;
Wd9=2*Wd8;
Wd10=3*Wd8;
Wd11=4*Wd8;
Wd12=5*Wd8;
Wd13=6*Wd8;
Wd14=7*Wd8;
W_DORIS= [Wa;Wsa;Wd8;Wd9;Wd10;Wd11;Wd12;Wd13;Wd14];
td = M_DEd(:,1);



%E
%creation d'une matrice des noms et matrice des lieux
%on remplit avec ceux de la staion doris
lieux = [lieud];
noms = [nomd];
%%on écrit le nom et le code de location de la station doris dans resultat_e
fprintf(resultat_e,'%s\t', lieud);
fprintf(resultat_e,'%s\n', nomd);

%lecuture de tout le fichier doris
SigmaDORISE=[M_DEd(:,3)]; %%vecteur des écarts-types
PDORISE=1./SigmaDORISE.^2;
Pi = PDORISE;

YDORISE = M_DEd(:,2);
YE = YDORISE; %%vecteur des observations

A_DORISE = []; %%partie doris de la matrice modèle

cd ..
cd gnss
NDORIS=size(M_DEd);
for k=1:NDORIS
   Wt_DORIS=td(k)*W_DORIS;
	ca = cos(td(k)*W_DORIS(1));
	sa = sin(td(k)*W_DORIS(1));
	csa = cos(td(k)*W_DORIS(2));
	ssa = sin(td(k)*W_DORIS(2));
	
   ME=[ca,sa,csa,ssa];
   
   for i=3:9
       cs=cos(Wt_DORIS(i));
       sn=sin(Wt_DORIS(i));
       ME=[ME,cs,sn];
   end
	for j=1:c*14
		ME = [ME,0];
	end

   A_DORISE=[A_DORISE;ME];

end


%on s'occupe ensuite des fichiers gnss
A_GPSE = [];

for m=1:c %%boucle sur les fichiers gnss ayant le même code de location que le fichier doris en entrée
	type(cell2mat(l_gnss(m)))
	[M_DEg,M_DNg,M_DHg, lieug, nomg] = lecture(cell2mat(l_gnss(m)));
	fprintf(resultat_e,'%s\t', lieug);
	fprintf(resultat_e,'%s\n', nomg);


	NGPS=size(M_DEg);

	an = [];
	Draco = [];

%crÃ©ation de matrices vides qui nous aideront à construire la matrice modèle des mesures gnss
		vide1 = zeros(NGPS(1),14*m);
		vide2 = zeros(NGPS(1),14*(c-m));



	
	for j=1:NGPS
	
	tg = M_DEg(:,1);
   	Wt_GPS=tg(j)*W_GPS;
	ca = cos(tg(j)*W_DORIS(1));
	sa = sin(tg(j)*W_DORIS(1));
	csa = cos(tg(j)*W_DORIS(2));
	ssa = sin(tg(j)*W_DORIS(2));
	lAN = [ca,sa,csa,ssa];

   	an=[an;lAN];


	
	ME = [];
	for k=1:7
       		cs=cos(Wt_GPS(k));
       		sn=sin(Wt_GPS(k));
       		ME=[ME,cs,sn];
	end
	Draco = [Draco;ME];
	
	end
	Af_GPS = [an,vide1,Draco,vide2];
	A_GPSE = [A_GPSE;Af_GPS];


	YGPSE = M_DEg(:,2);
	YE = [YE;YGPSE]; %%matrice des observations globales
	SigmaGPSE=[M_DEg(:,3)];
	PGPSE=1./SigmaGPSE.^2;
	Pi = [Pi;PGPSE];
	
end
PE = diag(Pi); %%matrice de poids globale
AE = [A_DORISE;A_GPSE]; %%matrice modèle
	

%calcul des coefficients par moindres carrÃ©es
X_DE=inv(AE'*PE*AE)*AE'*PE*YE;
VE=AE*X_DE-YE; %%vecteur des résidus

ecart_tE=(VE'*PE*VE)/(size(YE,1)-size(X_DE,1)); %%facteur unitaire de variance

var_XE=ecart_tE*inv(AE'*PE*AE); %%matrice de covariance

X_DN = [];
X_DH = [];


corr_xE = []; %%matrice de corrélation
for i=1:14*(c+1)+4
	corrE = [];
	for j=1:14*(c+1)+4
		correlation = var_XE(i,j)/sqrt(var_XE(i,i)*var_XE(j,j));
		corrE = [corrE,correlation];
		if correlation >= 0.5
			disp(['correlation : ' num2str(correlation)]);		
			
			
		end	
 	end
	corr_xE = [corr_xE;corrE];
end

%ecriture dans le fichier des résultats
fo = repmat('%11.8f ',1, 14*(c+1)+4);
fo = [fo,'\n'];


fprintf(resultat_e,'%f\n', ecart_tE);
fprintf(resultat_e,'\n');

fprintf(resultat_e,fo,X_DE);
fprintf(resultat_e,'\n');

fprintf(resultat_e,fo,var_XE);
fprintf(resultat_e,'\n');

fprintf(resultat_e,fo,corr_xE);
fprintf(resultat_e,'\n');
fprintf(resultat_e,'\n');
fprintf(resultat_e,'\n');









%N
%creation d'une matrice des noms et matrice des lieux
%on remplit avec ceux de la staion doris
lieux = [lieud];
noms = [nomd];
fprintf(resultat_n,'%s\t', lieud);
fprintf(resultat_n,'%s\n', nomd);

%lecuture de tout le fichier doris
SigmaDORISN=[M_DNd(:,3)];
PDORISN=1./SigmaDORISN.^2;
PiN = PDORISN;

YDORISN = M_DNd(:,2);
YN = YDORISN;

A_DORISN = [];

cd ..
cd gnss
NDORIS=size(M_DNd);
for k=1:NDORIS
   Wt_DORIS=td(k)*W_DORIS;
	ca = cos(td(k)*W_DORIS(1));
	sa = sin(td(k)*W_DORIS(1));
	csa = cos(td(k)*W_DORIS(2));
	ssa = sin(td(k)*W_DORIS(2));
	
   MN=[ca,sa,csa,ssa];
   
   for i=3:9
       cs=cos(Wt_DORIS(i));
       sn=sin(Wt_DORIS(i));
       MN=[MN,cs,sn];
   end
	for j=1:c*14
		MN = [MN,0];
	end

   A_DORISN=[A_DORISN;MN];

end


%on s'occupe ensuite des fichiers gnss
A_GPSN = [];

for m=1:c

	[M_DEg,M_DNg,M_DHg, lieug, nomg] = lecture(cell2mat(l_gnss(m)));
	fprintf(resultat_n,'%s\t', lieug);
	fprintf(resultat_n,'%s\n', nomg);


	NGPS=size(M_DNg);

	an = [];
	Draco = [];



		vide1 = zeros(NGPS(1),14*m);



		vide2 = zeros(NGPS(1),14*(c-m));



	
	for j=1:NGPS
	
	tg = M_DNg(:,1);
   	Wt_GPS=tg(j)*W_GPS;
	ca = cos(tg(j)*W_DORIS(1));
	sa = sin(tg(j)*W_DORIS(1));
	csa = cos(tg(j)*W_DORIS(2));
	ssa = sin(tg(j)*W_DORIS(2));
	lAN = [ca,sa,csa,ssa];

   	an=[an;lAN];


	MN = [];
	for k=1:7
       		cs=cos(Wt_GPS(k));
       		sn=sin(Wt_GPS(k));
       		MN=[MN,cs,sn];
	end
	Draco = [Draco;MN];
	
	
	end
	Af_GPSN = [an,vide1,Draco,vide2];
	A_GPSN = [A_GPSN;Af_GPSN];


	
	YGPSN = M_DNg(:,2);
	YN = [YN;YGPSN];
	SigmaGPSN=[M_DNg(:,3)];
	PGPSN=1./SigmaGPSN.^2;
	PiN = [PiN;PGPSN];
	
end
PN = diag(PiN);
AN = [A_DORISN;A_GPSN];
	

%calcul des moindres carrÃ©es
size(AN)
size(PN)
size(YN)
X_DN=inv(AN'*PN*AN)*AN'*PN*YN;

VN=AN*X_DN-YN;

ecart_tN=(VN'*PN*VN)/(size(YN,1)-size(X_DN,1));

var_XN=ecart_tN*inv(AN'*PN*AN);



corr_xN = [];
for i=1:14*(c+1)+4
	corrN = [];
	for j=1:14*(c+1)+4
		correlation = var_XN(i,j)/sqrt(var_XN(i,i)*var_XN(j,j));
		corrN = [corrN,correlation];
		if correlation >= 0.5
			disp(['correlation : ' num2str(correlation)]);		
			
			
		end	
 	end
	corr_xN = [corr_xN;corrN];
end

%ecriture dans le fichier


fprintf(resultat_n,'%f\n', ecart_tN);
fprintf(resultat_n,'\n');

fprintf(resultat_n,fo,X_DN);
fprintf(resultat_n,'\n');

fprintf(resultat_n,fo,var_XN);
fprintf(resultat_n,'\n');

fprintf(resultat_n,fo,corr_xN);
fprintf(resultat_n,'\n');
fprintf(resultat_n,'\n');
fprintf(resultat_n,'\n');











%H
%creation d'une matrice des noms et matrice des lieux
%on remplit avec ceux de la staion doris
lieux = [lieud];
noms = [nomd];
fprintf(resultat_h,'%s\t', lieud);
fprintf(resultat_h,'%s\n', nomd);

%lecuture de tout le fichier doris
SigmaDORISH=[M_DHd(:,3)];
PDORISH=1./SigmaDORISH.^2;
PiH = PDORISH;

YDORISH = M_DHd(:,2);
YH = YDORISH;

A_DORISH = [];

cd ..
cd gnss
NDORIS=size(M_DHd);
for k=1:NDORIS
   Wt_DORIS=td(k)*W_DORIS;
	ca = cos(td(k)*W_DORIS(1));
	sa = sin(td(k)*W_DORIS(1));
	csa = cos(td(k)*W_DORIS(2));
	ssa = sin(td(k)*W_DORIS(2));
	
   MH=[ca,sa,csa,ssa];
   
   for i=3:9
       cs=cos(Wt_DORIS(i));
       sn=sin(Wt_DORIS(i));
       MH=[MH,cs,sn];
   end
	for j=1:c*14
		MH = [MH,0];
	end

   A_DORISH=[A_DORISH;MH];

end


%on s'occupe ensuite des fichiers gnss
A_GPSH = [];

for m=1:c

	[M_DEg,M_DNg,M_DHg, lieug, nomg] = lecture(cell2mat(l_gnss(m)));
	fprintf(resultat_h,'%s\t', lieug);
	fprintf(resultat_h,'%s\n', nomg);


	NGPS=size(M_DHg);

	an = [];
	Draco = [];



		vide1 = zeros(NGPS(1),14*m);


		vide2 = zeros(NGPS(1),14*(c-m));



	
	for j=1:NGPS
	
	tg = M_DHg(:,1);
   	Wt_GPS=tg(j)*W_GPS;
	ca = cos(tg(j)*W_DORIS(1));
	sa = sin(tg(j)*W_DORIS(1));
	csa = cos(tg(j)*W_DORIS(2));
	ssa = sin(tg(j)*W_DORIS(2));
	lAN = [ca,sa,csa,ssa];

   	an=[an;lAN];


	
	MH = [];
	for k=1:7
       		cs=cos(Wt_GPS(k));
       		sn=sin(Wt_GPS(k));
       		MH=[MH,cs,sn];
	end
	Draco = [Draco;MH];
	
	end
	Af_GPSH = [an,vide1,Draco,vide2];
	A_GPSH = [A_GPSH;Af_GPSH];


	
	YGPSH = M_DHg(:,2);
	YH = [YH;YGPSH];
	SigmaGPSH=[M_DHg(:,3)];
	PGPSH=1./SigmaGPSH.^2;
	PiH = [PiH;PGPSH];
	
end
PH = diag(PiH);
AH = [A_DORISH;A_GPSH];
	

%calcul des moindres carrÃ©es
X_DH=inv(AH'*PH*AH)*AH'*PH*YH;
VH=AH*X_DH-YH;

ecart_tH=(VH'*PH*VH)/(size(YH,1)-size(X_DH,1));

var_XH=ecart_tH*inv(AH'*PH*AH);



corr_xH = [];
for i=1:14*(c+1)+4
	corrH = [];
	for j=1:14*(c+1)+4
		correlation = var_XH(i,j)/sqrt(var_XH(i,i)*var_XH(j,j));
		corrH = [corrH,correlation];
		if correlation >= 0.5
			disp(['correlation : ' num2str(correlation)]);		
			
			
		end	
 	end
	corr_xH = [corr_xH;corrH];
end

%ecriture dans le fichier


fprintf(resultat_h,'%f\n', ecart_tH);
fprintf(resultat_h,'\n');

fprintf(resultat_h,fo,X_DH);
fprintf(resultat_h,'\n');

fprintf(resultat_h,fo,var_XH);
fprintf(resultat_h,'\n');

fprintf(resultat_h,fo,corr_xH);
fprintf(resultat_h,'\n');
fprintf(resultat_h,'\n');
fprintf(resultat_h,'\n');

cd ..
cd doris


else %%s'il n'y a pas de fichier gnss ayant le même code de location que le fichier doris en entrée
X_DE = [];
X_DN = [];
X_DH = [];

end

end

