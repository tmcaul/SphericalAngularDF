folders={'Local_NMF','Local'};

i=1;
j=1;
fol=folders{i};
fol=[fol,num2str(j)];
cd('E:\Tom\GammaPrime_Data\V208C\yprime3_analysis')
cd(fol)
Av.NMF_Pats=h5read('LabelPatterns.h5','/NMF');
Av.NMF_Map=h5read('LabelPatterns.h5','/NMF_region');
Av.NMF_Pats=rot90((Av.NMF_Pats),1);
Av.NMF_Map=rot90(Av.NMF_Map,1);

i=2;
j=1;
fol=folders{i};
fol=[fol,num2str(j)];
cd('E:\Tom\GammaPrime_Data\V208C\yprime3_analysis')
cd(fol)
Av.PCA_Pats=h5read('LabelPatterns.h5','/PCA');
Av.PCA_Map=h5read('LabelPatterns.h5','/PCA_region');
Av.PCA_Pats=rot90((Av.PCA_Pats),1);
Av.PCA_Map=rot90(Av.PCA_Map,1);

clear i j fol folders
cd(homefol);

