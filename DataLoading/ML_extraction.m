folders={'Local_NMF','Local'};

ML.scores_NMF=zeros(72,95,3,5);
ML.coeffs_NMF=zeros(300,300,3,5);

i=1;
for j=1:5
    fol=folders{i};
    fol=[fol,num2str(j)];

    cd('E:\Tom\GammaPrime_Data\V208C\yprime3_analysis')
    cd(fol)
    coeff=h5read('scorescoeffs.h5','/coeffs');
    score=h5read('scorescoeffs.h5','/scores');
    
%     if mean(score(:))<0
%         coeff=-coeff;
%         score=-score;
%     end
    
    ML.scores_NMF(:,:,:,j)=rot90((score),1);
    ML.coeffs_NMF(:,:,:,j)=rot90((coeff),1);

end

%%
%row,col,component,spatial
ML.scores_PCA=zeros(72,95,5,5);
ML.coeffs_PCA=zeros(300,300,5,5);

i=2;
for j=1:5
    fol=folders{i};
    fol=[fol,num2str(j)];

    cd('E:\Tom\GammaPrime_Data\V208C\yprime3_analysis')
    cd(fol)
    coeff=h5read('scorescoeffs.h5','/coeffs');
    score=h5read('scorescoeffs.h5','/scores');
    
%     if mean(score(:))<0
%         coeff=-coeff;
%         score=-score;
%     end
    
    ML.scores_PCA(:,:,:,j)=rot90((score),1);
    ML.coeffs_PCA(:,:,:,j)=rot90((coeff),1);
    ML.explained_PCA(:,j)=h5read('scorescoeffs.h5','/eigenvalues');
    ML.totalvar_PCA(j)=h5read('scorescoeffs.h5','/total_var');
    
    ML.explainedprop_PCA(:,j)=ML.explained_PCA(:,j)./sum(ML.explained_PCA(:,j));
    ML.coefflengths_PCA(:,j)=sqrt(ML.explained_PCA(:,j)*6839); %equivalent to singular values
end

clear i j fol score coeff
cd(homefol);

