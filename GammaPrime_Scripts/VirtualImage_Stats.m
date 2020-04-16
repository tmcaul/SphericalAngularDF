VirtualImage

%% PDF Histograms

close all
for phase=1:4;
    data=reshape(VImage(:,:,phase),1,[]);
    h=histogram(data,40,'EdgeColor','None','EdgeAlpha',0,'Normalization','pdf','FaceColor',colors{phase}(end-100,:));
    hold on
    
    dist=fitdist(data','Normal');
    
    x=[h.BinLimits(1)-2*h.BinWidth,h.BinEdges,h.BinLimits(2)+2*h.BinWidth]+h.BinWidth/2;
    
    %plot(x,pdf(dist,x),'Color',colors{phase}(end,:))
    
end

ylim([0,200])
xlim([0,0.4])

pbaspect([5,2,1])

legend({'{111}','{110}','{100}','{131}'})


%% Get error in contrast value
% Load the ML data
InputUser.resultsfolder=[InputUser.ResultsDir,'/Results_',InputUser.HDF5_file(1:end-3)];
cd(InputUser.resultsfolder);
ML_extraction;
AvPat_extraction;
AE=Autoencoder_extraction('C:\Users\tpm416\Documents\GammaPrime_Autoencoders\Matlab_Autoencoder\AE_results_270320_kernel3_moresparsity');

nmfmap=flipud(Av.NMF_Map);
pcamap=flipud(Av.PCA_Map);
aemap=flipud(Av.AE_Map);

mlmaps={pcamap,nmfmap,aemap};

%final index is spatial
pcascores=zeros(size(ML.scores_PCA(:,:,:,1)));
for i=1:5
    pcascores(:,:,i)=flipud(ML.scores_PCA(:,:,i,3));
end

pca_separator=sum(pcascores(:,:,1:3),3);
nmf_separator=flipud(ML.scores_NMF(:,:,3,3));
ae_separator=reshape(AE.scores_AE(:,:,1)',72,95);

%%
qmax=99;

contrast=zeros(qmax,3,4);
error_contrast=zeros(qmax,3,4);

for q=1:qmax

quantileval=q/100;
    
pca_quantile=quantile(pca_separator(:),quantileval);
nmf_quantile=quantile(nmf_separator(:),quantileval);
ae_quantile=quantile(ae_separator(:),quantileval);

pcamap=pca_separator<pca_quantile;
nmfmap=nmf_separator<nmf_quantile;
aemap=ae_separator<ae_quantile;

mlmaps={pcamap,nmfmap,aemap};

means=zeros(3,2,4); %three clustering approaches, two phases, four bands
std_errors=zeros(3,2,4); % as above

%1 is gp, 0 is gamma so output is 1 gamma, 2 gp
for cluster=1:3
    for phase=1:2
        
        for band=1:4
            bandimage=VImage(:,:,band);
            vals=bandimage(mlmaps{cluster}==phase-1);
            
            means(cluster,phase,band)=mean(vals);
            std_errors(cluster,phase,band)=std(vals)./sqrt(length(vals));
            
            %savedvals{cluster,phase,band}=vals;
            
        end
    end
end

%contrast factor
for cluster=1:3
    for band=1:4
        mean_diff=(means(cluster,2,band)-means(cluster,1,band)); %gp = g
        error_diff=sqrt(std_errors(cluster,2,band).^2+std_errors(cluster,1,band).^2);
        contrast(q,cluster,band)=mean_diff/means(cluster,1,band);% (gp - g)/g
        error_contrast(q,cluster,band)=sqrt((error_diff/mean_diff).^2+(std_errors(cluster,1,band)/means(cluster,1,band)).^2);
    end
end

end

%% Simulation contrast factors

for band=1:4
    for phasepair=1:2
        contrast_sim(band,phasepair)=(braggsum(phasepair*2,band)-braggsum(phasepair*2-1,band))/braggsum(phasepair*2,band);
    end
end

%%
ClusterCol=[222,71,58;95,173,86;236,132,19]./255;
ind=[0.01:0.01:0.99];

ylims=[-0.1,0.05;-0.08,0.04;-0.15,0;-0.2,0.2];


for band=1:4
    subplot(2,2,band)
    for cluster=1:3
        %errorbar([0.001:0.001:0.999],contrast(:,cluster,band),error_contrast(:,cluster,band),'o','Color','k','MarkerFaceColor',ClusterCol(cluster,:),'MarkerSize',10)
        
        hold on
        ub=contrast(:,cluster,band)+error_contrast(:,cluster,band);
        lb=contrast(:,cluster,band)-error_contrast(:,cluster,band);
        
        %plot(ind,lb,'Color',ClusterCol(cluster,:))
        %plot(ind,ub,'Color',ClusterCol(cluster,:))
        
        ind2=[ind,fliplr(ind)];
        inbetween=[lb',fliplr(ub')];
        
        keep = ~isnan(ind2) & ~isnan(inbetween);
        ind2=ind2(keep);
        inbetween=inbetween(keep);
        
        fill(ind2,inbetween,ClusterCol(cluster,:),'FaceAlpha',0.15)
        
        
        xlim([0.2,0.8])
        ylim(ylims(band,:))
        
    end
    
    for cluster=1:3
        ub=contrast(:,cluster,band)+error_contrast(:,cluster,band);
        lb=contrast(:,cluster,band)-error_contrast(:,cluster,band);
        
        plot(ind,lb,'Color',ClusterCol(cluster,:),'LineWidth',2)
        plot(ind,ub,'Color',ClusterCol(cluster,:),'LineWidth',2)
        %plot(ind,contrast(:,cluster,band),'LineWidth',1,'Color',ClusterCol(cluster,:))
    end
    
    plot(ind,repmat(contrast_sim(cluster,1),size(ind,2)),'k--','LineWidth',1.5)
    %plot(ind,repmat(0,size(ind,2)),'k-')
    %plot(ind,repmat(contrast_sim(cluster,2),size(ind,2)),'b--')
    
    %ax=gca;
    %ax.YScale='log';
%     ax.XAxisLocation='origin';
    
    title(num2str(band))
    %xlim([0,4])
    
    y_range=ylim;
    y_midpoint=mean(y_range);
    y_range2=[y_midpoint-0.75*(y_range(2)-y_range(1)),y_midpoint+0.75*(y_range(2)-y_range(1))];
    ylim(y_range2)

end

%%
ylims=[-0.1,0.05;-0.1,0.05;-0.1,0.05;-0.1,0.05];
ind=[0:1:5];
for band=1:4
    subplot(2,2,band)
    hold on
    plot(ind,repmat(contrast_sim(cluster,1),size(ind,2)),'k--','LineWidth',1.5)
    plot(ind,repmat(0,size(ind,2)),'k-')
    plot(ind,repmat(contrast_sim(cluster,2),size(ind,2)),'b--')
    
    title(num2str(band))
    
    
    for cluster=1:3
        errorbar(cluster,contrast(37,cluster,band),error_contrast(37,cluster,band),'o','Color',ClusterCol(cluster,:),'LineWidth',1.5,'MarkerEdgeColor',ClusterCol(cluster,:),'MarkerFaceColor',ClusterCol(cluster,:),'MarkerSize',10)
    end
    xlim([0,4])
    
%     y_range=ylim;
%     y_midpoint=mean(y_range);
%     y_range2=[y_midpoint-2*(y_range(2)-y_range(1)),y_midpoint+2*(y_range(2)-y_range(1))];
%     ylim(y_range2)
    
   
    ylim([ylims(band,:)])
    
end





