%% Load the ML data
%cd(InputUser.resultsfolder);
ML_extraction;
%row,col,comp,spatial

%% Histograms of scores for different ML regimes
cd(InputUser.resultsfolder);
mkdir('BinningHistograms')
cd('BinningHistograms')

histfol=pwd;
mkdir('NMF')
cd('NMF')

for factor=1:3;
    clear freq bins
    figure
    for i=1:5 %loop over spatial localisations
    d=ML.scores_NMF(:,:,factor,i);
    
    m_nmf(factor,i)=mean(d(:));%calculate mean scores
    s_nmf(factor,i)=std(d(:));%calculate standard deviation of scores
    
    [freq(i,:),bins(i,:)]=hist(d(:),30);
    plot(bins(i,:),freq(i,:),'LineWidth',1.5)
    hold on
    end
    hold off
    legend('r = 1', 'r = 2', 'r = 3', 'r = 4', 'r = 5')
    xlabel('NMF score')
    ylabel('Frequency density')
    print(['NMF_SpatialHistograms_',num2str(factor)],'-dpng','-r300')
    
end

figure
for factor=1:3
    plot(s_nmf(factor,:),'LineWidth',1.5)
    hold on
end
ylabel('Standard deviation of scores')
xlabel('Spatial kernel, r')
legend('Factor 1','Factor 2','Factor 3','Location','south','Orientation','horizontal','FontSize',10)
print('Stdevs_NMF','-dpng','-r300')

figure
for factor=1:3
    plot(m_nmf(factor,:),'LineWidth',1.5)
    hold on
end
ylabel('Mean score')
xlabel('Spatial kernel, r')
legend('Factor 1','Factor 2','Factor 3','Location','south','Orientation','horizontal','FontSize',10)
print('Means_NMF','-dpng','-r300')

%%
cd(histfol)
mkdir('PCA')
cd('PCA')

for factor=1:5;
    clear freq bins
    
    figure
    for i=1:5 %loop over spatial localisations
    d=ML.scores_PCA(:,:,factor,i);
    d2=d./ML.explained_PCA(factor,i); %normalise by eigenvalue
    
    m(factor,i)=mean(d(:));%calculate mean scores
    s(factor,i)=std(d(:));%calculate standard deviation of scores
    
    dist=fitdist(d2(:),'Normal');
    m_normal(factor,i)=dist.mu;
    s_normal(factor,i)=dist.sigma;
    
    [freq(i,:),bins(i,:)]=hist(d2(:),30);
    plot(bins(i,:),freq(i,:),'LineWidth',1.5)
    hold on
    end
    hold off
    legend('r = 1', 'r = 2', 'r = 3', 'r = 4', 'r = 5')
    xlabel(['PC ',num2str(factor), ' score'])
    xlim([-30,30])
    ylabel('Frequency density')
    print(['PCA_SpatialHistograms_',num2str(factor)],'-dpng','-r300')
     
end


% Plot means and standard deviations
% fig = figure;
% left_color = [0 0 0];
% right_color = [0.9 0 0];
% set(fig,'defaultAxesColorOrder',[left_color; right_color]);
% 
% yyaxis left
% %s=std(freq,0,2);
% plot(s,'LineWidth',1.5)
% ylabel('Standard deviation of scores')
% 
% yyaxis right
% %m=mean(freq,2);
% plot(m,'LineWidth',1.5)
% ylabel('Mean score')
% 
% xlabel('r')
% xticklabels([1,2,3,4,5])
% xlim([1,5])

figure
for factor=1:5
    plot(s(factor,:),'LineWidth',1.5)
    hold on
end
ylabel('Standard deviation of scores')
xlabel('Spatial kernel, r')
legend('PC 1','PC 2','PC 3','PC 4','PC 5','Location','south','Orientation','horizontal','FontSize',10)
print('Stdevs','-dpng','-r300')

figure
for factor=1:5
    plot(s_normal(factor,:),'LineWidth',1.5)
    hold on
end
ylabel('Sigma of fitted normal')
xlabel('Spatial kernel, r')
legend('PC 1','PC 2','PC 3','PC 4','PC 5','Location','NorthWest','Orientation','vertical','FontSize',10)
print('Sigmas','-dpng','-r300')

%print(['PCA_SpatialStatistics_',num2str(factor)],'-dpng','-r300')

cd(histfol)

%% Compare singular values and eigenvalues for differing spatial regimes (PCA)
close all

cd(histfol)
cd('PCA')

%Scree plot (eigenvalues)
figure
for j=1:5
    plot(100*ML.explainedprop_PCA(:,j),'LineWidth',1.5)
    hold on
end
legend('r = 1','r = 2','r = 3','r = 4','r = 5');
%set(gca,'YScale','log')
xticklabels([1,2,3,4,5])
%yticklabels([10,100])
xlabel('Principal component')
ylabel('Explained variance / %')
print('Scree','-dpng','-r300')


%Singular values
figure
for j=1:5
    plot(ML.coefflengths_PCA(:,j),'LineWidth',1.5)
    hold on
end
legend('r = 1','r = 2','r = 3','r = 4','r = 5');
%set(gca,'YScale','log')
xticklabels([1,2,3,4,5])
xlabel('Principal component')
ylabel('Component L2 norm')
print('SingularValues','-dpng','-r300')

%total variance explained
figure
plot(ML.totalvar_PCA,'LineWidth',1.5)
set(gca,'YScale','log')
ylim([0.9,1.01])
xticklabels([1,2,3,4,5])
xlabel('Spatial kernel, r')
ylabel('Total variance explained')

%% PDF investigations
% factor=2;i=1;
% d=ML.scores_NMF(:,:,factor,i);
% distribution=fitdist(d(:),'Lognormal');
% [h,p]=chi2gof(d(:),'CDF',distribution,'Alpha',0.01)

