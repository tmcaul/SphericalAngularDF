function [psi_pats_experiment,psi_pats_Rot_experiment]=Spherical_Phases(ExpPats,InputUser,RTI,MapData,PC,cs,Sphere)

for experiment=1:size(ExpPats,3)

%cd(InputUser.resultsfolder)
homefol=pwd;
expfolder=['Pattern_',num2str(experiment)];
mkdir(expfolder)
cd(expfolder)
resultsfolder=pwd;
    
cs;
h = Miller({1,0,0},{1,1,0},{1,3,1},{1,1,1},cs);%{2,0,0},{1,1,2},{1,2,3},cs);
hnames = {'{1,0,0}','{1,1,0}','{1,3,1}','{1,1,1}'};%,'{2,0,0}','{1,1,2}','{1,2,3}'};

Sphere.h=h;

PC_X=PC(1);
PC_Y=PC(2);
PC_Z=PC(3);

%%
pat=ExpPats(:,:,experiment);

cd(resultsfolder)
figure
imagesc(pat)
colormap('gray')
pbaspect([1,1,1])
print(['Pat'],'-dpng','-r300')
cd(homefol)

%% Get the correct orientation for MTEX
new_eulers=conv_G_to_EA(Sphere.rottoplot);
ori = orientation.byEuler(-new_eulers(1),-new_eulers(2),-new_eulers(3),'ZYZ',cs);
%ori = orientation.byEuler(-phi1,-PHI,-phi2,'ZYZ',cs);
ori0=orientation.byEuler(0,0,0,cs);

%% Look at  pattern
det = detector(RTI.screensize,RTI.screensize,PC_Z,[1-PC_Y,1-PC_X]);
% define a cut off function for the detector
mask = det.S2CutOffMask(0);
maskHarm = S2FunHarmonic.quadrature(mask);
correction = max(maskHarm.radon,0.1);

%approximate the measured pattern by harmonic function
pHarm = det.pattern2Fun(pat,'bandwidth',256,'quadrature','delta',0.1);
%get baseline pHarm (zeros)
pHarm_z = det.pattern2Fun(ones(size(pat)),'bandwidth',256,'quadrature','delta',0.1);

%get radon transform
RadonPHarm = conv(pHarm,Sphere.profileHarm)./ correction;
RadonPHarm.bandwidth = 64;

% Peak detection -> Figure 5b
peaks = S2PeakDetection(RadonPHarm,det,15,-0.11);

figure
plotx2east
plot(RadonPHarm,'pcolor','resolution',0.25*degree,'upper')
colormap('gray')
annotate(peaks,'markerSize',10,'antipodal','MarkerEdgeColor','c','MarkerFaceColor','none')

cd(resultsfolder)
print(['Radon_Peaks'],'-dpng','-r300')
cd(homefol)

%% Re-project pattern onto sphere and highlight strongest bands
plot(pHarm,'pcolor','resolution',0.25*degree,'upper')
colormap('gray')
circle(peaks,'linecolor',[0.35,0.35,0.35])

cd(resultsfolder)
print(['Reprojected_StrongestBands'],'-dpng','-r300')
cd(homefol)

%% Plot empty pattern
plot(pHarm_z,'pcolor','resolution',0.25*degree,'upper')
colormap('gray')

cd(resultsfolder)
print(['EmptyPattern'],'-dpng','-r300')
cd(homefol)


%% Rotate the projected pattern to the correct place on the ref [0,0,0] sphere
pHarm_rot=rotate(pHarm,inv(ori));
RadonPHarm_rot = conv(pHarm_rot,Sphere.profileHarm)./ correction;
RadonPHarm_rot.bandwidth = 64;
peaks_rot = S2PeakDetection(RadonPHarm_rot,det,15,-0.11);

figure
plot(pHarm_rot,'pcolor','resolution',0.25*degree,'upper')
circle(peaks_rot,'linecolor',[0.35,0.35,0.35])
colormap('gray')

cd(resultsfolder)
print(['Reprojected_OnRefSphere'],'-dpng','-r300');
cd(homefol)

%% Do the same for the empty pattern
pHarm_rot_z=rotate(pHarm_z,inv(ori));
figure
plot(pHarm_rot_z,'pcolor','resolution',0.25*degree,'upper')
colormap('gray')

cd(resultsfolder)
print(['EmptyPattern_Rotated'],'-dpng','-r300')
cd(homefol)

%% Highlight some specific plane normals

colors = cbrewer('qual','Set1',length(h));

figure
plot(RadonPHarm,'pcolor','resolution',0.25*degree,'upper','complete')
colormap('gray')

hold on
for i = 1:length(h)
  plot(ori*h(i).symmetrise,'MarkerSize',25,'MarkerFaceColor','none','MarkerEdgeColor',colors(i,:))
  plot(ori*h(i).symmetrise,'MarkerSize',26,'MarkerFaceColor','none','MarkerEdgeColor',colors(i,:))
  plot(ori*h(i).symmetrise,'MarkerSize',27,'MarkerFaceColor','none','MarkerEdgeColor',colors(i,:))
  plot(ori*h(i).symmetrise,'MarkerSize',23,'MarkerFaceColor','none','MarkerEdgeColor',colors(i,:))
  plot(ori*h(i).symmetrise,'MarkerSize',24,'MarkerFaceColor','none','MarkerEdgeColor',colors(i,:))
end
title(num2str(experiment))
hold off

cd(resultsfolder)
print(['ReprojectedRadon_SpecificBands'],'-dpng','-r300')
cd(homefol)

%%

figure
plot(pHarm,'pcolor','resolution',0.25*degree,'upper')
colormap('gray')
hold on
for i = 1:length(h)
    l{i}=circle(ori*h(i).symmetrise,'linecolor',colors(i,:),'LineWidth',1.5);
end

hold off
cd(resultsfolder)
print(['Reprojected_SpecificBands_',num2str(experiment)],'-dpng','-r300')
cd(homefol)

%% Highlight some specific plane normals for rotated

figure
plot(RadonPHarm_rot,'pcolor','resolution',0.25*degree,'upper','complete')
colormap('gray')

hold on
for i = 1:length(h)
  plot(h(i).symmetrise,'MarkerSize',25,'MarkerFaceColor','none','MarkerEdgeColor',colors(i,:))
  plot(h(i).symmetrise,'MarkerSize',26,'MarkerFaceColor','none','MarkerEdgeColor',colors(i,:))
  plot(h(i).symmetrise,'MarkerSize',27,'MarkerFaceColor','none','MarkerEdgeColor',colors(i,:))
  plot(h(i).symmetrise,'MarkerSize',23,'MarkerFaceColor','none','MarkerEdgeColor',colors(i,:))
  plot(h(i).symmetrise,'MarkerSize',24,'MarkerFaceColor','none','MarkerEdgeColor',colors(i,:))
end
title(num2str(experiment))
hold off

cd(resultsfolder)
print(['ReprojectedRadon_rot_SpecificBands'],'-dpng','-r300')
cd(homefol)

%%
figure
plot(pHarm_rot,'pcolor','resolution',0.25*degree,'upper')
colormap('gray')
hold on
for i = 1:length(h)
    circle(h(i).symmetrise,'linecolor',colors(i,:),'LineWidth',1.5)
end
hold off
cd(resultsfolder)
print(['Reprojected_rot_SpecificBands'],'-dpng','-r300')
cd(homefol)

%% Integrate along an unrotated empty sphere
% for i = 1:length(h)
%   [~,psi_pat_Baseline{i}] = pHarm_z.symmetrise(h(i));
%   psi_pat_Baseline{i}.A(1) = 0;
% end

%% Integrate along a rotated empty sphere
% for i = 1:length(h)
%   [~,psi_pat_Baseline_rot{i}] = pHarm_rot_z.symmetrise(h(i));
%   psi_pat_Baseline_rot{i}.A(1) = 0;
% end

%% Integrate along rotated normals on unrotated sphere

%Add up peak intensities for all bands of a given type
for i = 1:length(h)
  g=h(i).symmetrise;
  for j = 1:length(g(i).symmetrise)  
  [~,psi_pat{j}] = pHarm.symmetrise(ori*g(j));
  psi_pat{j}.A(1) = 0;
  end
  
  psi_pats_experiment{experiment,i}=psi_pat;
  clear psi_pat;
end

%% Integrate along unrotated normals on rotated sphere
for i = 1:length(h)
  g=h(i).symmetrise;
  for j=1:length(g(i).symmetrise)
  [~,psi_pat_Rot{j}] = pHarm_rot.symmetrise(g(j));
  psi_pat_Rot{j}.A(1) = 0;
  end
  
  psi_pats_Rot_experiment{experiment,i}=psi_pat_Rot;
  clear psi_pat_Rot;
end

%% Integrate along unrotated normals on unrotated master pattern
%only need to do for one because they're all the same...
%domain=domain.*180/pi;
close all
end