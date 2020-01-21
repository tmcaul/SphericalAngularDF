function [psi_full,psi_pats,psi_pats_Rot]=Spherical_Phases(InputUser,RTI,MapData,PC,crystals,Sphere)

for phase=1:length(InputUser.Phases)

%cd(InputUser.resultsfolder)
homefol=pwd;
expfolder=['Phase_',num2str(phase)];
mkdir(expfolder)
cd(expfolder)
resultsfolder=pwd;


    
cs=crystals{phase};
h = Miller({1,0,0},{1,1,0},{1,3,1},{1,1,1},cs);%{2,0,0},{1,1,2},{1,2,3},cs);
hnames = {'{1,0,0}','{1,1,0}','{1,3,1}','{1,1,1}'};%,'{2,0,0}','{1,1,2}','{1,2,3}'};
    
%% Define some interesting bands for summing
InputUser.Phase_Input=InputUser.Phases(phase);

PC_X=PC(1);
PC_Y=PC(2);
PC_Z=PC(3);

%% Generate a pattern

% Master patterns
[screen_int,facedata] = Cube_Generate(InputUser.BinFiles{phase},0);
% as a function handle
master = S2FunHandle(@(v) Cube_Sample(v.x(:),v.y(:),v.z(:),screen_int,0));
% as hamonic expansion
%masterHarm = S2FunHarmonicSym.quadrature(master,cs,'bandwidth',50); %run to initialise
masterHarm = S2FunHarmonicSym.quadrature(master,cs,'bandwidth',512);

[ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,num_Phases, RTI_info ] = Phase_Builder_RTI( InputUser.Phase_Input,RTI.Phase_Folder, RTI.Bin_loc );
[screen_int,facedata] = Cube_Generate(RTI_info.bin_file,RTI_info.isHex);

PatternInfo.ScreenWidth=RTI.screensize;
PatternInfo.ScreenHeight=RTI.screensize;
%locs=1;
[ EBSP_pat ] = EBSP_Gnom( PatternInfo,PC);
[ pat ] = EBSP_gen( EBSP_pat,Sphere.rottoplot,screen_int,RTI_info.isHex );

pat=(pat-mean(pat(:)))/std(pat(:));

    %save the gamma and gamma prime pats
%     if phase==1
%         patg=pat;
%     elseif phase==2
%         patgp=pat;
%     end
    

%% Get the correct orientation for MTEX
new_eulers=conv_G_to_EA(Sphere.rottoplot);
ori = orientation.byEuler(-new_eulers(1),-new_eulers(2),-new_eulers(3),'ZYZ',cs);
%ori = orientation.byEuler(-phi1,-PHI,-phi2,'ZYZ',cs);
ori0=orientation.byEuler(0,0,0,cs);

%% Generate MASTER pattern
if phase < 3
[screen_int,facedata] = Cube_Generate(InputUser.BinFiles{phase},0);
% as a function handle
master = S2FunHandle(@(v) Cube_Sample(v.x(:),v.y(:),v.z(:),screen_int,0));
% as hamonic expansion
masterHarm = S2FunHarmonicSym.quadrature(master,cs,'bandwidth',512);

plot(master)
cd(resultsfolder)
colormap('gray')
print(['Master_',num2str(phase)],'-dpng','-r300')
cd(homefol)

figure
plot(masterHarm.radon,'pcolor','resolution',0.25*degree,'upper','complete')
annotate(ori0*h().symmetrise,'markerSize',15,'antipodal','MarkerEdgeColor','c','MarkerFaceColor','none','Marker','o')
colormap('gray')
cd(resultsfolder)
print(['MasterRadon_',num2str(phase)],'-dpng','-r300')
cd(homefol)
end

%% Look at  pattern
det = detector(RTI.screensize,RTI.screensize,MapData.DD(1),[1-PC_Y,1-PC_X]);
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
print(['Radon_Peaks_',num2str(phase)],'-dpng','-r300')
cd(homefol)

%% Re-project pattern onto sphere and highlight strongest bands
plot(pHarm,'pcolor','resolution',0.25*degree,'upper')
colormap('gray')
circle(peaks,'linecolor',[0.35,0.35,0.35])

cd(resultsfolder)
print(['Reprojected_StrongestBands_',num2str(phase)],'-dpng','-r300')
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
print(['Reprojected_OnRefSphere_',num2str(phase)],'-dpng','-r300')
cd(homefol)

%% Do the same for the empty pattern
pHarm_rot_z=rotate(pHarm_z,inv(ori));
figure
plot(pHarm_rot_z,'pcolor','resolution',0.25*degree,'upper')
colormap('gray')

cd(resultsfolder)
print(['EmptyPattern_Rotated',num2str(phase)],'-dpng','-r300')
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
title(num2str(phase))
hold off

cd(resultsfolder)
print(['ReprojectedRadon_SpecificBands_',num2str(phase)],'-dpng','-r300')
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
print(['Reprojected_SpecificBands_',num2str(phase)],'-dpng','-r300')
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
title(num2str(phase))
hold off

cd(resultsfolder)
print(['ReprojectedRadon_rot_SpecificBands_',num2str(phase)],'-dpng','-r300')
cd(homefol)

%%

h = Miller({1,0,0},{1,1,0},{1,3,1},{1,1,1},cs);%{2,0,0},{1,1,2},{1,2,3},cs);
hnames = {'{1,0,0}','{1,1,0}','{1,3,1}','{1,1,1}'};%,'{2,0,0}','{1,1,2}','{1,2,3}'};

figure
plot(pHarm_rot,'pcolor','resolution',0.25*degree,'upper')
colormap('gray')
hold on
for i = 1:length(h)
    circle(h(i).symmetrise,'linecolor',colors(i,:),'LineWidth',1.5)
end
hold off
cd(resultsfolder)
print(['Reprojected_rot_SpecificBands_',num2str(phase)],'-dpng','-r300')
cd(homefol)

%% Integrate along an unrotated empty sphere
for i = 1:length(h)
  [~,psi_pat_Baseline{i}] = pHarm_z.symmetrise(h(i));
  psi_pat_Baseline{i}.A(1) = 0;
end

%% Integrate along a rotated empty sphere
for i = 1:length(h)
  [~,psi_pat_Baseline_rot{i}] = pHarm_rot_z.symmetrise(h(i));
  psi_pat_Baseline_rot{i}.A(1) = 0;
end

%% Integrate along rotated normals on unrotated sphere

%Add up peak intensities for all bands of a given type
for i = 1:length(h)
  g=h(i).symmetrise;
  for j = 1:length(g(i).symmetrise)  
  [~,psi_pat{j}] = pHarm.symmetrise(ori*g(j));
  psi_pat{j}.A(1) = 0;
  end
  
  psi_pats{phase,i}=psi_pat;
  clear psi_pat;
end

%% Integrate along unrotated normals on rotated sphere
for i = 1:length(h)
  g=h(i).symmetrise;
  for j=1:length(g(i).symmetrise)
  [~,psi_pat_Rot{j}] = pHarm_rot.symmetrise(g(j));
  psi_pat_Rot{j}.A(1) = 0;
  end
  
  psi_pats_Rot{phase,i}=psi_pat_Rot;
  clear psi_pat_Rot;
end

%% Integrate along unrotated normals on unrotated master pattern
%only need to do for one because they're all the same...
for i = 1:length(h)
  [~,psi_full{phase,i}] = masterHarm.symmetrise(h(i));
  psi_full{phase,i}.A(1) = 0;
end
%domain=domain.*180/pi;
end