InputUser.Phases={'Ni','Ni3Al'};
BinFiles={'C:\Users\tpm416\Documents\GitHub\RTM_indexing\masterpatterns\Ni_1024.bin','C:\Users\tpm416\Documents\GitHub\RTM_indexing\masterpatterns\Ni3Al_1024.bin'}; 
cifnames={'C:\Users\tpm416\Documents\GitHub\AstroEBSD\phases\Ni.cif','C:\Users\tpm416\Documents\GitHub\AstroEBSD\phases\Ni3Al.cif'};
isHexes={0,0};

Settings_Cor.gfilt=1; %use a high pass filter (1)
Settings_Cor.gfilt_s=4; %low pass filter sigma (4)
Settings_Cor.radius=0; %use a radius mask (0)
Settings_Cor.radius_frac=0.85; %fraction of the pattern width to use as the mask (0.85)
Settings_Cor.hotpixel=1; %hot pixel correction (1)
Settings_Cor.hot_thresh=1000; %hot pixel threshold (1000)
Settings_Cor.resize=1; %resize correction (1)
Settings_Cor.size=300; %image height
Settings_Cor.SplitBG=1; %split chip fix (1)
Settings_Cor.RealBG=0; %use a real BG (0)
%Settings_Cor.EBSP_bgnum=30; (30)
Settings_Cor.Square=1; %square crop (1)
Settings_Cor.SquareCrop=1; %(1)
Settings_Cor.LineError=1; %(1)
Settings_Cor.MeanCentre=1; %(1)

RTI.Bin_loc='C:\Users\tpm416\Documents\GitHub\RTM_indexing\masterpatterns';
RTI.Phase_Folder='C:\Users\tpm416\Documents\GitHub\AstroEBSD\phases';
RTI.screensize=Settings_Cor.size;
RTI.Sampling_Freq=8; %of SO3 space
RTI.XCF_type=2; %will be removed
RTI.iterations=4;
RTI.LPTsize=500;

run('C:\Users\tpm416\Documents\GitHub\AstroEBSD\start_AstroEBSD')
run('C:\Users\tpm416\Documents\GitHub\RTM_indexing\start_RTI')
run('C:\Users\tpm416\Documents\GitHub\EBSD-EDS_PCA\start_PCA')
run('C:\Communal_MatlabPlugins\mtex-5.2.beta2\startup_mtex')

addpath('C:\Users\tpm416\Documents\GitHub\chebfun-master')
%% Load and index an experimental pattern
InputUser.HDF5_file='yprime3.h5';
InputUser.HDF5_folder='C:\Users\tpm416\Documents\Data\V208C';
[ MapData,MicroscopeData,~,EBSPData ]=bReadHDF5( InputUser );
[Data_InputMap] = EBSD_Map(MapData,MicroscopeData);

% get reference pattern
[ RefPat ] = bReadEBSP(EBSPData,1);
[ RefPatCor ] = EBSP_BGCor(RefPat,Settings_Cor);

% ... and index it
phase=1;
InputUser.Phase_Input=InputUser.Phases(phase);
cs=loadCIF(InputUser.Phases{phase});
[ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,num_Phases, RTI_info ] = Phase_Builder_RTI( InputUser.Phase_Input,RTI.Phase_Folder, RTI.Bin_loc );
[screen_int,facedata] = Cube_Generate(RTI_info.bin_file,RTI_info.isHex);

MapSize1=max(MapData.YBeam)-min(MapData.YBeam)+1; % number of rows
MapSize2=max(MapData.XBeam)-min(MapData.XBeam)+1; % number of columns
SEMImage=transpose(MicroscopeData.SEMImage(:,:,1));

time1=clock;
screensize=RTI.screensize;
Sampling_Freq=RTI.Sampling_Freq;
XCF_type=RTI.XCF_type;
iterations=RTI.iterations;
LPTsize=RTI.LPTsize;
Bin_loc=RTI.Bin_loc;
Phase_Folder=RTI.Phase_Folder;

% Set up some things the same way RTI does
Data_InputMap_Start=Data_InputMap;
EBSD_DataInfo=EBSPData;
Settings_CorX=Settings_Cor;

% Get refined PC
%Index_Pat
%PC_refined=PC;

%PC_refined=[0.515,0.27,Data_InputMap.DD(3)];
PC_refined=[0.4569,0.2651,0.6358];

InputPats=RefPatCor;
for phasenum=1:2
    InputUser.Phase_Input=InputUser.Phases(phasenum);
    RTI_run

    RTI.peakheight(phasenum,:)=TemData(7,:);
    RTI.Phi1(phasenum,:)=TemData(1,:);
    RTI.PHI(phasenum,:)=TemData(2,:);
    RTI.Phi2(phasenum,:)=TemData(3,:);
end

clear Crystal_Family Crystal_LUT Crystal_UCell cs Data_InputMap_Start Detector_tilt EBSD_av EBSP_av facedata pat phasenum
clear phaseassign_val phaseassign_loc Euler1 Euler2 Euler3 LibNo RTI_peakheight TemData a screensize Bin_loc Sampling_Freq XCF_type iterations LPTsize Phase_Folder
clear library_G num_Phases PC_av PC_X PC_Y PC_Z PCAPats phase PhaseData RTI_info Rx Ry Rz screen_int screen_int_par Settings_LUT SettingsXCF SettingsXCF2
clear template_library template_library_par RTIOutput time1 librarypats matchedpats
clear Euler1_assign Euler2_assign Euler3_assign loc correction
clear exist_components exist_inputpats ComponentMap InputPats PCXs PCYs PMaps DDs Settings_CorX skipRead xi ret_comps Watershed_selcomps
clear testlib_size
clear MapSize1 MapSize 2 ans

Mtex_GenPat_comparison

%% Define a band profile
% modified Gaussian profile
profile = @(x) 1.5*exp(-(acos(x)-90*degree).^2./(3*degree).^2) - ...
  exp(-(acos(x)-87*degree).^2./(2*degree).^2) - ...
  exp(-(acos(x)-93*degree).^2./(2*degree).^2);

% expand into a Legendre series
profileHarm = S2Kernel.quadrature(profile);

%% Set up candidate orientation & PC
phi1=RTI.Phi1(1);
PHI=RTI.PHI(1);
phi2=RTI.Phi2(1);

tilt=MicroscopeData.TotalTilt;
%tilt=0;

PC=PC_refined;

PC_X=PC(1);
PC_Y=PC(2);
PC_Z=PC(3);

%Define all rotation matrices needed in the code
GMat_test=conv_EA_to_G([phi1,PHI,phi2]);
Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
Detector_tilt = Rx(tilt);
rottoplot=GMat_test*Detector_tilt;

%% Set saving directory

homefol=pwd;
resultsfolder=['Eulers_',num2str(phi1*180/pi,3),'_',num2str(PHI*180/pi,3),'_',num2str(phi2*180/pi,3),'_tilt_',num2str(90+tilt*180/pi)];
mkdir(resultsfolder)

%% 

for phase=1:4

%% Define some interesting bands for summing
if phase < 3

InputUser.Phase_Input=InputUser.Phases(phase);
cs=loadCIF(InputUser.Phases{phase});
h = Miller({1,0,0},{1,1,0},{1,3,1},{4,3,1},cs);%{2,0,0},{1,1,2},{1,2,3},cs);
hnames = {'{1,0,0}','{1,1,0}','{1,3,1}','{4,3,1}'};%,'{2,0,0}','{1,1,2}','{1,2,3}'};

%% Generate a pattern

% Master patterns
[screen_int,facedata] = Cube_Generate(BinFiles{phase},0);
% as a function handle
master = S2FunHandle(@(v) Cube_Sample(v.x(:),v.y(:),v.z(:),screen_int,0));
% as hamonic expansion
masterHarm = S2FunHarmonicSym.quadrature(master,cs,'bandwidth',512);

[ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,num_Phases, RTI_info ] = Phase_Builder_RTI( InputUser.Phase_Input,RTI.Phase_Folder, RTI.Bin_loc );
[screen_int,facedata] = Cube_Generate(RTI_info.bin_file,RTI_info.isHex);

PatternInfo.ScreenWidth=RTI.screensize;
PatternInfo.ScreenHeight=RTI.screensize;
%locs=1;
[ EBSP_pat ] = EBSP_Gnom( PatternInfo,PC);
[ pat ] = EBSP_gen( EBSP_pat,rottoplot,screen_int,RTI_info.isHex );

pat=(pat-mean(pat(:)))/std(pat(:));

    %save the gamma and gamma prime pats
    if phase==1
        patg=pat;
    elseif phase==2
        patgp=pat;
    end
    
elseif phase==3
    
pat=RefPatCor;

elseif phase==4    
pat=patgp-patg;
end


%% Get the correct orientation for MTEX
new_eulers=conv_G_to_EA(rottoplot);
ori = orientation.byEuler(-new_eulers(1),-new_eulers(2),-new_eulers(3),'ZYZ',cs);
%ori = orientation.byEuler(-phi1,-PHI,-phi2,'ZYZ',cs);
ori0=orientation.byEuler(0,0,0,cs);

%% Generate MASTER pattern
if phase < 3
[screen_int,facedata] = Cube_Generate(BinFiles{phase},0);
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
RadonPHarm = conv(pHarm,profileHarm)./ correction;
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
RadonPHarm_rot = conv(pHarm_rot,profileHarm)./ correction;
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

figure
plot(pHarm,'pcolor','resolution',0.25*degree,'upper')
colormap('gray')
hold on
for i = 1:length(h)
    l{i}=circle(ori*h(i).symmetrise,'linecolor',colors(i,:),'LineWidth',0.75);
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

figure
plot(pHarm_rot,'pcolor','resolution',0.25*degree,'upper')
colormap('gray')
hold on
for i = 1:length(h)
    circle(h(i).symmetrise,'linecolor',colors(i,:),'LineWidth',0.75)
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
for i = 1:length(h)
  [~,psi_pat{phase,i}] = pHarm.symmetrise(ori*h(i));
  psi_pat{phase,i}.A(1) = 0;
end

%% Integrate along unrotated normals on rotated sphere
for i = 1:length(h)
  [~,psi_pat_Rot{phase,i}] = pHarm_rot.symmetrise(h(i));
  psi_pat_Rot{phase,i}.A(1) = 0;
end

%% Integrate along unrotated normals on unrotated master pattern
for i = 1:length(h)
  [~,psi_full{phase,i}] = masterHarm.symmetrise(h(i));
  psi_full{phase,i}.A(1) = 0;
end
domain=domain.*180/pi;
end
% Plot band profiles
%for single, unrotated pat, rotated planes
for i = 1:length(h)
    figure
    title(hnames{i})
    for phase = 1:3
        if phase==1
        plot(psi_pat{phase,i},'DisplayName',char(h(i)),'linewidth',1,'LineStyle','-')
        hold on
        elseif phase==2
        plot(psi_pat{phase,i},'DisplayName',char(h(i)),'linewidth',1,'LineStyle','-')
        hold on
        elseif phase==3
        plot(psi_pat{phase,i},'DisplayName',char(h(i)),'linewidth',1,'LineStyle','-')
        hold on  
        end
    end
    legend('Ni','Ni3Al','Measured')
    xlim([80,100])
    ylim([-0.3,0.3])
    hold off
    title(hnames{i})
    
    cd(resultsfolder)
    print(['Integral_RotPlanes_',hnames{i}],'-dpng','-r300')
    cd(homefol)
    
end

%for single, rotated pat, unrotated planes
for i = 1:length(h)
    figure
    title(hnames{i})
    for phase = 1:3
        if phase==1
        plot(psi_pat_Rot{phase,i},'DisplayName',char(h(i)),'linewidth',1,'LineStyle','-')
        hold on
        elseif phase==2
        plot(psi_pat_Rot{phase,i},'DisplayName',char(h(i)),'linewidth',1,'LineStyle','-')
        hold on
        elseif phase==3
        plot(psi_pat_Rot{phase,i},'DisplayName',char(h(i)),'linewidth',1,'LineStyle','-')
        hold on  
        end
    end
    legend('Ni','Ni3Al','Measured')
    xlim([80,100])
    ylim([-0.3,0.3])
    hold off
    title(hnames{i})
    
    cd(resultsfolder)
    print(['Integral_RotPattern_',hnames{i}],'-dpng','-r300')
    cd(homefol)
    
end

%for zero 
for i = 1:length(h)
    figure
    title(hnames{i})
    plot(psi_pat_Baseline{i},'DisplayName',char(h(i)),'linewidth',1,'LineStyle','-')
    hold on

    xlim([80,100])
    ylim([-0.3,0.3])
    hold off
    title(hnames{i})
    
    cd(resultsfolder)
    print(['Integral_Zeros_',hnames{i}],'-dpng','-r300')
    cd(homefol) 
end

%for full library
for i = 1:length(h)
    figure
    title(hnames{i})
    for phase = 1:2
        if phase==1
        plot(psi_full{phase,i},'DisplayName',char(h(i)),'linewidth',1,'LineStyle','-')
        hold on
        elseif phase==2
        plot(psi_full{phase,i},'DisplayName',char(h(i)),'linewidth',1,'LineStyle','-')
        hold on  
        end
    end
    legend('Ni','Ni3Al')
    xlim([80,100])
    ylim([-0.3,0.3])
    hold off
    title(hnames{i})
    
    cd(resultsfolder)
    print(['MasterIntegral_',hnames{i}],'-dpng','-r300')
    cd(homefol)
end

% plot difference
for i = 1:length(h)
    figure
    plot(psi_pat{4,i},'DisplayName',char(h(i)),'linewidth',1,'LineStyle','-')
    xlim([80,100])
    ylim([-0.05,0.05])
    title(hnames{i})
    
    cd(resultsfolder)
    print(['Integral_RotPlanes_diff_',hnames{i}],'-dpng','-r300')
    cd(homefol)
    
end

close all