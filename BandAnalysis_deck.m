%% TPM 2020
%set(groot, 'defaultAxesFontName','Helvetica')
%set(groot, 'defaultAxesFontSize',12)

InputUser.Phases={'CoNi_HR','Co3AlW_HR','Ni_HR','Ni3Al_HR'};
InputUser.PhaseLabels={'CoNi','Co_3AlW','Ni','Ni_3Al'};

InputUser.ResultsDir='C:\Users\tpm416\Documents\GitHub\BandAnalysis\Results';
InputUser.BA_Dir='C:\Users\tpm416\Documents\GitHub\BandAnalysis';

%% Settings
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

%RTI settings
RTI.Phase_Folder='C:\Users\tpm416\Documents\GitHub\AstroEBSD\phases';
RTI.screensize=Settings_Cor.size;
RTI.Sampling_Freq=8; %of SO3 space
RTI.iterations=4;
RTI.LPTsize=500;

%add bits to path
run('C:\Users\tpm416\Documents\GitHub\AstroEBSD\start_AstroEBSD')
run('C:\Communal_MatlabPlugins\mtex-5.2.beta2\startup_mtex')

%need the chebyshev functions
addpath('C:\Users\tpm416\Documents\GitHub\chebfun-master')

%add Band analysis directories to path
addpath(addpath([InputUser.BA_Dir,'\PlotProfiles']))
addpath(addpath([InputUser.BA_Dir,'\SphericalAnalysis']))

%% Load a single experimental pattern for PC optimisation
InputUser.HDF5_file='yprime3.h5';
InputUser.HDF5_folder='E:\Tom\GammaPrime_Data\V208C';
[ MapData,MicroscopeData,~,EBSPData ]=bReadHDF5( InputUser );
[Data_InputMap] = EBSD_Map(MapData,MicroscopeData);

% get reference pattern
[ RefPat ] = bReadEBSP(EBSPData,1);
[ RefPatCor ] = EBSP_BGCor(RefPat,Settings_Cor);

%% Set up crystal symmetry and bands to index

for i=1:length(InputUser.Phases)
    [ ~,~,~,~,~, RTI_info ] = Phase_Builder_RTM(InputUser.Phases(i),RTI.Phase_Folder);
    InputUser.BinFiles{1,i}=RTI_info.bin_file;
    InputUser.cifnames{1,i}=RTI_info.cif_file;
    crystals{i}=loadCIF(RTI_info.cif_file);
end


%% Single PC refinement (nb. TPM original version, not TBB version - both are in Astro_v2) - this needs updating 
PC_start=[0.4760,0.2540,0.6136]; %initial value for PC
Eulers=[17.4*degree,60*degree,32.6*degree]; %initial value for Eulers

Refine.ss=0.08; %initial probe volume
Refine.p=2; %order of polynomial to fit to tetrahedron
Refine.n_its=150;
Refine.reindex=1;
Refine.print=1;

t1=clock;
pTime('Starting PC refinement',t1);
[Refine]=PC_refine(Eulers,Refine,PC_start,RefPatCor,MicroscopeData,RTI,InputUser.Phases{1});
pTime(['PC refined: ',num2str(100* Refine.Increase),'% increase'],t1);

clear Eulers

%% Define a band profile

%  Grab a symmetry and set of planes
cs=crystals{2};
h = Miller({1,0,0},{1,1,0},{1,3,1},{1,1,1},cs);

% modified Gaussian profile
Sphere.profile = @(x) 1.5*exp(-(acos(x)-90*degree).^2./(3*degree).^2) - ...
  exp(-(acos(x)-87*degree).^2./(2*degree).^2) - ...
  exp(-(acos(x)-93*degree).^2./(2*degree).^2);

% expand into a Legendre series
Sphere.profileHarm = S2Kernel.quadrature(Sphere.profile);

%% Set up candidate orientation & PC
phi1=Refine.Eulers_out(1);PHI=Refine.Eulers_out(2);phi2=Refine.Eulers_out(3);
tilt=MicroscopeData.TotalTilt;
%PC=Refine.PC_out;

GMat_test=conv_EA_to_G([phi1,PHI,phi2]);
Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
Detector_tilt = Rx(tilt);
Sphere.rottoplot=GMat_test*Detector_tilt;

%second test orientation
GMat_test=conv_EA_to_G([2,0.8,1]);
rot2=GMat_test*Detector_tilt;
clear Detector_tilt GMat_test Rx

%% Set saving directory
cd(InputUser.ResultsDir)
homefol=pwd;
InputUser.resultsfolder=[InputUser.ResultsDir,'/Results_',InputUser.HDF5_file(1:end-3)];
mkdir(InputUser.resultsfolder)
clear PHI phi1 phi2 PC_start tilt

%% -------CALCULATION------- %%
%% Firstly, run the spherical analysis for simulated patterns

% - NEED TO RUN, CANCEL, RUN AGAIN FOR THIS TO WORK!! - %

%Test all four crystal structures against all plane families in h

%Spherical_Phases function takes performs a global and singe template simulation + reprojection for each member of InputUser.Phases.
    %Band profile integrations are done for each plane family in h 
    %Saves in directory with suffix given by final input
    %Three outputs: Band profile integrations for global sim ; integrations for sphere rotated onto pattern ; integrations for pattern rotated onto sphere

%Spherical_Phases_Diffs is the same but for difference patterns between phases - output is integrations of pattern differences 1-2 ; 3-4 etc


cd(InputUser.resultsfolder)
[Sphere.psi_full,Sphere.psi_pats,Sphere.psi_pats_Rot]=Spherical_Phases(h,InputUser,RTI,MapData,Refine.PC_out,crystals,Sphere,Sphere.rottoplot,'Experimental');

%Test a second orientation
[~,Sphere.psi_pats_secondori,~]=Spherical_Phases(h,InputUser,RTI,MapData,Refine.PC_out,crystals(1:2),Sphere,rot2,'SecondOri');

%Test differences
[~,Sphere.psi_pats_diff,~]=Spherical_Phases_Diffs(h,InputUser,RTI,MapData,Refine.PC_out,crystals,Sphere,Sphere.rottoplot,'Differences');

%Outputs are eg. 2*4 cells - first dim is phases, second dim is plane
%families. Within each cell is eg 1*6 cell, corresponding to results from 6
%symmetric equivalent members of that {hkl}.

%% Run the band analysis for the experimental dataset

ExpPats=zeros(300,300,1);
ExpPats(:,:,1)=RefPatCor;

%Spherical_Experiments performs analysis for all members of an input set of
%patterns

%Run on experimental pattern
cd(InputUser.resultsfolder)
mkdir('Experimental')
cd('Experimental')
[Sphere.psi_pats_experiment,~]=Spherical_Experiments(ExpPats,h,InputUser,RTI,MapData,Refine.PC_out,cs,Sphere,Sphere.rottoplot);


%% -------PLOTTING-------
%% Colourmap
colours=cbrewer('qual','Paired',10);
colours(3,:)=colours(7,:);
colours(4,:)=colours(8,:);

%% Plot band profiles of simulated patterns
cd(InputUser.resultsfolder)

%Linestyles
Settings.LS={'-','-','-','-'};
%Linecolours
Settings.LC={colours(1,:),colours(2,:),colours(3,:),colours(4,:)};

%Plot this set of patterns for all present plane families
%Additional inputs are labels, xlims, and saving
%Plots each member of h separately. Can combine multiple plots into second
%arg (see below)
PlotPatterns(h,Sphere.psi_pats,'Simulation_',InputUser.PhaseLabels,[-0.3,0.3],1,Settings)
PlotPatterns(h,Sphere.psi_full,'GlobalSimulation_',InputUser.PhaseLabels,[-0.3,0.3],1,Settings)
close all


%% Compare orientations

cd(InputUser.resultsfolder)
mkdir(['OrientationComparison_',num2str(1)])
cd(['OrientationComparison_',num2str(1)])

%Can provide a set of band profiles to compare when feeding into plot
%patterns
%here- plot on the same graphs in order: Co gamma ori1; Co gamma prime
%ori1 ; Co gamma ori2; Co gamma prime ori2.
kernels=[Sphere.psi_pats(1,:);Sphere.psi_pats(2,:);Sphere.psi_pats_secondori(1,:);Sphere.psi_pats_secondori(2,:)]; %gamma, gamma prime (Co) for two ori's

%Labels
list_names={InputUser.PhaseLabels{1:2},[InputUser.PhaseLabels{1},' 2'],[InputUser.PhaseLabels{2},' 2']};
%Line spec
Settings.LS={'-','-','-','-'};
Settings.LC={colours(1,:),colours(2,:),colours(3,:),colours(4,:),[0,0,0]};

%Plot a specific comparison
PlotPatterns(h,kernels,['AvPCA_SimComp_Diff'],list_names,[-0.3,0.3],1,Settings)

%%
cd(InputUser.resultsfolder)
save('Results.mat')

cd('C:\Users\tpm416\Documents\GitHub\BandAnalysis')
SpatialStatistics
cd(homefol)
close all

