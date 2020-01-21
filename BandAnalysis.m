InputUser.Phases={'Ni_HR','Ni3Al_HR'};
InputUser.BinFiles={'C:\Users\tpm416\Documents\GitHub\RTM_indexing\masterpatterns\Ni_1024.bin','C:\Users\tpm416\Documents\GitHub\RTM_indexing\masterpatterns\Ni3Al_1024.bin'}; 
InputUser.cifnames={'C:\Users\tpm416\Documents\GitHub\AstroEBSD\phases\Ni.cif','C:\Users\tpm416\Documents\GitHub\AstroEBSD\phases\Ni3Al.cif'};
InputUser.ResultsDir='C:\Users\tpm416\Documents\GitHub\GammaPrime_FA\Results';
InputUser.BA_Dir='C:\Users\tpm416\Documents\GitHub\GammaPrime_FA';
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

%add Band analysis directories to path
addpath(addpath([InputUser.BA_Dir,'\DataLoading']))
addpath(addpath([InputUser.BA_Dir,'\EulerAngle_checking']))
addpath(addpath([InputUser.BA_Dir,'\PlotProfiles']))
addpath(addpath([InputUser.BA_Dir,'\SphericalAnalysis']))

%% Load, index and refine an experimental pattern
InputUser.HDF5_file='yprime3.h5';
InputUser.HDF5_folder='E:\Tom\GammaPrime_Data\V208C';
[ MapData,MicroscopeData,~,EBSPData ]=bReadHDF5( InputUser );
[Data_InputMap] = EBSD_Map(MapData,MicroscopeData);

% get reference pattern
[ RefPat ] = bReadEBSP(EBSPData,1);
[ RefPatCor ] = EBSP_BGCor(RefPat,Settings_Cor);

%%
PC_start=[0.4760,0.2540,0.6136]; %initial value for PC
Eulers=[17.4*degree,60*degree,32.6*degree]; %initial value for Eulers

Refine.ss=0.08; %initial probe volume
Refine.p=2; %order of polynomial to fit to tetrahedron
Refine.n_its=150;
Refine.reindex=1;
Refine.print=1;

t1=clock;
pTime('Starting PC refinement',t1);
[Refine]=PC_refine(Eulers,Refine,PC_start,RefPatCor,MicroscopeData,RTI,InputUser);
pTime(['PC refined: ',num2str(100* Refine.Increase),'% increase'],t1);

clear Eulers
%% Define a band profile
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

%Define all rotation matrices needed in the code
GMat_test=conv_EA_to_G([phi1,PHI,phi2]);
Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
Detector_tilt = Rx(tilt);
Sphere.rottoplot=GMat_test*Detector_tilt;
clear Detector_tilt GMat_test Rx

%% Set saving directory
cd(InputUser.ResultsDir)
homefol=pwd;
InputUser.resultsfolder=[InputUser.ResultsDir,'/Results_',InputUser.HDF5_file(1:end-3)];
mkdir(InputUser.resultsfolder)
clear PHI phi1 phi2 PC_start tilt

%% Run the spherical stuff for pure y / y' projections
crystals{1}=loadCIF(InputUser.Phases{1});
crystals{2}=loadCIF(InputUser.Phases{2});

cd(InputUser.resultsfolder)
[Sphere.psi_full,Sphere.psi_pats,Sphere.psi_pats_Rot]=Spherical_Phases(InputUser,RTI,MapData,Refine.PC_out,crystals,Sphere);

%% Load the ML data
cd(InputUser.resultsfolder);
ML_extraction;
AvPat_extraction;

%do some manual flips (currently done in VS Code)
% ML.coeffs_PCA(:,:,1,:)=-ML.coeffs_PCA(:,:,1,:);
% ML.coeffs_PCA(:,:,2,:)=-ML.coeffs_PCA(:,:,2,:);
% %ML.coeffs_PCA(:,:,3,:)=ML.coeffs_PCA(:,:,3,:);
% ML.coeffs_PCA(:,:,4,1)=-ML.coeffs_PCA(:,:,4,1);
% ML.coeffs_PCA(:,:,5,:)=-ML.coeffs_PCA(:,:,5,1);

%% Run the band analysis for different sets of patterns (experimental, NMF, PCA)
cs=crystals{2};

ExpPats=zeros(300,300,1);
ExpPats(:,:,1)=RefPatCor;

%Run on experimental pattern
cd(InputUser.resultsfolder)
mkdir('Experimental')
cd('Experimental')
[Sphere.psi_pats_experiment,~]=Spherical_Experiments(ExpPats,InputUser,RTI,MapData,Refine.PC_out,cs,Sphere);

%%
for i=1:5 %spatial localisation (NMF)
cd(InputUser.resultsfolder)

%Normalise stdev to 1, mean to 0
ExpPats=squeeze(ML.coeffs_NMF(:,:,:,i));
NormPats=ExpPats-ones(size(ExpPats)).*mean(ExpPats,[1,2]);
NormPats=NormPats./std(ExpPats,1,[1,2]);

mkdir(['NMF_spatial_',num2str(i)])
cd(['NMF_spatial_',num2str(i)])
[Sphere.(['psi_pats_experiment_NMF_',num2str(i)]),~]=Spherical_Experiments(NormPats,InputUser,RTI,MapData,Refine.PC_out,cs,Sphere);
end

for i=1:5 %spatial localisation (PCA)
cd(InputUser.resultsfolder)
mkdir(['PCA_spatial_',num2str(i)])
cd(['PCA_spatial_',num2str(i)])

%Normalise stdev to 1, mean to 0
ExpPats=squeeze(ML.coeffs_PCA(:,:,:,i));
NormPats=ExpPats;
%NormPats=ExpPats-ones(size(ExpPats)).*mean(ExpPats,[1,2]); mean should already be zero
NormPats=NormPats./std(ExpPats,1,[1,2]);

%Now multiply by the singular values (lengths of coefficient vectors)
%NormPats=NormPats.*reshape((ML.coefflengths_PCA(:,i)),1,1,5);%all variances, i'th spatial localisation

[Sphere.(['psi_pats_experiment_PCA_',num2str(i)]),~]=Spherical_Experiments(NormPats,InputUser,RTI,MapData,Refine.PC_out,cs,Sphere);
end

%% Average patterns
cd(InputUser.resultsfolder)
cd('Experimental')
mkdir(['PCA_Av_',num2str(1)])
cd(['PCA_Av_',num2str(1)])

NormPats=Av.PCA_Pats;
%NormPats=ExpPats-ones(size(ExpPats)).*mean(ExpPats,[1,2]); mean should already be zero
NormPats=NormPats./std(NormPats,1,[1,2]);
%Now multiply by the singular values (lengths of coefficient vectors)
%NormPats=NormPats.*reshape((ML.coefflengths_PCA(:,i)),1,1,5);%all variances, i'th spatial localisation
[Sphere.(['psi_pats_Average_PCA_',num2str(1)]),~]=Spherical_Experiments(NormPats,InputUser,RTI,MapData,Refine.PC_out,cs,Sphere);

cd(InputUser.resultsfolder)
cd('Experimental')
mkdir(['NMF_av_',num2str(1)])
cd(['NMF_av_',num2str(1)])

NormPats=Av.NMF_Pats;
NormPats=NormPats-ones(size(NormPats)).*mean(NormPats,[1,2]); %mean should already be zero
NormPats=NormPats./std(NormPats,1,[1,2]);
%Now multiply by the singular values (lengths of coefficient vectors)
%NormPats=NormPats.*reshape((ML.coefflengths_PCA(:,i)),1,1,5);%all variances, i'th spatial localisation
[Sphere.(['psi_pats_Average_NMF_',num2str(1)]),~]=Spherical_Experiments(NormPats,InputUser,RTI,MapData,Refine.PC_out,cs,Sphere);

%% Differences
NormPats=Av.NMF_Pats(:,:,1)-Av.NMF_Pats(:,:,2);
NormPats=NormPats-ones(size(NormPats)).*mean(NormPats,[1,2]); %mean should already be zero
NormPats=NormPats./std(NormPats,1,[1,2]);
%Now multiply by the singular values (lengths of coefficient vectors)
%NormPats=NormPats.*reshape((ML.coefflengths_PCA(:,i)),1,1,5);%all variances, i'th spatial localisation
[Sphere.(['psi_pats_Average_NMF_Diff_',num2str(1)]),~]=Spherical_Experiments(NormPats,InputUser,RTI,MapData,Refine.PC_out,cs,Sphere);

NormPats=Av.PCA_Pats(:,:,1)-Av.PCA_Pats(:,:,2);
%NormPats=NormPats-ones(size(NormPats)).*mean(NormPats,[1,2]); %mean should already be zero
NormPats=NormPats./std(NormPats,1,[1,2]);
%Now multiply by the singular values (lengths of coefficient vectors)
%NormPats=NormPats.*reshape((ML.coefflengths_PCA(:,i)),1,1,5);%all variances, i'th spatial localisation
[Sphere.(['psi_pats_Average_PCA_Diff_',num2str(1)]),~]=Spherical_Experiments(NormPats,InputUser,RTI,MapData,Refine.PC_out,cs,Sphere);

%% Plot band profiles of simulated patterns
h = Miller({1,0,0},{1,1,0},{1,3,1},{1,1,1},cs);%{2,0,0},{1,1,2},{1,2,3},cs);
hnames = {'{1,0,0}','{1,1,0}','{1,3,1}','{1,1,1}'};%,'{2,0,0}','{1,1,2}','{1,2,3}'};

cd(InputUser.resultsfolder)
PlotPatterns(h,hnames,Sphere.psi_pats,'Simulation_',{'Ni','Ni_3Al'},[-0.3,0.3],1)
PlotPatterns(h,hnames,Sphere.psi_full,'GlobalSimulation_',{'Ni','Ni_3Al'},[-0.3,0.3],1)
close all

%% Plot band profiles of some experimental sets - NMF (spatial kernel = 1)
%kernels=Sphere.psi_pats_experiment([list],:);
list_names={'NMF1','NMF2','NMF3'};
PlotPatterns(h,hnames,Sphere.psi_pats_experiment_NMF_1,'NMF_local1_',list_names,[-0.3,0.3],1)

%% Plot band profiles of some experimental sets - NMF (all NMF factor 2, 110, for kernels 1:5)
kernels=[Sphere.psi_pats_experiment_NMF_1(2,:);Sphere.psi_pats_experiment_NMF_2(2,:);Sphere.psi_pats_experiment_NMF_3(2,:);...
    Sphere.psi_pats_experiment_NMF_4(2,:);Sphere.psi_pats_experiment_NMF_5(2,:)];
list_names={'Kernel 1','Kernel 2','Kernel 3','Kernel 4', 'Kernel 5'};
PlotPatterns(h,hnames,kernels,'NMF_variedKernel_',list_names,[-0.3,0.3],1)

%% Plot band profiles of some experimental sets - PCA (spatial kernel = 1)
list_names={'PCA1','PCA2','PCA3','PCA4','PCA5'};
PlotPatterns(h,hnames,Sphere.psi_pats_experiment_PCA_1,'PCA_local1_',list_names,[-0.3,0.3],1)

%% Plot band profiles of some experimental sets - PCA (all PC2s, 110, for kernels 1:5)
kernels=[Sphere.psi_pats_experiment_PCA_1(2,:);Sphere.psi_pats_experiment_PCA_2(2,:);Sphere.psi_pats_experiment_PCA_3(2,:);...
    Sphere.psi_pats_experiment_PCA_4(2,:);Sphere.psi_pats_experiment_PCA_5(2,:)];
list_names={'Kernel 1','Kernel 2','Kernel 3','Kernel 4', 'Kernel 5'};
PlotPatterns(h,hnames,kernels,'PCA_variedKernel_',list_names,[-0.3,0.3],1)

%% Compare NMF to simulation
for n=1:5
for i=1:3
    cd(InputUser.resultsfolder)
    mkdir(['NMF_SimulationCom_spatial',num2str(n)])
    cd(['NMF_SimulationCom_spatial',num2str(n)])
    SpherComp=Sphere.(['psi_pats_experiment_NMF_',num2str(n)]);
    kernels=[Sphere.psi_pats(1,:);Sphere.psi_pats(2,:);SpherComp(i,:)]; %gamma, gamma prime, %NMF factor
    list_names={'Ni','Ni_3Al', ['NMF factor ',num2str(i)]};
    PlotPatterns(h,hnames,kernels,['NMF_SimComp_Factor',num2str(i)],list_names,[-0.3,0.3],1)
end
end
cd(InputUser.resultsfolder)
close all

%% Compare PCA to simulation
for n=1:5
for i=1:5
    cd(InputUser.resultsfolder)
    mkdir(['PCA_SimulationCom_spatial',num2str(n)])
    cd(['PCA_SimulationCom_spatial',num2str(n)])
    SpherComp=Sphere.(['psi_pats_experiment_PCA_',num2str(n)]);
    kernels=[Sphere.psi_pats(1,:);Sphere.psi_pats(2,:);SpherComp(i,:)]; %gamma, gamma prime, %NMF factor
    list_names={'Ni','Ni_3Al', ['PC ',num2str(i)]};
    PlotPatterns(h,hnames,kernels,['PCA_SimComp_Factor',num2str(i)],list_names,[-0.3,0.3],1)
end
end
cd(InputUser.resultsfolder)
close all

%% Compare averages to simulation
n=1;
%PCA
cd(InputUser.resultsfolder)
mkdir(['AvPats_SimulationCom_spatial',num2str(n)])
cd(['AvPats_SimulationCom_spatial',num2str(n)])
SpherComp=Sphere.(['psi_pats_Average_PCA_',num2str(n)]);
kernels=[Sphere.psi_pats(1,:);Sphere.psi_pats(2,:);SpherComp(1,:);SpherComp(2,:)]; %gamma, gamma prime, %NMF factor
list_names={'Ni','Ni_3Al', 'Ni_3Al Exp ','Ni Exp'};
PlotPatterns(h,hnames,kernels,['AvPCA_SimComp'],list_names,[-0.3,0.3],1)

cd(InputUser.resultsfolder)
mkdir(['AvPats_SimulationCom_spatial',num2str(n)])
cd(['AvPats_SimulationCom_spatial',num2str(n)])
SpherComp=Sphere.(['psi_pats_Average_NMF_',num2str(n)]);
kernels=[Sphere.psi_pats(1,:);Sphere.psi_pats(2,:);SpherComp(1,:);SpherComp(2,:)]; %gamma, gamma prime, %NMF factor
list_names={'Ni','Ni_3Al','Ni_3Al Exp ','Ni Exp'};
PlotPatterns(h,hnames,kernels,['AvNMF_SimComp'],list_names,[-0.3,0.3],1)

cd(InputUser.resultsfolder)
mkdir(['AvPats_SimulationCom_spatial',num2str(n)])
cd(['AvPats_SimulationCom_spatial',num2str(n)])
SpherComp=Sphere.(['psi_pats_Average_NMF_Diff_',num2str(n)]);
kernels=[Sphere.psi_pats(1,:);Sphere.psi_pats(2,:);SpherComp(1,:)]; %gamma, gamma prime, %NMF factor
list_names={'Ni','Ni_3Al','Diff'};
PlotPatterns(h,hnames,kernels,['AvNMF_SimComp_Diff'],list_names,[-0.3,0.3],1)

cd(InputUser.resultsfolder)
mkdir(['AvPats_SimulationCom_spatial',num2str(n)])
cd(['AvPats_SimulationCom_spatial',num2str(n)])
SpherComp=Sphere.(['psi_pats_Average_PCA_Diff_',num2str(n)]);
kernels=[Sphere.psi_pats(1,:);Sphere.psi_pats(2,:);SpherComp(1,:)]; %gamma, gamma prime, %NMF factor
list_names={'Ni','Ni_3Al','Diff'};
PlotPatterns(h,hnames,kernels,['AvPCA_SimComp_Diff'],list_names,[-0.3,0.3],1)

cd(InputUser.resultsfolder)
close all

%%
cd(InputUser.resultsfolder)
save('Results.mat')
SpatialStatistics
cd(homefol)
close all
