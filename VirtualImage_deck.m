%% TPM 23/01/20 - virtual imaging

InputUser.Phases={'CoNi_HR','Co3AlW_HR','Ni_HR','Ni3Al_HR'};
InputUser.PhaseLabels={'CoNi','Co_3AlW','Ni','Ni_3Al'};

InputUser.ResultsDir='C:\Users\tpm416\Documents\GitHub\BandAnalysis\Results';
InputUser.BA_Dir='C:\Users\tpm416\Documents\GitHub\BandAnalysis';

% Set up and go to a saving location
cd('C:\Users\tpm416\Documents\GitHub\BandAnalysis\Results\Results_yprime3')
mkdir('VirtualImages')
cd('VirtualImages')

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
[Refine]=PC_refine(Eulers,Refine,PC_start,RefPatCor,MicroscopeData,RTI,'Ni');
pTime(['PC refined: ',num2str(100* Refine.Increase),'% increase'],t1);

clear Eulers

%% Get orientation matrix
GMat_test=conv_EA_to_G([Refine.Eulers_out(1),Refine.Eulers_out(2),Refine.Eulers_out(3)]);
Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
OriMatrix=GMat_test*Rx(MicroscopeData.TotalTilt);

for i=1:length(InputUser.Phases)
    [ ~,~,~,~,~, RTI_info ] = Phase_Builder_RTM(InputUser.Phases(i),RTI.Phase_Folder);
    InputUser.BinFiles{1,i}=RTI_info.bin_file;
    InputUser.cifnames{1,i}=RTI_info.cif_file;
    crystals{i}=loadCIF(RTI_info.cif_file);
end

cs=crystals{2};

%% Need to run this, cancel, then run again.
[screen_int,~] = Cube_Generate(InputUser.BinFiles{1},0);
% as a function handle
master = S2FunHandle(@(v) Cube_Sample(v.x(:),v.y(:),v.z(:),screen_int,0));
% as hamonic expansion
masterHarm = S2FunHarmonicSym.quadrature(master,cs,'bandwidth',50); %run to initialise

%% - VIRTUAL IMAGING - %%
%% Virtual imaging for a group of plane families and functions:

delta=5; %acquistion window - this is what will be collected for subsequent analysis (eg. Bragg windowing)

SamplingNum=10000; %number of discrete points sampled within delta from band centre.

%Define bands to operate on
h1=Miller({1,1,1},cs);
h2=Miller({1,1,0},cs);
h3=Miller({1,0,0},cs);
h4=Miller({1,3,1},cs);
h_set=[h1,h2,h3,h4];

% get electron wavelength for Bragg angle calc later
lambda=6.626e-24/(2*9.109e-31*1.602e-19*20000)^0.5; % in angstroms

%% - COLLECTION - %%
%% Get the profiles for a given set of bands for the full map
%NOTE - this does NOT reindex for each scan point. It assumes the bands are
%stationary and only varying with intensity (as was approximately the case for yprime
%dataset)

% If printing = 0 this is quite fast for a single pattern
% Core function is AngularAperture(Pattern, delta, SamplingNum, h, symmetry, PC, Orientation, printing, lambda)
% This collects the integrated profile on the band specified by h for
% Pattern

for testno=1:length(h_set)
disp(['Starting analysis ',num2str(testno),' of ',num2str(length(h_set))]);

%Grab the right band
h=h_set(testno);

%% Plot for one pattern and show the highlighted region
[ RefPat ] = bReadEBSP(EBSPData,1);
[ RefPat_cor] = EBSP_BGCor( RefPat,Settings_Cor );
[I,sample]=AngularAperture(RefPat_cor,delta,SamplingNum,h,cs,Refine.PC_out,OriMatrix,1,lambda);

%% Load the data and do the angular aperturing, then apply the defined function (PARALLEL)
rows=Data_InputMap.ypts;cols=Data_InputMap.xpts;
Profile=zeros(rows,cols,SamplingNum);
t1=clock;

pTime('Starting profile calculation...',t1);
for i=1:rows
    parfor j=1:cols
        pattern_number=Data_InputMap.PMap(i,j);
        [ RefPat ] = bReadEBSP(EBSPData,pattern_number);
        [ RefPat_cor] = EBSP_BGCor( RefPat,Settings_Cor );
        
        %Get the profile
        [Profile(i,j,:),sample]=AngularAperture(RefPat_cor,delta,SamplingNum,h,cs,Refine.PC_out,OriMatrix,0,lambda); %get the angular profile - only needs lambda for plotting
        
    end
    pTime(['Completed row ',num2str(i),' of ',num2str(rows)],t1);
end

band=char(h);
Results.(['Profile_',band(2:end-1)])=Profile;
clear VImage band MinimaMaps MaximaMaps
end
save('Results.mat','Results','delta','Refine','-v7.3');

%%
cd('C:\Users\tpm416\Documents\GitHub\BandAnalysis\Results\Results_yprime3\VirtualImages')
load('Results.mat')

%% - POST ANALYSIS - %% - this is where operations are performed, eg. summing within a Bragg window
%% Apply a specific functions (RAW)

%Define the operating function and the dataset to be used
VIm_fn=@(I) sum(I);

Band={'111','110','100','131'}; %label the bands
bragg_delta=asin(lambda./(2.*h_set.dspacing))*180/pi; % in degrees
widths=bragg_delta./(delta/SamplingNum); %number of points in the Bragg width

colors={cbrewer('seq','Purples',1000),cbrewer('seq','Blues',1000),cbrewer('seq','Reds',1000),cbrewer('seq','Greens',1000)};

% Apply the defined function
VImage=zeros(rows,cols,4);

for bandno = 1:4
I=Results.(['Profile_',Band{bandno}]);
rows=size(I,1);
cols=size(I,2);

for i =1:rows
    for j=1:cols
        %VImage(i,j)=VIm_fn(I(i,j,round(5000-widths(bandno)):round(5000+widths(bandno))));
        %apply function and normalise wrt bin width
        VImage(i,j,bandno)=VIm_fn(I(i,j,sample>90-bragg_delta(1)&sample<90+bragg_delta(1)))*(sample(2)-sample(1));
    end
end

cmap=colors{bandno};

imagesc(VImage(:,:,bandno))
colormap(cmap)
colorbar
axis image off;

print(gcf,['BraggSum_',Band{bandno}],'-dpng','-r300')
end

%% Map width of separation (can be useful - may be bugged)

% I=Results.Profile_131;
% n_max=2;
% 
% VImage=zeros(rows,cols);
% for i=1:rows
%     for j=1:cols
%         
%         [p,l]=findpeaks(squeeze(-I(i,j,:)),'SortStr','descend','NPeaks',n_max);
%         diff=abs(l(2)-l(1));
%         VImage(i,j)=diff; %save the locations
%     
%     end
% end
% 
% cmap=cbrewer('seq','YlOrRd',1000);
% imagesc(VImage)
% colormap(cmap)
% axis image off;
% print(gcf,'PeakSep_131','-dpng','-r300')


%% SIMULATION COMPARISON %% - all of below is only for comparing measurement to simulations
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

Sphere.profile = @(x) 1.5*exp(-(acos(x)-90*degree).^2./(3*degree).^2) - ...
  exp(-(acos(x)-87*degree).^2./(2*degree).^2) - ...
  exp(-(acos(x)-93*degree).^2./(2*degree).^2);

% expand into a Legendre series
Sphere.profileHarm = S2Kernel.quadrature(Sphere.profile);

%Define all rotation matrices needed in the code
GMat_test=conv_EA_to_G([phi1,PHI,phi2]);
Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
Detector_tilt = Rx(tilt);
Sphere.rottoplot=GMat_test*Detector_tilt;

%second test orientation
GMat_test=conv_EA_to_G([2,0.8,1]);
rot2=GMat_test*Detector_tilt;
clear Detector_tilt GMat_test Rx

%Collect the band profiles
[Sphere.psi_full,Sphere.psi_pats,Sphere.psi_pats_Rot]=Spherical_Phases(h_set,InputUser,RTI,MapData,Refine.PC_out,crystals,Sphere,Sphere.rottoplot,'Experimental');

%% Calculation of the simulated profile
%generate the inputs to the spherical evaluation
x = linspace(-1,1,SamplingNum);
input = acos(x)./degree;

%check if there is only one symmetrically equivalent band

for phase=1:4
    for plane=1:4
        psi_pat=Sphere.psi_pats{phase,plane};

        if size(psi_pat,2) ~= 1
            for i=1:size(psi_pat,2) %and loop through the symmetric equivalents
                y(i,:) =(psi_pat{i}.eval(x));
            end
            y_av=mean(y,1); %then average
        else
            y_av=psi_pat.eval(x); %or just take the original
        end

        % Now have intensity profile as function of input angle, so can perform function

        %sample points
        sample=linspace(90-delta,90+delta,SamplingNum);
        I=interp1(input,y_av,sample);

        braggsum(phase,plane)=sum(I(sample>90-bragg_delta(1)&sample<90+bragg_delta(1)))*(sample(2)-sample(1));
    end
end







