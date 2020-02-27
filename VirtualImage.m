%% TPM 23/01/20
InputUser.BinFiles={'C:\Users\tpm416\Documents\GitHub\RTM_indexing\masterpatterns\Ni_1024.bin','C:\Users\tpm416\Documents\GitHub\RTM_indexing\masterpatterns\Ni3Al_1024.bin'}; 
InputUser.cifnames={'C:\Users\tpm416\Documents\GitHub\AstroEBSD\phases\Ni.cif','C:\Users\tpm416\Documents\GitHub\AstroEBSD\phases\Ni3Al.cif'};
InputUser.ResultsDir='C:\Users\tpm416\Documents\GitHub\BandAnalysis\Results';
InputUser.BA_Dir='C:\Users\tpm416\Documents\GitHub\BandAnalysis';

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
[Refine]=PC_refine(Eulers,Refine,PC_start,RefPatCor,MicroscopeData,RTI,'Ni');
pTime(['PC refined: ',num2str(100* Refine.Increase),'% increase'],t1);

clear Eulers

%% Get orientation matrix
GMat_test=conv_EA_to_G([Refine.Eulers_out(1),Refine.Eulers_out(2),Refine.Eulers_out(3)]);
Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
OriMatrix=GMat_test*Rx(MicroscopeData.TotalTilt);

cs=loadCIF('Ni');

%% Need to run this then cancel the first time...
[screen_int,~] = Cube_Generate(InputUser.BinFiles{1},0);
% as a function handle
master = S2FunHandle(@(v) Cube_Sample(v.x(:),v.y(:),v.z(:),screen_int,0));
% as hamonic expansion
masterHarm = S2FunHarmonicSym.quadrature(master,cs,'bandwidth',50); %run to initialise

%% Virtual imaging for a group of plane families and functions:
%first delta was 2.5
delta=5; %acquistion window

SamplingNum=10000;

%Define bands to operate on
h1=Miller({1,1,1},cs);
h2=Miller({1,1,0},cs);
h3=Miller({1,0,0},cs);
h4=Miller({1,3,1},cs);
h_set=[h1,h2,h3,h4];

lambda=6.626e-24/(2*9.109e-31*1.602e-19*20000)^0.5; % in angstroms

%% Get the profiles for a given set of bands for the full map
for testno=1:length(h)
disp(['Starting analysis ',num2str(testno),' of ',num2str(length(h_set))]);

%Grab the right band
h=h_set(testno);

%% Plot for one pattern and show the highlighted region
cd('C:\Users\tpm416\Documents\GitHub\BandAnalysis\Results\Results_yprime3')
mkdir('VirtualImages')
cd('VirtualImages')

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
        [Profile(i,j,:),sample]=AngularAperture(RefPat_cor,delta,SamplingNum,h,cs,Refine.PC_out,OriMatrix,0,lambda); %get the angular profile
        
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
%% Apply a specific functions

%Define the operating function and the dataset to be used
VIm_fn=@(I) sum(I);

Band={'111','110','100','131'};
bragg_delta=asin(lambda./(2.*h_set.dspacing))*180/pi; % in degrees
widths=bragg_delta./(5/10000); %number of points in the Bragg width

for bandno = 1:4
I=Results.(['Profile_',Band{bandno}]);
rows=size(I,1);
cols=size(I,2);

% Apply the defined function
VImage=zeros(rows,cols);
for i =1:rows
    for j=1:cols
        VImage(i,j)=VIm_fn(I(i,j,round(5000-widths(bandno)):round(5000+widths(bandno))));
    end
end

cmap=cbrewer('seq','YlOrRd',1000);
imagesc(VImage)
colormap(cmap)
axis image off;
print(gcf,['BraggSum_',Band{bandno}],'-dpng','-r300')
end

%% Map width of separation (used on 111)

I=Results.Profile_131;
n_max=2;

VImage=zeros(rows,cols);
for i=1:rows
    for j=1:cols
        
        [p,l]=findpeaks(squeeze(-I(i,j,:)),'SortStr','descend','NPeaks',n_max);
        diff=abs(l(2)-l(1));
        VImage(i,j)=diff; %save the locations
    
    end
end

cmap=cbrewer('seq','YlOrRd',1000);
imagesc(VImage)
colormap(cmap)
axis image off;
print(gcf,'PeakSep_131','-dpng','-r300')

%%
I=squeeze(Results.Profile_111(1,1,:));
L=length(I);
F=abs(fft(I)/L);
F1 = F(1:L/2+1);
F1(2:end-1) = 2*F1(2:end-1);




