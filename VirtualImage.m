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
for testno=1:12
disp(['Starting test ',num2str(testno),' of 12']);
%% Image specific settings

delta=2.5; %search window

%Choose which band to operate on
switch testno
    case {1,2,3}
        h=Miller({1,1,1},cs); %2h 40min each
    case {4,5,6}
        h=Miller({1,1,0},cs); %3h 56min each
    case {7,8,9}
        h=Miller({1,0,0},cs); %2h each
    case {10,11,12}
        h=Miller({1,3,1},cs);
end

%Choose the operating function
switch testno
    case {1,4,7,10}
        VIm_fn=@(I) max(I);
        label='max';
    case {2,5,8,11}
        VIm_fn=@(I) range(I);
        label='range';
    case {3,6,9,12}
        VIm_fn=@(I) sum(I);
        label='sum';
end

%% Plot for one pattern and show the highlighted region
cd('C:\Users\tpm416\Documents\GitHub\BandAnalysis\Results\Results_yprime3')
mkdir('VirtualImages')
cd('VirtualImages')

[ RefPat ] = bReadEBSP(EBSPData,1);
[ RefPat_cor] = EBSP_BGCor( RefPat,Settings_Cor );
[I,sample]=AngularAperture(RefPat_cor,delta,h,cs,Refine.PC_out,OriMatrix,1);

%% Load the data and do the angular aperturing, then apply the defined function (PARALLEL)
rows=Data_InputMap.ypts;cols=Data_InputMap.xpts;
VImage=zeros(rows,cols);
t1=clock;

pTime('Starting virtual image calculation...',t1);
for i=1:rows
    parfor j=1:cols
        pattern_number=Data_InputMap.PMap(i,j);
        [ RefPat ] = bReadEBSP(EBSPData,pattern_number);
        [ RefPat_cor] = EBSP_BGCor( RefPat,Settings_Cor );
        [I,sample]=AngularAperture(RefPat_cor,delta,h,cs,Refine.PC_out,OriMatrix,0); %get the angular profile
        VImage(i,j)=VIm_fn(I); %apply the function
    end
    pTime(['Completed row ',num2str(i),' of ',num2str(rows)],t1);
end

figure
imagesc(VImage)

band=char(h);
Results.(['Image_',band(2:end-1),'_',label])=VImage;
clear VImage band
end
%Ims.[]=VImage;

%% Single virtual image
cd('C:\Users\tpm416\Documents\GitHub\BandAnalysis\Results\Results_yprime3')
mkdir('VirtualImages')
cd('VirtualImages')

h=Miller({1,1,1},cs); %Band to operate on
VIm_fn=@(I) max(I); %Function
delta=2.5; %search window

[ RefPat ] = bReadEBSP(EBSPData,1);
[ RefPat_cor] = EBSP_BGCor( RefPat,Settings_Cor );
[I,sample]=AngularAperture(RefPat_cor,delta,h,cs,Refine.PC_out,OriMatrix,1);

% Load the data and do the angular aperturing, then apply the defined function (PARALLEL)
rows=Data_InputMap.ypts;cols=Data_InputMap.xpts;
VImage=zeros(rows,cols);
t1=clock;

pTime('Starting virtual image calculation...',t1);
for i=1:rows
    parfor j=1:cols
        pattern_number=Data_InputMap.PMap(i,j);
        [ RefPat ] = bReadEBSP(EBSPData,pattern_number);
        [ RefPat_cor] = EBSP_BGCor( RefPat,Settings_Cor );
        [I,sample]=AngularAperture(RefPat_cor,delta,h,cs,Refine.PC_out,OriMatrix,0); %get the angular profile
        VImage(i,j)=VIm_fn(I); %apply the function
    end
    pTime(['Completed row ',num2str(i),' of ',num2str(rows)],t1);
end

figure
imagesc(VImage)


