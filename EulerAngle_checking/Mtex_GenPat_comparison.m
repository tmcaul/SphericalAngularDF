InputUser.Phases={'Ni','Ni3Al'};
BinFiles={'C:\Users\tpm416\Documents\GitHub\RTM_indexing\masterpatterns\Ni_1024.bin','C:\Users\tpm416\Documents\GitHub\RTM_indexing\masterpatterns\Ni3Al_1024.bin'}; 
cifnames={'C:\Users\tpm416\Documents\GitHub\AstroEBSD\phases\Ni.cif','C:\Users\tpm416\Documents\GitHub\AstroEBSD\phases\Ni3Al.cif'};
isHexes={0,0};

phase=1;

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

%%
InputUser.Phase_Input=InputUser.Phases(phase);
cs=loadCIF(InputUser.Phases{phase});

InputUser.HDF5_file='yprime3.h5';
InputUser.HDF5_folder='E:\Tom\GammaPrime_Data\V208C';
[ MapData,MicroscopeData,~,EBSPData ]=bReadHDF5( InputUser );
[Data_InputMap] = EBSD_Map(MapData,MicroscopeData);

%%
tilt=MicroscopeData.TotalTilt;

phi1=RTI.Phi1(1);
PHI=RTI.PHI(1);
phi2=RTI.Phi2(1);

PC_X = PC_refined(1);
PC_Y = PC_refined(2);
PC_Z = PC_refined(3);

%%

% get reference pattern
[ RefPat ] = bReadEBSP(EBSPData,1);
[ RefPatCor ] = EBSP_BGCor(RefPat,Settings_Cor);

[ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,num_Phases, RTI_info ] = Phase_Builder_RTI( InputUser.Phase_Input,RTI.Phase_Folder, RTI.Bin_loc );
[screen_int,facedata] = Cube_Generate(RTI_info.bin_file,RTI_info.isHex);

GMat_test=conv_EA_to_G([phi1,PHI,phi2]);

%Define all rotation matrices needed in the code
Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
Detector_tilt = Rx(tilt);
rottoplot=GMat_test*Detector_tilt;

%rottoplot

new_eulers=conv_G_to_EA(rottoplot);

PatternInfo.ScreenWidth=RTI.screensize;
PatternInfo.ScreenHeight=RTI.screensize;

%locs=1;
PC_pat=[PC_X,PC_Y,PC_Z];
[ EBSP_pat ] = EBSP_Gnom( PatternInfo,PC_pat);
[ pat ] = EBSP_gen( EBSP_pat,rottoplot,screen_int,RTI_info.isHex );

% Simulate MTEX pattern
%ori = orientation.byEuler(-phi1,-PHI,-phi2,'ZYZ',cs);
ori = orientation.byEuler(-new_eulers(1),-new_eulers(2),-new_eulers(3),'ZYZ',cs);
det = detector(RTI.screensize,RTI.screensize,PC_Z,[1-PC_Y,1-PC_X]);
[screen_int,facedata] = Cube_Generate(BinFiles{phase},0);
% as a function handle
master = S2FunHandle(@(v) Cube_Sample(v.x(:),v.y(:),v.z(:),screen_int,0));
pat_mtex = det.simulatePattern(master,(ori),1000,1);

% Plot both together
subplot(1,3,1)
imagesc(pat)
title('GenPat')
pbaspect([1,1,1])
subplot(1,3,2)
imagesc(pat_mtex)
title('MTEX')
pbaspect([1,1,1])
subplot(1,3,3)
imagesc(RefPatCor)
pbaspect([1,1,1])
title('Measured')

