InputUser.HDF5_folder='E:\Tom\GammaPrime_Data\V208C\';
InputUser.HDF5_file='yprime3.h5';
[ MapData,MicroscopeData,PhaseData,EBSPData ]=bReadHDF5( InputUser );
[Data_InputMap] = EBSD_Map(MapData,MicroscopeData);

%%Map Settings
pat=RefPatCor;
delta=2.5;
cs;
h=Miller({1,0,0},cs);
PC=Refine.PC_out;
OriMatrix=Sphere.rottoplot;

%define the function you want to operate with (eg: maximum)
VIm_fn=@(I) max(I);

%% Plot for one pattern and show the highlighted region
[ RefPat ] = bReadEBSP(EBSPData,1);
[ RefPat_cor] = EBSP_BGCor( RefPat,Settings_Cor );
[I,sample]=AngularAperture(pat,delta,h,cs,PC,OriMatrix,1);

%% Load the data and do the angular aperturing, then apply the defined function (SERIES)
% rows=Data_InputMap.ypts;cols=Data_InputMap.xpts;
% VImage=zeros(rows,cols);
% wb=waitbar(0,'Generating virtual image...');
% for i=1:rows
%     for j=1:cols
%         pattern_number=Data_InputMap.PMap(i,j);
%         [ RefPat ] = bReadEBSP(EBSPData,i);
%         [ RefPat_cor] = EBSP_BGCor( RefPat,Settings_Cor );
%         [I,sample]=AngularAperture(RefPat_cor,delta,h,cs,PC,OriMatrix,0); %get the angular profile
%         VImage(i,j)=VIm_fn(I); %apply the function
%         
%         %update the waitbar
%         waitbar(((i-1)*cols+j)/(rows*cols),wb,['Evaluating pattern ',num2str(((i-1)*cols+j)),' of ',num2str(rows*cols)])
%         
%     end
% end

%% Load the data and do the angular aperturing, then apply the defined function (PARALLEL)
rows=Data_InputMap.ypts;cols=Data_InputMap.xpts;
VImage=zeros(rows,cols);
t1=clock;

pTime('Starting virtual image calculation...',t1);
for i=1:rows
    parfor j=1:cols
        pattern_number=Data_InputMap.PMap(i,j);
        [ RefPat ] = bReadEBSP(EBSPData,i);
        [ RefPat_cor] = EBSP_BGCor( RefPat,Settings_Cor );
        [I,sample]=AngularAperture(RefPat_cor,delta,h,cs,PC,OriMatrix,0); %get the angular profile
        VImage(i,j)=VIm_fn(I); %apply the function
    end
    pTime(['Completed row ',num2str(i),' of ',num2str(rows)],t1);
end


