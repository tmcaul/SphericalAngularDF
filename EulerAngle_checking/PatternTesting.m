
setMTEXpref('xAxisDirection','west');
setMTEXpref('zAxisDirection','outOfPlane');

det = detector(RTI.screensize,RTI.screensize,MapData.DD(1),[0.5,0.5]);
%det = detector(RTI.screensize,RTI.screensize,MapData.DD(1),[0,0]);

%[det2,pattern]=GenPat(InputUser,RTI,1,Data_InputMap,MicroscopeData);

% define a cut off function for the detector
mask = det.S2CutOffMask(0);
maskHarm = S2FunHarmonic.quadrature(mask);
correction = max(maskHarm.radon,0.1);

%eulers_RTImtex=[4.3034,0.7404,2.4879];
%eulers_RTI=[0.3524,1.0007,0.5073];
%eulers_bruker=[-115.1227*degree,43.1124*degree,144.5219*degree];

eulers_b2=[245.6*degree,43*degree,144.1*degree];
%eulers_b2=[144.1*degree,43*degree,245.6*degree];
e=eulers_b2;
GMat_test=conv_EA_to_G(e);
Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
Detector_tilt = Rx(MicroscopeData.TotalTilt);
rottoplot=(GMat_test*Detector_tilt);
ori=orientation.byMatrix(rottoplot,cs,cs);

eulers_empirical=[41*degree,-22*degree,2*degree];
%EQUIV to [221,22,182]
e2=(eulers_empirical);
GMat_empirical=conv_EA_to_G(e2);
ori2 = orientation('Euler',e2);

pattern_empirical = det.simulatePattern(master,(ori2),1000,1);

pattern_known = det.simulatePattern(master,(ori),1000,1);

figure
subplot(1,2,1)
imagesc(pattern_empirical)
pbaspect([1,1,1])
subplot(1,2,2)
imagesc(pattern_known)
pbaspect([1,1,1])