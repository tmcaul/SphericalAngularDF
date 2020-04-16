%TPM
function [I,sample]=AngularAperture(pat,delta,Sampling,h,cs,PC,OriMatrix,plotting,lambda)

PC_X=PC(1);
PC_Y=PC(2);
PC_Z=PC(3);

new_eulers=conv_G_to_EA(OriMatrix);
ori = orientation.byEuler(-new_eulers(1),-new_eulers(2),-new_eulers(3),'ZYZ',cs);

det = detector(size(pat,1),size(pat,1),PC_Z,[1-PC_Y,1-PC_X]);

%approximate the measured pattern by harmonic function
pHarm = det.pattern2Fun(pat,'bandwidth',256,'quadrature','delta',0.1);

% %get radon transform
% RadonPHarm = conv(pHarm,profileHarm)./ correction;
% RadonPHarm.bandwidth = 64;

%Integrate along the profile
g=h.symmetrise;
for j = 1:length(g(1).symmetrise)  
[~,psi_pat{j}] = pHarm.symmetrise(ori*g(j));
psi_pat{j}.A(1) = 0;
end

% Calculation of the profile
%generate the inputs to the spherical evaluation
x = linspace(-1,1,Sampling);
input = acos(x)./degree;

%check if there is only one symmetrically equivalent band
if size(psi_pat,2) ~= 1
    for i=1:size(psi_pat,2) %and loop through the symmetric equivalents
        y(i,:) =(psi_pat{i}.eval(x));
    end
y_av=mean(y,1); %then average
else
y_av=psi_pat.eval(x); %or just take the original
end

%% Now have intensity profile as function of input angle, so can perform function

%sample points
sample=linspace(90-delta,90+delta,Sampling);
I=interp1(input,y_av,sample);

%% Plot the boundaries if required
if plotting==1
    
    bragg_delta=asin(lambda./(2.*h.dspacing))*180/pi; % in degrees
    width=bragg_delta./(delta/Sampling); %number of points in the Bragg width

    figure
    plot(pHarm,'pcolor','resolution',0.25*degree,'upper','complete')
    colormap('gray')
    circle(ori*h.symmetrise,'linecolor','r','LineWidth',1.5)
    print(gcf,['Bands_',char(h)],'-dpng','-r300')
    
    minval=-0.3;
    maxval=0.3;

    % Highlight a region
    figure
    plot(input,y_av,'LineWidth',2)
    xlim([80,100])
    ylim([-0.3,0.3])

    % h1 = line([90-delta 90-delta],[minval maxval]);
    % h2 = line([90+delta 90+delta],[minval maxval]);
    % set([h1 h2],'Color','k','LineWidth',0)
    patch([90-delta 90+delta 90+delta 90-delta],[-0.3 -0.3 0.3 0.3],'r','EdgeColor','None')
    alpha(0.5)
    set(gca,'children',flipud(get(gca,'children')))
    hold on
    lower=round(Sampling/2-width/2);
    upper=round(Sampling/2+width/2);
    plot([sample(lower);sample(upper)],[I(lower);I(upper)],'ko','MarkerSize',5,'MarkerFaceColor','k')
    
    
    print(gcf,['AngularAperture_',char(h),'_',strrep(num2str(delta),'.',',')],'-dpng','-r300')
end
    
end