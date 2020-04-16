%TPM 2020
function []=PlotPatterns(h,psi_pats,label,leg,ylimits,saving,varargin)

if size(varargin,1)>0
    Settings=varargin{:};
    LS=Settings.LS;
    LC=Settings.LC;
end

%Grab hnames
for i=1:length(h)
    bandname=char(h(i));
    bandname=bandname(2:end-1);
    hnames{i}=['{',bandname,'}'];
end

% Plot band profiles
%for single, unrotated pat, rotated planes
for i = 1:length(h)
    %choose which band in the symmetrised set to plot
%     if i~=4
%         a=1;
%     else
%         a=3;
%     end
    
    figure
    title(hnames{i})
    
    for phase = 1:size(psi_pats,1)
        
        p=psi_pats{phase,i};
        
        if iscell(p)
            SymmAvPlot(p,'DisplayName',char(h(i)),'linewidth',1.5,'LineStyle','-')
        else
            SymmAvPlot(p,'DisplayName',char(h(i)),'linewidth',1.5,'LineStyle','-')
        end
        hold on
        
        if exist('LS')==1
            line = findobj(gcf, 'type', 'line');
            set(line(1),'LineStyle',LS{phase});
        end
        
        if exist('LC')==1
            line = findobj(gcf, 'type', 'line');
            set(line(1),'Color',LC{phase});
        end
        
    end
    
    legend(leg)
    xlim([80,100])
    xlabel('Band opening angle')
    ylim(ylimits)
    ylabel('Normalised profile intensity')
    hold off
    title(hnames{i})
    
    %cd(resultsfolder)
    if saving==1
    print([label,hnames{i}],'-dpng','-r300')
    end
    %cd(homefol) 
end