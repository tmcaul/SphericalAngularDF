
function []=PlotPatterns(h,hnames,psi_pats,label,leg,ylimits,saving)

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
            SymmAvPlot(p,'DisplayName',char(h(i)),'linewidth',1,'LineStyle','-')
        else
            SymmAvPlot(p,'DisplayName',char(h(i)),'linewidth',1,'LineStyle','-')
        end
        hold on
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