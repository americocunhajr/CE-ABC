% -----------------------------------------------------------------
%  graph_QoIcum.m
% ----------------------------------------------------------------- 
function fig = graph_QoIcum(time,CumH,D,graphobj)
    
    % check number of arguments
    if nargin < 4
        error('Too few inputs.')
    elseif nargin > 4
        error('Too many inputs.')
    end

    % check arguments
    if length(time) ~= length(CumH)
        error('time and CumH vectors must be same length')
    end
    
    if length(time) ~= length(D)
        error('time and D vectors must be same length')
    end
    
    % convert to row vectors so fliplr can work
    if find( size(time) == max(size(time)) ) < 2
        time=time';
    end
    
    if find( size(CumH) == max(size(CumH)) ) < 2
        CumH=CumH';
    end
    
    if find( size(D) == max(size(D)) ) < 2
        D=D';
    end
    
    % custom color
    brown  = [101  33 33]/256;
    
    fig = figure('Name',graphobj.gname,'NumberTitle','off');
    
    figCumH = plot(time,CumH,'DisplayName',graphobj.labelCumH,'Color',brown);
    hold on
    figD    = plot(time,D   ,'DisplayName',graphobj.labelD   ,'Color','k');
    set(gcf,'color','white');
    set(gca,'position',[0.2 0.2 0.7 0.7]);
    set(gca,'Box','on');
    set(gca,'TickDir','out','TickLength',[.02 .02]);
    set(gca,'XMinorTick','off','YMinorTick','off');
    set(gca,'XGrid','off','YGrid','off');
    set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
    set(gca,'FontName','Helvetica');
    set(gca,'FontSize',18);
    
    leg = [figCumH; figD];
    leg = legend(leg,'Location','Best');
    set(leg,'FontSize',16);

    if ( strcmp(graphobj.xmin,'auto') || strcmp(graphobj.xmax,'auto') )
        xlim('auto');
    else
        xlim([graphobj.xmin graphobj.xmax]);
    end
    
    if ( strcmp(graphobj.ymin,'auto') || strcmp(graphobj.ymax,'auto') )
        ylim('auto');
    else
        ylim([graphobj.ymin graphobj.ymax]);
    end
    
    labX = xlabel(graphobj.xlab,'FontSize',20,'FontName','Helvetica');
    labY = ylabel(graphobj.ylab,'FontSize',20,'FontName','Helvetica');
    
    set(figCumH,'LineWidth',2);
    set(figD   ,'LineWidth',2);
    
    hold off
	grid on
    
	title(graphobj.gtitle,'FontSize',20,'FontName','Helvetica');
    
    if ( strcmp(graphobj.flag,'eps') )
        saveas(gcf,graphobj.gname,'epsc2');
    end

end
% -----------------------------------------------------------------