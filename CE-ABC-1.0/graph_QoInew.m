% -----------------------------------------------------------------
%  graph_QoInew.m
% ----------------------------------------------------------------- 
function fig = graph_QoInew(time,H,NewD,graphobj)
    
    % check number of arguments
    if nargin < 4
        error('Too few inputs.')
    elseif nargin > 4
        error('Too many inputs.')
    end

    % check arguments
    if length(time) ~= length(NewD)
        error('time and NewD vectors must be same length')
    end
    
    if length(time) ~= length(H)
        error('time and H vectors must be same length')
    end
    
    % convert to row vectors so fliplr can work
    if find( size(time) == max(size(time)) ) < 2
        time=time';
    end
    
    if find( size(NewD) == max(size(NewD)) ) < 2
        NewD=NewD';
    end
    
    if find( size(H) == max(size(H)) ) < 2
        H=H';
    end
    
    % custom color
    brown  = [101  33 33]/256;
    
    fig = figure('Name',graphobj.gname,'NumberTitle','off');
    
    figH    = plot(time,H   ,'DisplayName',graphobj.labelH   ,'Color',brown);
    hold on
    figNewD = plot(time,NewD,'DisplayName',graphobj.labelNewD,'Color','k');
    set(gcf,'color','white');
    set(gca,'position',[0.2 0.2 0.7 0.7]);
    set(gca,'Box','on');
    set(gca,'TickDir','out','TickLength',[.02 .02]);
    set(gca,'XMinorTick','off','YMinorTick','off');
    set(gca,'XGrid','off','YGrid','off');
    set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
    set(gca,'FontName','Helvetica');
    set(gca,'FontSize',18);
    
    leg = [figH; figNewD];
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
    
    set(figH   ,'LineWidth',2);
    set(figNewD,'LineWidth',2);
    
    hold off
	grid on
    
	title(graphobj.gtitle,'FontSize',20,'FontName','Helvetica');
    
    if ( strcmp(graphobj.flag,'eps') )
        saveas(gcf,graphobj.gname,'epsc2');
    end

end
% -----------------------------------------------------------------