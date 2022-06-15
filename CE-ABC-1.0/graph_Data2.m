% -----------------------------------------------------------------
%  graph_Data2.m
% ----------------------------------------------------------------- 
function fig = graph_Data2(time,yraw1,yraw2,graphobj)
    
    % check number of arguments
    if nargin < 4
        error('Too few inputs.')
    elseif nargin > 4
        error('Too many inputs.')
    end

    % check arguments
    if length(time) ~= length(yraw1)
        error('time and yraw1 vectors must be same length')
    end
    
    if length(time) ~= length(yraw2)
        error('time and yraw2 vectors must be same length')
    end
    
    % convert to row vectors so fliplr can work
    if find( size(time) == max(size(time)) ) < 2
        time=time';
    end
    
    if find( size(yraw1) == max(size(yraw1)) ) < 2
        yraw1=yraw1';
    end
    
    if find( size(yraw2) == max(size(yraw2)) ) < 2
        yraw2=yraw2';
    end
    
    fig = figure('Name',graphobj.gname,'NumberTitle','off');
    
    fh1 = plot(time,yraw1,'o--','Color',graphobj.color1);
    hold on
    fh2 = plot(time,yraw2,'*--','Color',graphobj.color2);
    set(gcf,'color','white');
    set(gca,'position',[0.2 0.2 0.7 0.7]);
    set(gca,'Box','on');
    set(gca,'TickDir','out','TickLength',[.02 .02]);
    set(gca,'XMinorTick','off','YMinorTick','off');
    set(gca,'XGrid','off','YGrid','off');
    set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
    set(gca,'FontName','Helvetica');
    set(gca,'FontSize',18);
    datetick('x',28,'keeplimits');
    
    leg = legend(graphobj.leg1,graphobj.leg2,'Location','Best');
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
    
    grid on;
    
    set(fh1,'LineWidth',0.5)
    set(fh2,'LineWidth',0.5)
    set(gca,'XTickLabelRotation',30)
    
    hold off
    
	title(graphobj.gtitle,'FontSize',20,'FontName','Helvetica');
    
    if ( strcmp(graphobj.flag,'eps') )
        saveas(gcf,graphobj.gname,'epsc2');
    end

end
% -----------------------------------------------------------------