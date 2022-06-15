% -----------------------------------------------------------------
%  graph_QoI_UQ.m
% ----------------------------------------------------------------- 
function fig = graph_QoI_UQ(time1,yraw1,time2,yraw2,samples,...
                            time,yopt,ybest,ymedian,ylow,yupp,graphobj)
    
    % check number of arguments
    if nargin < 12
        error('Too few inputs.')
    elseif nargin > 12
        error('Too many inputs.')
    end

    % check arguments
    if length(time1) ~= length(yraw1)
        error('time1 and yraw1 vectors must be same length')
    end
    
    if length(time2) ~= length(yraw2)
        error('time2 and yraw2 vectors must be same length')
    end

    if length(time) ~= length(yopt)
        error('time and yopt vectors must be same length')
    end
    
    if length(time) ~= length(ybest)
        error('time and ybest vectors must be same length')
    end
    
    if length(time) ~= length(ymedian)
        error('time and ymedian vectors must be same length')
    end
    
    if length(ylow) ~= length(yupp)
        error('ylow and yupp vectors must be same length')
    end
    
    if length(time) ~= length(ylow)
        error('time and ylow vectors must be same length')
    end
    
    % convert to row vectors so fliplr can work
    if find( size(time1) == max(size(time1)) ) < 2
        time1=time1';
    end
    
    % convert to row vectors so fliplr can work
    if find( size(time2) == max(size(time2)) ) < 2
        time2=time2';
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
    
    if find( size(yopt) == max(size(yopt)) ) < 2
        yopt=yopt';
    end
    
    if find( size(ybest) == max(size(ybest)) ) < 2
        ybest=ybest';
    end
    
    if find( size(ymedian) == max(size(ymedian)) ) < 2
        ymedian=ymedian';
    end
    
    if find( size(ylow) == max(size(ylow)) ) < 2
        ylow = ylow';
    end
    
    if find( size(yupp) == max(size(yupp)) ) < 2
        yupp = yupp';
    end
    
    % custom color
    MyGray  = [0.9 0.9 0.9];
    
    fig = figure('Name',graphobj.gname,'NumberTitle','off');
    
    fh7 = plot(time ,samples,'-'  );
    hold all
    fh6 = fill([time fliplr(time)],[yupp fliplr(ylow)],graphobj.colorT);
    fh5 = plot(time ,ymedian,'-'  );
    fh4 = plot(time ,ybest  ,'--' );
    fh3 = plot(time ,yopt   ,'-.' );
    fh2 = plot(time2,yraw2  ,'*--');
    fh1 = plot(time1,yraw1  ,'o--');
    
    set(gcf,'color'      ,'white'                           );
    set(gca,'position'   ,[0.2 0.2 0.7 0.7]                 );
    set(gca,'Box'        ,'on'                              );
    set(gca,'TickDir'    ,'out'     ,'TickLength',[.02 .02] );
    set(gca,'XMinorTick' ,'off'     ,'YMinorTick','off'     );
    set(gca,'XGrid'      ,'off'     ,'YGrid'     ,'off'     );
    set(gca,'XColor'     ,[.3 .3 .3],'YColor'    ,[.3 .3 .3]);
    set(gca,'FontName'   ,'Helvetica'                       );
    set(gca,'FontSize'   ,18                                );
    
    datetick('x',28,'keeplimits');
    
    set(fh1,'DisplayName',graphobj.leg1);
    set(fh2,'DisplayName',graphobj.leg2);
    set(fh3,'DisplayName',graphobj.leg3);
    set(fh4,'DisplayName',graphobj.leg4);
    set(fh5,'DisplayName',graphobj.leg5);
    set(fh6,'DisplayName',graphobj.leg6);
    set(fh7,'DisplayName',graphobj.leg7);
    leg         = [fh6; fh7(1); fh5; fh4; fh3; fh1; fh2];
    leg =  legend(leg  ,'Location','Best');
    %icons       = findobj(icons,'Type'    ,'Patch')
    %set(icons,'FaceAlpha',0.25)
    set(leg  ,'FontSize' ,16);

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
    
    hold off
    
    set(fh1,'Color','m'           ,'LineWidth',0.5)
    set(fh2,'Color','c'           ,'LineWidth',0.5)
    set(fh3,'Color',graphobj.color,'LineWidth',1)
    set(fh4,'Color',graphobj.color,'LineWidth',2)
    set(fh5,'Color',graphobj.color,'LineWidth',3)
    set(fh7,'Color',MyGray        ,'LineWidth',0.5)
    set(fh6,'EdgeColor','w')
    set(gca,'XTickLabelRotation',30)
    
	title(graphobj.gtitle,'FontSize',20,'FontName','Helvetica');
    
    if ( strcmp(graphobj.flag,'eps') )
        saveas(gcf,graphobj.gname,'epsc2');
    end

end
% -----------------------------------------------------------------