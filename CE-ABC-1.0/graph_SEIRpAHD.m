% -----------------------------------------------------------------
%  graph_SEIRpAHD.m
% ----------------------------------------------------------------- 
function fig = graph_SEIRpAHD(time,S,E,I,R,A,H,D,N,graphobj)
    
    % check number of arguments
    if nargin < 10
        error('Too few inputs.')
    elseif nargin > 10
        error('Too many inputs.')
    end

    % check arguments
    if length(time) ~= length(S)
        error('time and S vectors must be same length')
    end

    if length(time) ~= length(E)
        error('time and E vectors must be same length')
    end
    
    if length(time) ~= length(I)
        error('time and I vectors must be same length')
    end
    
    if length(time) ~= length(R)
        error('time and R vectors must be same length')
    end
    
    if length(time) ~= length(A)
        error('time and A vectors must be same length')
    end
    
    if length(time) ~= length(H)
        error('time and H vectors must be same length')
    end
    
    if length(time) ~= length(D)
        error('time and D vectors must be same length')
    end

    if length(time) ~= length(N)
        error('time and N vectors must be same length')
    end
    
    % convert to row vectors so fliplr can work
    if find( size(time) == max(size(time)) ) < 2
        time=time';
    end
    
    if find( size(S) == max(size(S)) ) < 2
        S=S';
    end
    
    if find( size(E) == max(size(E)) ) < 2
        E=E';
    end
    
    if find( size(I) == max(size(I)) ) < 2
        I=I';
    end
    
    if find( size(R) == max(size(R)) ) < 2
        R=R';
    end
    
    if find( size(A) == max(size(A)) ) < 2
        A=A';
    end
    
    if find( size(H) == max(size(H)) ) < 2
        H=H';
    end
    
    if find( size(D) == max(size(D)) ) < 2
        D=D';
    end
    
    if find( size(N) == max(size(N)) ) < 2
        N=N';
    end
    
    % custom color
    yellow = [255 204  0]/256;
    orange = [256 128  0]/256;
    brown  = [101  33 33]/256;
    
    fig = figure('Name',graphobj.gname,'NumberTitle','off');
    
    figS = plot(time,S,'DisplayName',graphobj.labelS,'Color','b');
    hold on
    figE = plot(time,E,'DisplayName',graphobj.labelE,'Color',yellow);
    figI = plot(time,I,'DisplayName',graphobj.labelI,'Color','r');
    figR = plot(time,R,'DisplayName',graphobj.labelR,'Color','g');
    figA = plot(time,A,'DisplayName',graphobj.labelA,'Color',orange);
    figH = plot(time,H,'DisplayName',graphobj.labelH,'Color',brown);
    figD = plot(time,D,'DisplayName',graphobj.labelD,'Color','k');
    figN = plot(time,N,'DisplayName',graphobj.labelN,'Color','m');
    set(gcf,'color','white');
    set(gca,'position',[0.2 0.2 0.7 0.7]);
    set(gca,'Box','on');
    set(gca,'TickDir','out','TickLength',[.02 .02]);
    set(gca,'XMinorTick','off','YMinorTick','off');
    set(gca,'XGrid','off','YGrid','off');
    set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
    set(gca,'FontName','Helvetica');
    set(gca,'FontSize',18);
    
    leg = [figS; figE; figI; figR; figA; figH; figD; figN];
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
    
    set(figS,'LineWidth',2);
    set(figE,'LineWidth',2);
    set(figI,'LineWidth',2);
    set(figR,'LineWidth',2);
    set(figA,'LineWidth',2);
    set(figH,'LineWidth',2);
    set(figD,'LineWidth',2);
    set(figN,'LineWidth',2);
    
    hold off
	grid on
    
	title(graphobj.gtitle,'FontSize',20,'FontName','Helvetica');
    
    if ( strcmp(graphobj.flag,'eps') )
        saveas(gcf,graphobj.gname,'epsc2');
    end

end
% -----------------------------------------------------------------