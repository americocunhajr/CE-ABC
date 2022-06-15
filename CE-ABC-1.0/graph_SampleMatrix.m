% -----------------------------------------------------------------
%  graph_SampleMatrix.m
% ----------------------------------------------------------------- 
function fig = graph_SampleMatrix(ABCobj,graphobj)
    
    % check number of arguments
    if nargin < 2
        error('Too few inputs.')
    elseif nargin > 2
        error('Too many inputs.')
    end
    
    % custom color
    MyPurple = [255 0 127]/256;
    
    fig = figure('Name',graphobj.gname,'NumberTitle','off');
    
    [fh1,ax] = gplotmatrix([ABCobj.x_accept(1,:);...
                            ABCobj.x_accept(2,:);...
                            ABCobj.x_accept(3,:);...
                            ABCobj.x_accept(4,:);...
                            ABCobj.x_accept(5,:);...
                            ABCobj.x_accept(6,:);...
                            ABCobj.x_accept(7,:);...
                            ABCobj.x_accept(8,:);...
                            ABCobj.x_accept(9,:);...
                            ABCobj.x_accept(10,:);...
                            ABCobj.x_accept(11,:);...
                            ABCobj.x_accept(12,:)]',...
                            [],[],MyPurple);
    set(gcf,'color','white');
    set(gca,'position',[0.2 0.2 0.7 0.7]);
    set(gca,'Box','on');
    set(gca,'TickDir','out','TickLength',[.02 .02]);
    set(gca,'XMinorTick','off','YMinorTick','off');
    set(gca,'XGrid','off','YGrid','off');
    set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
    set(gca,'FontName','Helvetica');
    set(gca,'FontSize',10);
    set(ax,'XTick',[]);
    set(ax,'YTick',[]);
    set(ax,'Xticklabel',[]);
    set(ax,'Yticklabel',[]) ;

%     if ( strcmp(graphobj.xmin,'auto') || strcmp(graphobj.xmax,'auto') )
%         xlim('auto');
%     else
%         xlim([graphobj.xmin graphobj.xmax]);
%     end
%     
%     if ( strcmp(graphobj.ymin,'auto') || strcmp(graphobj.ymax,'auto') )
%         ylim('auto');
%     else
%         ylim([graphobj.ymin graphobj.ymax]);
%     end
    
	title(graphobj.gtitle,'FontSize',20,'FontName','Helvetica');
    
    if ( strcmp(graphobj.flag,'eps') )
        saveas(gcf,graphobj.gname,'epsc2');
    end

end
% -----------------------------------------------------------------