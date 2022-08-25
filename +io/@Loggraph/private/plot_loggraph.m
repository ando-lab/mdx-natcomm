function plot_loggraph(LogTable,graphNumber)

cols = LogTable.Graphs(graphNumber).cols;

xdata = LogTable.data(:,cols(1));
ydata = LogTable.data(:,cols(2:end));

axisLimits = LogTable.Graphs(graphNumber).axisLimits;
axisMode = LogTable.Graphs(graphNumber).axisMode;
graphName = LogTable.Graphs(graphNumber).name;
curveNames = LogTable.columns(cols(2:end));
xAxisName = LogTable.columns{cols(1)};

switch LogTable.type
    case {'graphs'}
        plotMarker = '-';
    case {'scatter'}
        plotMarker = 'o';
    otherwise
        plotMarker = '-'; % default: use line plot
end
    
    plot(xdata,ydata,plotMarker);
    
    legend(curveNames,'Interpreter','None');
    xlabel(xAxisName,'Interpreter','None');
    title(graphName,'Interpreter','None');
    ax = gca;

% set axes limits
switch axisMode
    case {'auto',''}
        % do nothing
    case 'nought'
        ax.YLim(1) = 0;
    case 'manual'
        try
        ax.XLim = axisLimits(1:2);
        ax.YLim = axisLimits(3:4);
        catch
            % sometimes this fails...
        end
end


end
