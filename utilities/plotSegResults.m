function plotSegResults(X, posGrandTruth, posSet, XtickStep)
%PLOTSEGRESULTS 
% This routine plots the multivariate time series segmentation results
% Inputs:
%   X - multivariate time series
%   posGrandTruth - the grand truth segmentation posistions
%   posSet - a Set of segmentation positions
%       posSet{i,1} - the segmentation posistions 
%       posSet{i,2} - the name of the algorithm
%   XtickStep - the step of Xtick

[n, m] = size(X);
xx = 1:n;
xxlim = n;
xxtick = XtickStep*(0:floor(n/XtickStep));
M = size(posSet, 1); % the number of algorithms

% color order: red, yellow, purple, lightblue, darkred, blue
% newcolors = {'#F00','#A0F','#0B0','#50F','#A0F','#F00'};
newcolors = {'red','c','m','k','b','y'};
% newcolors = {'b','c','m','b','k','y'};

figure;
for i=1:m
    % plot the i-th dimension
    subaxis(m+1,1,i, 'SpacingVert',0);
    plot(xx, X(:,i), 'LineWidth', 0.5)
    xlim([0 xxlim]);

    yy=get(gca,'ylim');
    
    % draw vertical lines at the grand truth segmentation positions
    for j=1:length(posGrandTruth)
        hold on
        plot([posGrandTruth(j) posGrandTruth(j)], yy, 'color', 'green', ...
            'LineWidth', 3)
    end
       
    % draw vertical lines at segmentation positions found by algorithms
    for I=1:M
        pos = posSet{I,1};
        for k=1:length(pos)
            hold on
            plot([pos(k) pos(k)], yy, 'color', newcolors{I}, 'LineWidth', 1)
        end
    end
    set(gca,'xtick',[])
    set(gca,'ytick',[])
end

xlabel('Time');
set(gca,'xtick',xxtick);
% sgtitle('Time series segmentation       ');

% Add legend
subaxis(m+1,1,m+1, 'SpacingVert',0);
h = zeros(M+1,1);
h(1) = plot(NaN, NaN, '-g');
for i=1:M
    hold on
    h(i+1) = plot(NaN, NaN, '-', 'Color', newcolors{i});
end
% legend(h, 'Grand truth','DiPCASeg');
labelArray = posSet(:,2);
labelArray = ['Grand Truth'; labelArray];
legend(h, labelArray);
axis off

end

