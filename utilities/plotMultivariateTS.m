function plotMultivariateTS(X, posGrandTruth, XtickStep)
%plotMultivariateTS
% This routine plots the multivariate time series with ground truth
% segmentation locations.
%
% Inputs:
%   X - multivariate time series, each row is an observation
%   posGrandTruth - the grand truth segmentation posistions
%   XtickStep - the step of Xtick

[n, m] = size(X);
xx = 1:n;
xxlim = n;
xxtick = XtickStep*(0:floor(n/XtickStep));

newcolors = {'red','c','m','k','b','y'};

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
       
    set(gca,'xtick',[])
    set(gca,'ytick',[])
end

xlabel('Time');
set(gca,'xtick',xxtick);
% sgtitle('Time series segmentation       ');

end

