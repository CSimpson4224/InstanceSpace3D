function scriptfcn
% -------------------------------------------------------------------------
% scriptfcn.m
% -------------------------------------------------------------------------
%
% By: Mario Andres Munoz Acosta
%     School of Mathematics and Statistics
%     The University of Melbourne
%     Australia
%     2020
%
% -------------------------------------------------------------------------

writeArray2CSV = @(data,colnames,rownames,filename) writetable(array2table(data,'VariableNames',colnames,...
                                                                                'RowNames',rownames),...
                                                               filename,'WriteRowNames',true);
writeCell2CSV = @(data,colnames,rownames,filename) writetable(cell2table(data,'VariableNames',colnames,...
                                                                              'RowNames',rownames),...
                                                              filename,'WriteRowNames',true);
makeBndLabels = @(data) arrayfun(@(x) strcat('bnd_pnt_',num2str(x)),1:size(data,1),'UniformOutput',false);
colorscale  = @(data) round(255.*bsxfun(@rdivide, bsxfun(@minus, data, min(data,[],1)), range(data)));
colorscaleg = @(data) round(255.*bsxfun(@rdivide, bsxfun(@minus, data, min(data(:))), range(data(:))));

assignin('caller','writeArray2CSV',writeArray2CSV);
assignin('caller','writeCell2CSV',writeCell2CSV);
assignin('caller','makeBndLabels',makeBndLabels);
assignin('caller','colorscale',colorscale);
assignin('caller','colorscaleg',colorscaleg);
assignin('caller','drawSources',@drawSources);
assignin('caller','drawScatter',@drawScatter);
assignin('caller','drawPortfolioSelections',@drawPortfolioSelections);
assignin('caller','drawPortfolioFootprint',@drawPortfolioFootprint);
assignin('caller','drawGoodBadFootprint',@drawGoodBadFootprint);
assignin('caller','drawFootprint',@drawFootprint);
assignin('caller','drawBinaryPerformance',@drawBinaryPerformance);

end
% =========================================================================
% SUBFUNCTIONS
% =========================================================================
function handle = drawSources(Z, S)

ubound = ceil(max(Z));
lbound = floor(min(Z));
sourcelabels = cellstr(unique(S));
nsources = length(sourcelabels);
clrs = flipud(lines(nsources));
handle = zeros(nsources,1);
if size(Z,2) == 3
    for i=nsources:-1:1
        line(Z(S==sourcelabels{i},1), ...
             Z(S==sourcelabels{i},2), ...
             Z(S==sourcelabels{i},3), ...
             'LineStyle', 'none', ...
             'Marker', '.', ...
             'Color', clrs(i,:), ...
             'MarkerFaceColor', clrs(i,:), ...
             'MarkerSize', 4);
        handle(i) = patch([0 0],[0 0],[0 0], clrs(i,:), 'EdgeColor','none');
    end
    view([1,1,1])
    grid on
    xlabel('z_{1}'); ylabel('z_{2}'); zlabel('z_{3}'); title('Sources');
else
    for i=nsources:-1:1
        line(Z(S==sourcelabels{i},1), ...
             Z(S==sourcelabels{i},2), ...
             'LineStyle', 'none', ...
             'Marker', '.', ...
             'Color', clrs(i,:), ...
             'MarkerFaceColor', clrs(i,:), ...
             'MarkerSize', 4);
        handle(i) = patch([0 0],[0 0], clrs(i,:), 'EdgeColor','none');
    end
    xlabel('z_{1}'); ylabel('z_{2}'); title('Sources');
end
legend(handle, sourcelabels, 'Location', 'NorthEastOutside');
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square; axis([lbound(1)-1 ubound(1)+1 lbound(2)-1 ubound(2)+1]);

end
% =========================================================================
function handle = drawScatter(Z, X, titlelabel)

if size(Z,2) == 2 

    ubound = ceil(max(Z));
    lbound = floor(min(Z));
    handle = scatter(Z(:,1), Z(:,2), 8, X, 'filled');
    caxis([0,1])
    xlabel('z_{1}'); ylabel('z_{2}'); title(titlelabel);
    set(findall(gcf,'-property','FontSize'),'FontSize',12);
    set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
    axis square; axis([lbound(1)-1 ubound(1)+1 lbound(2)-1 ubound(2)+1]);
    colorbar('EastOutside');
    
elseif size(Z,2) == 3 
    
    ubound = ceil(max(Z));
    lbound = floor(min(Z));
    handle = scatter3(Z(:,1),Z(:,2),Z(:,3), 4, X, 'filled');
    caxis([0,1])
    view([1,1,1])
    grid on
    xlabel('z_{1}'); ylabel('z_{2}'); zlabel('z_{3}'); title(titlelabel);
    set(findall(gcf,'-property','FontSize'),'FontSize',12);
    set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
    axis square; %axis([lbound(1)-1 ubound(1)+1 lbound(2)-1 ubound(2)+1]);
    colorbar('EastOutside');
    
end

end
% =========================================================================
function drawPortfolioSelections(Z, P, algolabels, titlelabel)

ubound = ceil(max(Z));
lbound = floor(min(Z));
nalgos = length(algolabels);
algolbls = cell(1,nalgos+1);
h = zeros(1,nalgos+1);
isworthy = sum(bsxfun(@eq, P, 0:nalgos))~=0;
clr = flipud(lines(nalgos+1));
if size(Z,2) == 3
    for i=0:nalgos
        if ~isworthy(i+1)
            continue;
        end
        line(Z(P==i,1), Z(P==i,2), Z(P==i,3), 'LineStyle', 'none', ...
                                   'Marker', '.', ...
                                   'Color', clr(i+1,:), ...
                                   'MarkerFaceColor', clr(i+1,:), ...
                                   'MarkerSize', 4);
        h(i+1) = patch([0 0],[0 0],[0 0], clr(i+1,:), 'EdgeColor','none');
        if i==0
            algolbls{i+1} = 'None';
        else
            algolbls{i+1} = strrep(algolabels{i},'_',' ');
        end
    end
    xlabel('z_{1}'); ylabel('z_{2}'); zlabel('z_{3}'); title(titlelabel);
else
    for i=0:nalgos
        if ~isworthy(i+1)
            continue;
        end
        line(Z(P==i,1), Z(P==i,2), 'LineStyle', 'none', ...
                                   'Marker', '.', ...
                                   'Color', clr(i+1,:), ...
                                   'MarkerFaceColor', clr(i+1,:), ...
                                   'MarkerSize', 4);
        h(i+1) = patch([0 0],[0 0], clr(i+1,:), 'EdgeColor','none');
        if i==0
            algolbls{i+1} = 'None';
        else
            algolbls{i+1} = strrep(algolabels{i},'_',' ');
        end
    end
    xlabel('z_{1}'); ylabel('z_{2}'); title(titlelabel);
end
legend(h(isworthy), algolbls(isworthy), 'Location', 'NorthEastOutside');
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square; axis([lbound(1)-1 ubound(1)+1 lbound(2)-1 ubound(2)+1]);

end
% =========================================================================
function h = drawPortfolioFootprint(Z, best, P, algolabels)

% Color definitions
ubound = ceil(max(Z));
lbound = floor(min(Z));
nalgos = length(algolabels);
algolbls = cell(1,nalgos+1);
isworthy = sum(bsxfun(@eq, P, 0:nalgos))~=0;
clr = flipud(lines(nalgos+1));
h = zeros(1,nalgos+1);
if size(Z,2) == 3
    for i=0:nalgos
        if ~isworthy(i+1)
            continue;
        end
        line(Z(P==i,1), Z(P==i,2), Z(P==i,3), 'LineStyle', 'none', ...
                                   'Marker', '.', ...
                                   'Color', clr(i+1,:), ...
                                   'MarkerFaceColor', clr(i+1,:), ...
                                   'MarkerSize', 4);
        h(i+1) = patch([0 0],[0 0],[0 0], clr(i+1,:), 'EdgeColor','none');
        if i==0
            algolbls{i+1} = 'None';
        else
            drawFootprint(best{i}, clr(i+1,:), 0.3);
            algolbls{i+1} = strrep(algolabels{i},'_',' ');
        end
    end
    view([1,1,1])
    grid on
    xlabel('z_{1}'); ylabel('z_{2}'); zlabel('z_{3}'); title('Portfolio footprints');
else
    for i=0:nalgos
        if ~isworthy(i+1)
            continue;
        end
        line(Z(P==i,1), Z(P==i,2), 'LineStyle', 'none', ...
                                   'Marker', '.', ...
                                   'Color', clr(i+1,:), ...
                                   'MarkerFaceColor', clr(i+1,:), ...
                                   'MarkerSize', 4);
        h(i+1) = patch([0 0],[0 0], clr(i+1,:), 'EdgeColor','none');
        if i==0
            algolbls{i+1} = 'None';
        else
            drawFootprint(best{i}, clr(i+1,:), 0.3);
            algolbls{i+1} = strrep(algolabels{i},'_',' ');
        end
    end
    xlabel('z_{1}'); ylabel('z_{2}'); title('Portfolio footprints');
end
legend(h(isworthy), algolbls(isworthy), 'Location', 'NorthEastOutside');
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square; axis([lbound(1)-1 ubound(1)+1 lbound(2)-1 ubound(2)+1]);

end
% =========================================================================
function h = drawGoodBadFootprint(Z, good, Ybin, titlelabel,viewangle)

ubound = ceil(max(Z));
lbound = floor(min(Z));
orange = [1.0 0.6471 0.0];
blue = [0.0 0.0 1.0];
lbls = {'GOOD','BAD'};
h = zeros(1,2);
if size(Z,2) == 3
    if any(~Ybin)
        % drawFootprint(bad, orange, 0.2);
        line(Z(~Ybin,1), Z(~Ybin,2),Z(~Ybin,3), 'LineStyle', 'none', ...
                                     'Marker', '.', ...
                                     'Color', orange, ...
                                     'MarkerFaceColor', orange, ...
                                     'MarkerSize', 4);
        h(2) = patch([0 0],[0 0],[0 0], orange, 'EdgeColor','none');
    end
    if any(Ybin)
        line(Z(Ybin,1), Z(Ybin,2), Z(Ybin,3), 'LineStyle', 'none', ...
                                   'Marker', '.', ...
                                   'Color', blue, ...
                                   'MarkerFaceColor', blue, ...
                                   'MarkerSize', 4);
        h(1) = patch([0 0],[0 0],[0 0], blue, 'EdgeColor','none');
        drawFootprint(good, blue, 0.3);
    end
    if ~exist('viewangle','var')
        viewangle = [1,1,1];
    end 
    view(viewangle)
    grid on
    xlabel('z_{1}'); ylabel('z_{2}');zlabel('z_{3}'); title([titlelabel ' Footprints']);
else
    if any(~Ybin)
        % drawFootprint(bad, orange, 0.2);
        line(Z(~Ybin,1), Z(~Ybin,2), 'LineStyle', 'none', ...
                                     'Marker', '.', ...
                                     'Color', orange, ...
                                     'MarkerFaceColor', orange, ...
                                     'MarkerSize', 4);
        h(2) = patch([0 0],[0 0], orange, 'EdgeColor','none');
    end
    if any(Ybin)
        line(Z(Ybin,1), Z(Ybin,2), 'LineStyle', 'none', ...
                                   'Marker', '.', ...
                                   'Color', blue, ...
                                   'MarkerFaceColor', blue, ...
                                   'MarkerSize', 4);
        h(1) = patch([0 0],[0 0], blue, 'EdgeColor','none');
        drawFootprint(good, blue, 0.3);
    end
    xlabel('z_{1}'); ylabel('z_{2}'); title([titlelabel ' Footprints']);
end
legend(h(h~=0), lbls(h~=0), 'Location', 'NorthEastOutside');
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square; axis([lbound(1)-1 ubound(1)+1 lbound(2)-1 ubound(2)+1]);

end
% =========================================================================
function handle = drawFootprint(footprint, color, alpha)
% 
hold on;
if isempty(footprint) || isempty(footprint.polygon)
    if size(Z,2) == 3
        handle = patch([0 0],[0 0],[0 0], color, 'EdgeColor','none');
    else
        handle = patch([0 0],[0 0], color, 'EdgeColor','none');
    end
    return
end

handle = plot(footprint.polygon,'FaceColor', color, 'EdgeColor','none', 'FaceAlpha', alpha);
hold off;

end
% =========================================================================
function h = drawBinaryPerformance(Z, Ybin, titlelabel)

ubound = ceil(max(Z));
lbound = floor(min(Z));
orange = [1.0 0.6471 0.0];
blue = [0.0 0.0 1.0];
lbls = {'GOOD','BAD'};
h = zeros(1,2);
if size(Z,2) == 3
    if any(~Ybin)
        h(2) = patch([0 0],[0 0],[0 0], orange, 'EdgeColor','none');
        line(Z(~Ybin,1), Z(~Ybin,2), Z(~Ybin,3), 'LineStyle', 'none', ...
                                     'Marker', '.', ...
                                     'Color', orange, ...
                                     'MarkerFaceColor', orange, ...
                                     'MarkerSize', 4);
    end
    if any(Ybin)
        h(1) = patch([0 0],[0 0],[0 0], blue, 'EdgeColor','none');
        line(Z(Ybin,1), Z(Ybin,2), Z(Ybin,3), 'LineStyle', 'none', ...
                                   'Marker', '.', ...
                                   'Color', blue, ...
                                   'MarkerFaceColor', blue, ...
                                   'MarkerSize', 4);
    end
    view([1,1,1])
    grid on
    xlabel('z_{1}'); ylabel('z_{2}'); zlabel('z_{3}'); title(titlelabel);
else
    if any(~Ybin)
        h(2) = patch([0 0],[0 0], orange, 'EdgeColor','none');
        line(Z(~Ybin,1), Z(~Ybin,2), 'LineStyle', 'none', ...
                                     'Marker', '.', ...
                                     'Color', orange, ...
                                     'MarkerFaceColor', orange, ...
                                     'MarkerSize', 4);
    end
    if any(Ybin)
        h(1) = patch([0 0],[0 0], blue, 'EdgeColor','none');
        line(Z(Ybin,1), Z(Ybin,2), 'LineStyle', 'none', ...
                                   'Marker', '.', ...
                                   'Color', blue, ...
                                   'MarkerFaceColor', blue, ...
                                   'MarkerSize', 4);
    end
    xlabel('z_{1}'); ylabel('z_{2}'); title(titlelabel);
end
legend(h(h~=0), lbls(h~=0), 'Location', 'NorthEastOutside');
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square; axis([lbound(1)-1 ubound(1)+1 lbound(2)-1 ubound(2)+1]);

end
% =========================================================================