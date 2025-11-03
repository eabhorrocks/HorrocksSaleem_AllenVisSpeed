%% Supp Figure 2: Hierarchical sorting of tuning curves

%% Get tuning curves
load('goodUnits.mat')

% --- Thresholding and Data Preparation ---
r2_thresh  = 0.1;
r2p_thresh = 0.05;

allTuning = cat(1, goodUnits.tuning_stat, goodUnits.tuning_run);
allr2     = [[goodUnits.r2_stat], [goodUnits.r2_run]]';
allr2p    = [[goodUnits.r2pval_stat], [goodUnits.r2pval_run]]';

% Select units that pass criteria:
idx         = find(allr2 >= r2_thresh & allr2p < r2p_thresh);
allCurves   = allTuning(idx, :);
allCurves_z = zscore(allCurves, [], 2);   % z-score each curve (row-wise)
nCurves     = size(allCurves_z, 1);

r2tunedOnly = allr2(idx);

%% Sorting using optimalleaforder
% Compute Dissimilarity Matrix & Hierarchical Linkage
D  = pdist(allCurves_z, 'euclidean');
dm = squareform(D);
Z  = linkage(D, 'average');

% find the optimal leaf order
leafOrder = optimalleaforder(Z, D);

% re-order the dissimilarity matrix and tuning curves
dm_reordered = dm(leafOrder, leafOrder);
tuning_reordered = allCurves_z(leafOrder, :);
r2_reordered = r2tunedOnly(leafOrder,:);
tuning_orig_reordered = allCurves(leafOrder, :);

%% plots

figure
subplot(121)
imagesc(allCurves_z), caxis([-2.25 2.25]), axis xy
subplot(122)
imagesc(allCurves_z(leafOrder,:)),caxis([-2.25 2.25]), axis xy

figure
imagesc(dm_reordered)
caxis([0 4.9])
axis xy



%% plot examples

i2plot = [30 818 2075 3837 4450];

for i = i2plot
    figure
    plot(tuning_orig_reordered(i,:));
    xlim([0.5 7.5])
    ax=gca; ax.YLim(1) = ax.YLim(1)-0.5; ax.YLim(2) = ax.YLim(2)+0.5;
    
    defaultAxesProperties(gca, false)
        box on
        strplot = sprintf('exampleTuning_%s',num2str(i));
        % print(gcf,strplot,'-dsvg', '-painters')
        title(i)
end
