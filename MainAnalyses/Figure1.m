%% Visual speed tuning class characterisation
% Figure 1
addpath 'D:\Code\GenericFunctions'

% --- TEMPORARY CODE FOR DEPENDENCY ANALYSIS ---
% This block only runs if the variables weren't loaded, as during static analysis.
if ~exist('goodUnits', 'var')
    goodUnits(1) = struct; % Creates an empty table placeholder
    goodUnits(2) = struct;
end
% ---------------------------------------------

%% get excitatory and suppressive responses for each tuning curve
% not calculated in processing script so done here instead

FRscale = 5;

for iunit = 1:numel(goodUnits)
    
    % stationary
    goodUnits(iunit).sigChangeP_stat = nan(4,7);
    goodUnits(iunit).sigChangeBool_stat = nan(4,7);
    goodUnits(iunit).ExcBool_stat = nan(4,7);
    goodUnits(iunit).nSig_stat = nan(4,1);
    goodUnits(iunit).nExc_stat = nan(4,1);
    goodUnits(iunit).nSupp_stat = nan(4,1);

    % find if stimulus evoked FR is higher or lower than baseline
    goodUnits(iunit).ExcBool_stat = double(cellfun(@mean, goodUnits(iunit).spikeCounts_stat)>(cellfun(@mean, goodUnits(iunit).baselineSpikeCounts_stat)*FRscale));
    goodUnits(iunit).ExcBool_stat(goodUnits(iunit).ExcBool_stat==0)=-1; % change zeros to minus 1

    % use sign-rank test to get sig difference p value for baseline vs
    % stimulus evoked
    for irow = 1:size(goodUnits(iunit).baselineSpikeCounts_stat,1)
        for icol = 1:size(goodUnits(iunit).baselineSpikeCounts_stat,2)
            if ~isempty(goodUnits(iunit).spikeCounts_stat{irow,icol})
                goodUnits(iunit).sigChangeP_stat(irow,icol) = ...
                    signrank(goodUnits(iunit).spikeCounts_stat{irow,icol},...
                    goodUnits(iunit).baselineSpikeCounts_stat{irow,icol}.*FRscale);
            end
        end
    end

    % convert to boolean significance test using holmbonferroni corrected p<0.05
    for irow = 1:4
        goodUnits(iunit).sigChangeP_stat(irow,:) = holmbonferroni(goodUnits(iunit).sigChangeP_stat(irow,:));
    end
    goodUnits(iunit).sigChangeBool_stat = goodUnits(iunit).sigChangeP_stat<0.05; % bonferroni correction
    goodUnits(iunit).nSig_stat = sum(goodUnits(iunit).sigChangeBool_stat,2)'; % number of significant responses
    goodUnits(iunit).nExc_stat = sum(goodUnits(iunit).sigChangeBool_stat==1 & goodUnits(iunit).ExcBool_stat==1,2)'; % number of sig excitatory responses
    goodUnits(iunit).nSupp_stat = sum(goodUnits(iunit).sigChangeBool_stat==1 & goodUnits(iunit).ExcBool_stat==-1,2)'; % number of sig suppressive responses

    % locomotion
    goodUnits(iunit).sigChangeP_run = nan(4,7);
    goodUnits(iunit).sigChangeBool_run = nan(4,7);
    goodUnits(iunit).ExcBool_run = nan(4,7);
    goodUnits(iunit).nSig_run = nan(4,1);
    goodUnits(iunit).nExc_run = nan(4,1);
    goodUnits(iunit).nSupp_run = nan(4,1);

    % find if stimulus evoked FR is higher or lower than baseline
    goodUnits(iunit).ExcBool_run = double(cellfun(@mean ...
        , goodUnits(iunit).spikeCounts_run)>(cellfun(@mean, goodUnits(iunit).baselineSpikeCounts_run)*FRscale));
    goodUnits(iunit).ExcBool_run(goodUnits(iunit).ExcBool_run==0)=-1; % change zeros to minus 1

    % use sign-rank test to get sig difference p value for baseline vs
    % stimulus evoked
    for irow = 1:size(goodUnits(iunit).baselineSpikeCounts_run,1)
        for icol = 1:size(goodUnits(iunit).baselineSpikeCounts_run,2)
            if ~isempty(goodUnits(iunit).spikeCounts_run{irow,icol})
                goodUnits(iunit).sigChangeP_run(irow,icol) = ...
                    signrank(goodUnits(iunit).spikeCounts_run{irow,icol},...
                    goodUnits(iunit).baselineSpikeCounts_run{irow,icol}.*FRscale);
            end
        end
    end

    % convert to boolean significance test using bonferroni corrected p<0.05
    for irow = 1:4
        goodUnits(iunit).sigChangeP_run(irow,:) = holmbonferroni(goodUnits(iunit).sigChangeP_run(irow,:)); 
    end
    goodUnits(iunit).sigChangeBool_run = goodUnits(iunit).sigChangeP_run<0.05; % bonferroni correction
    goodUnits(iunit).nSig_run = sum(goodUnits(iunit).sigChangeBool_run,2)'; % number of significant responses
    goodUnits(iunit).nExc_run = sum(goodUnits(iunit).sigChangeBool_run==1 & goodUnits(iunit).ExcBool_run==1,2)'; % number of sig excitatory responses
    goodUnits(iunit).nSupp_run = sum(goodUnits(iunit).sigChangeBool_run==1 & goodUnits(iunit).ExcBool_run==-1,2)'; % number of sig suppressive responses


end


%% get values of basic properties for each tuning class (gc)

allTuning = cat(1,goodUnits.tuning_stat,goodUnits.tuning_run);
allSSI = cat(1,goodUnits.SSI_stat,goodUnits.SSI_run);
allr2 = [[goodUnits.r2_stat], [goodUnits.r2_run]]';
allr2p = [[goodUnits.r2pval_stat], [goodUnits.r2pval_run]]';
allGC = [[goodUnits.gaussChar_stat], [goodUnits.gaussChar_run]]';
allNSig = [[goodUnits.nSig_stat], [goodUnits.nSig_run]]';
allNExc = [[goodUnits.nExc_stat], [goodUnits.nExc_run]]';
allNSupp = [[goodUnits.nSupp_stat], [goodUnits.nSupp_run]]';
allPrefSpeed = [[goodUnits.prefSpeed_stat],[goodUnits.prefSpeed_run]]';

r2p_thresh = 0.05;
r2_thresh=0.1;

for igc = 1:4
    idx = find(allr2>=r2_thresh & allr2p<r2p_thresh & allGC==igc); % tuned cells of this class
    gc(igc).r2vals = allr2(idx);
    gc(igc).nSigvals = allNSig(idx);
    gc(igc).nExcvals = allNExc(idx);
    gc(igc).nSuppvals = allNSupp(idx);
    gc(igc).prefSpeed = allPrefSpeed(idx);
    gc(igc).allTuning = allTuning(idx,:);
    gc(igc).allSSI = allSSI(idx,:);
    [~,gc(igc).prefSSI] = max(gc(igc).allSSI,[],2);
end


%% Example tuning/SSI curves and population means

gaussFun =  @(params,xdata) params(1) + params(2).*exp(-(((xdata-params(3)).^2)/(2*(params(4).^2))));
IDs = [951005818,951167265,950997203,951002872];
plotOrder = [1 3 4 2]; % lowpass, bandpass, highpass, bandreject
idx=[];
IDs = IDs(plotOrder);
for iid = 1:numel(IDs)
    idx(iid)=find(cat(1,allUnits.ID)==IDs(iid));
end
idirs = [2 4 2 4];

figure
for iunit = 1:numel(idx)
subplot(4,4,iunit), hold on
    % plot trial spike counts for each speed
    for ispeed = 1:7
        plot(repelem(ispeed,1,numel(allUnits(idx(iunit)).spikeCounts_run{idirs(iunit),ispeed})),...
            allUnits(idx(iunit)).spikeCounts_run{idirs(iunit),ispeed}, 'k.', 'MarkerSize',10)
    end
    % plot best gaussian fit
    plot(1:0.01:7, feval(gaussFun, allUnits(idx(iunit)).gaussParams_run(idirs(iunit),:),  1:0.01:7), 'r-')
    ax = gca; ax.XTick = 1:7;
    defaultAxesProperties(gca, true)
    % plot corresponding SSI
    subplot(4,4,iunit+4)
    plot(1:7, allUnits(idx(iunit)).SSI_run(idirs(iunit),:), 'k-')
    ylim([0.5 1.9]); 
    ax= gca; ax.YTick = 0.5:0.2:1.9; ax.XTick = 1:7;
    defaultAxesProperties(gca, true)
    
end

% average tuning curves
ii=0;
for igc = plotOrder
    ii=ii+1;
    subplot(4,4,8+ii)
    shadedErrorBar(1:7, nanmean(gc(igc).allTuning), nansem(gc(igc).allTuning,1)*1.96)
    ylim([6 18]), ax = gca; ax.XTick= 1:7; ax.YTick = 6:2:18;
    defaultAxesProperties(gca, true)
end
% average SSI curves
ii=0;
for igc = plotOrder
    ii=ii+1;
    subplot(4,4,12+ii)
    shadedErrorBar(1:7, nanmean(gc(igc).allSSI), nansem(gc(igc).allSSI,1)*1.96)
    ylim([0.7 0.9]), ax = gca; ax.XTick= 1:7; ax.YTick = 0.7:0.1:0.9;
    defaultAxesProperties(gca, true)
end

%% Distribution of tuning classes

gcNums = cellfun(@numel, {gc.r2vals});
allChar = gcNums./sum(gcNums);
nTuned_forEachClass = gcNums(plotOrder)

figure
bar(1:4, allChar(plotOrder))
ax=gca; ax.XTick = 1:4; ax.XTickLabel = {'lowpass', 'bandpass', 'bandreject', 'highpass'};
ylabel('Probability')
defaultAxesProperties(gca,false)

%% excitation and suppression by tuning class

for igc = 1:4
    SigVals_mean(igc) = mean(gc(igc).nSigvals);
    SigVals_sem(igc) = sem(gc(igc).nSigvals);

    ExcVals_mean(igc) = mean(gc(igc).nExcvals);
    ExcVals_sem(igc) = sem(gc(igc).nExcvals);

    SuppVals_mean(igc) = mean(gc(igc).nSuppvals);
    SuppVals_sem(igc) = sem(gc(igc).nSuppvals);

    countTable(igc,:) = [sum(gc(igc).nExcvals), sum(gc(igc).nSuppvals)];
end

figure, hold on
bar(1:4, [ExcVals_mean(plotOrder); SuppVals_mean(plotOrder)],'stacked')
ylabel('# Significant responses')
legend({'excited', 'suppressed'})

% X2 test of independence
[X2stat, p, df] = chi2test_countTable(countTable)

%% p(SSI|pref speed) conditional distribution

allPrefSSI = cat(1,gc.prefSSI);
allPrefSpeeds = cat(1,gc.prefSpeed);

vals = histcounts2(allPrefSSI(:), allPrefSpeeds(:));
vals_norm = vals./sum(vals,1);
figure
imagesc(vals_norm), colorbar, colormap(cmocean('tempo')), axis xy
xlabel('Pref speed'), ylabel('Pref SSI'), caxis([0, 0.5])


%% distribution of preferred speeds

allPref = cat(1,gc.prefSpeed);
allChar = repelem(1:4, cellfun(@numel, {gc.prefSpeed}));

ii=0;
for ichar = plotOrder
    ii=ii+1;
    gcPrefs = allPref(allChar==ichar);
    for ispeed = 1:7
        gcSpeedVals(ispeed,ii)=sum(gcPrefs==ispeed);
    end
end

figure
bar(gcSpeedVals./sum(gcSpeedVals(:)),'stacked')
defaultAxesProperties(gca,false)
title('Preferred speeds')
ax=gca; ax.XTick = 1:7; ax.XTickLabel = [0 16 32 64 128 256 512]; 
xlabel('Speed (deg/s)'), ylabel('Probability')