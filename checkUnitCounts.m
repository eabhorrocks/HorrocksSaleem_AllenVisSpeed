%% check unit counts


areas = {'VISp', 'VISl', 'VISal', 'VISrl', 'VISam', 'VISpm', 'LGd', 'LP'};


for iarea = 1:8
    areaUnits = goodUnits(strcmp([goodUnits.ecephys_structure_acronym], areas(iarea)));

    allr2_stat = cat(1,areaUnits.r2_stat);
    allr2_run = cat(1,areaUnits.r2_run);

    nstat(iarea) = sum(any(~isnan(allr2_stat),2));
    nrun(iarea) = sum(any(~isnan(allr2_run),2));
end