%% Supp Figure 1: Locomotion speed distributions for both datasets

load('behaviourSampling_dg.mat') % drifting gratings dataset

idg=0;
for isesh = 1:numel(sesh)
    % if isesh==13 | isesh==26
        % continue
    % end
    for itf = 1:5
        theseVals = sesh(isesh).allRun(:,:,itf);

        if any(~isnan(theseVals))
            idg=idg+1;
            isesh
            itf

            [f2,x2] = ecdf(theseVals(:));
            figure(990), hold on
            plot(x2,f2,'c')
           
        end
    end    
end


clear sesh
load('behaviourSampling_dm.mat') % dot motion dataset
idm=0;
for isesh = 1:numel(sesh)
    if any(~isnan(sesh(isesh).allRun(:)))
    for idir = 1:4

        theseVals = sesh(isesh).allRun(:,:,idir);
        if any(~isnan(theseVals))
        idm=idm+1;
            [f2,x2] = ecdf(theseVals(:));
            figure(990), hold on
            plot(x2,f2,'m')

        end
    end
    end
end
