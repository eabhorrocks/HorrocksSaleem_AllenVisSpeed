function trials = getRunAndPupilInfo(trials, runInfo, pupilInfo, preTime, postTime)

for itrial = 1:numel(trials)
    startrunidx=[];
    stoprunidx=[];
    startpupilidx=[];
    stoppupilidx=[];
    
    % get mean run speed
    %[~, startrunidx] = min(abs(runInfo.runMidTimes-trials(itrial).start_time));
    %[~, stoprunidx] = min(abs(runInfo.runMidTimes-trials(itrial).stop_time));
    
    startrunidx = find(runInfo.runTimeInterp>trials(itrial).start_time-preTime,1,'first');
    stoprunidx = find(runInfo.runTimeInterp>trials(itrial).stop_time+postTime,1,'first');
    
    trials(itrial).runTrace = runInfo.smthRunSpeedInterp(startrunidx:stoprunidx);
    trials(itrial).runTime = runInfo.runTimeInterp(startrunidx:stoprunidx) - trials(itrial).start_time;
    trials(itrial).meanRunSpeed = mean(trials(itrial).runTrace);
    trials(itrial).minRunSpeed = min(trials(itrial).runTrace);
    trials(itrial).maxRunSpeed = max(trials(itrial).runTrace);
    
    if isstruct(pupilInfo) % some sessions don't have eye tracking
        
        pupilInfo.pupilArea = pi.*pupilInfo.pupilWidth.*pupilInfo.pupilHeight;

        [~, startpupilidx] = min(abs(pupilInfo.pupilTime-trials(itrial).start_time));
        [~, stoppupilidx] = min(abs(pupilInfo.pupilTime-trials(itrial).stop_time));
        trials(itrial).pupilTime = pupilInfo.pupilTime(startpupilidx:stoppupilidx);

        trials(itrial).pupilAreaTrace = pupilInfo.pupilArea(startpupilidx:stoppupilidx);
        trials(itrial).meanPupilArea = nanmean(trials(itrial).pupilAreaTrace);
        
        trials(itrial).pupilTrace = [pupilInfo.pupil_x(startpupilidx:stoppupilidx);...
                                        pupilInfo.pupil_y(startpupilidx:stoppupilidx)];
    else
        trials(itrial).pupilAreaTrace = NaN;
        trials(itrial).meanPupilArea = NaN;
        trials(itrial).pupilTrace = NaN;
    end
end

end