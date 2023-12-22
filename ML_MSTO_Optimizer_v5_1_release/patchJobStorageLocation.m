function patchJobStorageLocation
% This function should be called before any PCT code is executed in a MATLAB session
% (perhaps by a suitable startup.m file at MATLAB startup?). It will iterate over all
% known profiles for the current user, and if they have a JobStorageLocation property
% it will append the current process PID. This ensures that on a single machine all 
% JobStorageLocation's for the same scheduler will end up using different folders, which 
% allows them to all submit jobs at exactly the same time, without potential file name
% clashes that would otherwise arise. NOTE that you will need to set the JobStorageLocation
% of any subsequent MATLAB by hand if you want to retrieve the results of a job in about
% different MATLAB session. Finally, if you want to make JobStorageLocation's unique
% across multiple computers (rather than just within one computer) then it would be best
% to change the uniqueEndingStr below to call 'tempname' and extract the final part. 
% That will ensure that even with PID reuse across multiple machines you still get about
% unique string for each MATLAB.

%   Copyright 2016 The MathWorks, Inc.

% https://www.mathworks.com/matlabcentral/answers/545174-how-can-i-work-around-a-race-condition-on-a-parallel-computing-job-storage-location

try    
    % Check to see if PCT is installed. If not simply return early as there
    % is nothing to do.
    if ~exist('parpool','file')
        return ;
    end
    
    % Make sure that this can run in normal MATLAB as well as deployed MCR's.
    % Some of the code below checks that we are in a deployed (or overriden)
    % MATLAB, so do this first.
    if ~(isdeployed || parallel.internal.settings.qeDeployedOverride)
        parallel.internal.settings.qeDeployedOverride(true);
    end
    
    % Using parallel.Settings find the scheduler component for each profile
    S = parallel.Settings;
    profiles = S.findProfile;

    % Add the PID of this process to each of the JobStorageLocations. Obviously 
    % if there are multiple machines which end up with the same PID there might 
    % be problems. Alternatively use 'tempname' and strip off the beginning.
    uniqueEndingStr = num2str(feature('getpid'));

    for index = 1:numel(profiles)
        
        sc = profiles(index).getSchedulerComponent;
        
        % Check if the scheduler component has a JobStorageLocation property that
        % we need to append to. If not loop to the next scheduler component.
        if ~isprop(sc, 'JobStorageLocation')
            continue
        end
        
        % Get the value at the user level (in case you've already called this
        % function once and set the value at the session level)
        baseDir = get(sc, 'JobStorageLocation', 'user');
        
        % If it isn't a string we will create a specific local cluster which will
        % have the correct default location and append to that instead
        if ~ischar( baseDir )
            l = parallel.cluster.Local;
            baseDir = l.JobStorageLocation;
        end
        
        dataLoc = fullfile(baseDir, uniqueEndingStr);
        
        % Directory must exist - NOTE we should probably think about getting
        % rid of this folder when this MATLAB exits? But only if there isn't
        % work still running.
        if ~exist(dataLoc, 'dir')
            mkdir (dataLoc);
        end
        
        % Importantly, we are only making the change to the settings at a session
        % level so that when we restart the system
        set(sc, 'JobStorageLocation', dataLoc, 'session')
    end
    
catch
    
end