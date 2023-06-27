function [drivers] = sampleDrivers ( name, varargin )
% Calculate the drivers and plot the distribution of number of drivers
% the size of activity patterns.
% Creates the following files:
% (1) <name>.drivers: m x n matrix for m samples and n reactions, where
% (i,j) is 1 if reaction j is a driver in sample 1, 0 otherwise.
% (2) <name>.numdrivers: 1 x n vector for n reactions, where value j
% contains the average number of drivers for controlling j reactions.
% (3) <name>.alldrivers: m x n matrix for m samples and n reactions,
% where (i,j) contains the number of times reaction j is a driver for
% all random subsets of sample i.

    evalc('run initTomlab');
    
    % re-seed the random number generator
    rng('shuffle');

    mode = 1;

    % parse the arguments
    full = 0;
    arg = 1;
    while (arg+1 <= nargin)
        if (strcmp(inputname(arg+1), 'w'))
            w = varargin{arg};
        elseif (strcmp(inputname(arg+1), 'rev'))
            rev = varargin{arg};
        elseif (strcmp(inputname(arg+1), 'blocked'))
            blocked = varargin{arg};
        elseif (strcmp(inputname(arg+1), 'coupling'))
            coupling = varargin{arg};
        elseif (strcmp(inputname(arg+1), 'solver'))
            solver = varargin{arg};            
        elseif (ischar(varargin{arg}))
            if (strcmp(varargin{arg}, 'full'))
                full = 1;
            end
        elseif (isnumeric(varargin{arg}))
            mode = varargin{arg};
        end
        arg = arg+1;
    end
    
    charCount = 0;

    % load the data not passed as arguments
    if (~exist('blocked', 'var'))
        blocked = load(strcat(name, '.blocked'));
    end
    if (~exist('rev', 'var'))
        rev = load(strcat(name, '.rev'));
        rev = rev(~blocked);
    end
    if (~exist('coupling', 'var'))
    	coupling = load(strcat(name, '.coupling'));
    end
    if (~exist('w', 'var'))
        if (mode == 1)
            w = load(strcat(name, '.samples'));
        else
            w = load(strcat(name, '.samples', num2str(mode)));
        end
    end
    if (~exist('solver', 'var'))
        solver = 'cplex';
    end

    m = size(w, 1);
    n = size(w, 2);
    drivers = zeros(m, n);
    % calculate the subset drivers and all-drivers only for 100 samples
    numDrivers = zeros(min(m, 100), n);
    allDrivers = zeros(min(m, 100), n);
	% for allDrivers 101:1000
% 	allDrivers = zeros(900, n);
    driversTime = 0;
    
	% for allDrivers 101:1000 set to 101:m
    for sampleIndex=1:m
        % iterate over each sample
        if (charCount > 0)
            for k=1:charCount
                fprintf('\b');
            end
            charCount = 0;
        end
        message = sprintf('Calculating drivers for %s (%d reactions): sample %d of %d.', name, n, sampleIndex, m);
        fprintf(message);
        charCount = charCount + length(message);

        % uncomment for non-repetitive sampling
%         x = NaN(1, n);
%         perm = randperm(n);

        from = 1;
        % calculate the subset drivers only for 100 samples
		% for allDrivers 101:1000 comment this
        if (full || sampleIndex > 100)
            % only calculate the drivers for full flux configuratoins
            from = n;
        end
        for sampleSize=from:n
            % choose a new reaction randomly from sample m
%             indices = perm(1:sampleSize);
            % choose sampleSize reactions randomly from sample m
            specified = randperm(n, sampleSize);
            % determine the drivers for the specified indices in sample w
            [d, time] = Drivers(name, w(sampleIndex, :), specified, rev, blocked, coupling, solver);
            driversTime = driversTime + time;
			% for allDrivers 101:1000 uncomment this
% 			allDrivers(sampleIndex-100, :) = allDrivers(sampleIndex-100, :) + d';
            % store the number of drivers for each sample size
            if (~full && sampleIndex <= 100)
                numDrivers(sampleIndex, sampleSize) = length(find(d));
                if (numDrivers(sampleIndex, sampleSize) > sampleSize)
                    error('More drivers (%d) than size of the flux configuration (%d)!', numDrivers(sampleIndex, sampleSize), sampleSize);
                end
                % add the driver counts to the sample row for allDrivers
                allDrivers(sampleIndex, :) = allDrivers(sampleIndex, :) + d';
            end
            % save the drivers of the full sample
            if (sampleSize == n)
                % matlab automatically transposes d in this assignment
                drivers(sampleIndex, :) = d;
            end
        end
    end

    if (~full)
        % take the min, mean, and max number of drivers for each sample size
        % this mean should be more precise than mean(mean(drivers, 2))
        % because of rounding errors
        minima = min(numDrivers, [], 1);
        means = mean(numDrivers, 1);
        maxima = max(numDrivers, [], 1);
        stdevs = std(numDrivers, [], 1);
            
        % save the mean numbers of drivers and drivers of the full samples
        if (mode == 1)
			% for allDrivers 101:1000 comment this
            dlmwrite(strcat(name, '.mindrivers'), minima, '\t');
            dlmwrite(strcat(name, '.meandrivers'), means, '\t');
            dlmwrite(strcat(name, '.maxdrivers'), maxima, '\t');
            dlmwrite(strcat(name, '.stddrivers'), stdevs, '\t');
			% for allDrivers 101:1000 change output file name
            dlmwrite(strcat(name, '.alldrivers'), allDrivers, '\t');
        else
            dlmwrite(strcat(name, '.mindrivers.', num2str(mode)), minima, '\t');
            dlmwrite(strcat(name, '.meandrivers.', num2str(mode)), means, '\t');
            dlmwrite(strcat(name, '.maxdrivers.', num2str(mode)), maxima, '\t');
            dlmwrite(strcat(name, '.stddrivers.', num2str(mode)), stdevs, '\t');
            dlmwrite(strcat(name, '.alldrivers.', num2str(mode)), allDrivers, '\t');
        end
    end
    % write the drivers of the full configurations
	% for allDrivers 101:1000 comment this
    if (mode == 1)
        dlmwrite(strcat(name, '.drivers'), drivers, '\t');
    else
        dlmwrite(strcat(name, '.drivers.', num2str(mode)), drivers, '\t');
    end

    % print elapsed time
    t = driversTime;
    d = floor(t/86400);
    h = floor(t/3600);
    minutes = floor(t/60);
    s = rem(t, 60);
    if (d > 0)
        fprintf(' Time: %g days, %g hours.\n', d, h-(d*24));
    elseif (h > 0)
        fprintf(' Time: %g hours, %g minutes.\n', h, minutes-(h*60));
    else
        fprintf(' Time: %g minutes, %g seconds.\n', minutes, s);
    end

end
