function runAll( network, varargin )
% Run flux coupling analysis, create flux coupling graph, generate flux
% target states, and calculate the drivers, on a single metabolic network.
% 'network' specifies the name of the network files, e.g. 'ControlNetwork'.
% 'solver' specifies the optimization solver for Tomlab (optional). The
% default is 'cplex'.

% print date and process ID
fprintf(strcat(datestr(now), '. PID: %d.\n'), sscanf(evalc('feature getpid'), '%*s%*s%d'));


arg = 1;
while (arg+1 <= nargin)
    if (ischar(varargin{arg}))
        solver = varargin{arg};
    elseif (isnumeric(varargin{arg}))
        series = varargin{arg};
    end
    arg = arg+1;
end

if (~exist('solver', 'var'))
    solver = 'cplex';
end

t1 = tic;

if (~exist('series', 'var'))
    % run all calculations on a single metabolic network
    % blocked reactions, full, partial, and directional coupling
    [S, rev, Reactions, fctable, blocked] = FCA(network);

    % anti- and inhibitor coupling
    % TODO comment this line and uncomment the next to run all analysis
    coupling = Coupling(network, S, rev, fctable, blocked, solver, 'inhibitor');
%     coupling = Coupling(network, S, rev, fctable, blocked, solver, 'essential', 'anti', 'inhibitor');

    % create profiles and FC graph
    createFCGraph(network, rev, coupling, blocked, Reactions);

    % sample random reaction activity patterns
    w = randomFluxConfiguration(network, 1000, 1, S, rev, blocked, Reactions, solver);

    % determine the drivers ('full': only fully specified activity patterns)
    sampleDrivers(network, w, rev, blocked, coupling, solver, 'full');
    
else
    % run all calculations on a set of networks
    anticoupled = 0;
	errors = zeros(length(series), 1);
	
    for i=1:length(series)
		n = series(i);
        if (isnumeric(n))
            name = strcat(network, '.', num2str(n));
        else
            name = strcat(network, '.', n);
        end

    %	if (i > 0)
            % blocked reactions, full, partial, and directional coupling
            [S, rev, Reactions, fctable, blocked] = FCA(name);

        if (i > 9 && anticoupled == 0)
            % skip anti-coupling as we don't expect to see any more
            if (i == 10)
                fprintf('Skipping anti-coupling.\n');
            end
            coupling = Coupling(name, S, rev, fctable, blocked, 'inhibitor');
        else
            % anti- and inhibitor coupling
            coupling = Coupling(name, S, rev, fctable, blocked, 'essential', 'anti', 'inhibitor');
            anticoupled = anticoupled + length(find(coupling==5));
        end

        % check for import/export reactions error
        if (size(coupling, 1) == 1 && size(coupling, 2) == 1 && coupling(1,1) == -1)
            errors(i) = n;
            continue;
        end

        % create profiles and FC graph
        createFCGraph(name, rev, coupling, blocked, Reactions);

        % sample random flux configurations
        w = randomFluxConfiguration(name, 1000, 1, S, rev, blocked, Reactions);

        % determine the drivers
        sampleDrivers(name, w, rev, blocked, coupling, 'full');
        
        % free unused memory
        clear S; clear rev; clear Reactions; clear fctable; clear blocked;
        clear coupling; clear w;
    end
    
    if (length(find(errors)) > 0)
        fprintf('Errors in');
        fprintf(' %d', errors(find(errors)));
        fprintf('.\n');
    end
end

% print elapsed time
t1 = toc(t1);
d = floor(t1/86400);
h = floor(t1/3600);
minutes = floor(t1/60);
s = rem(t1, 60);
if (d > 0)
    fprintf('Total time for %s: %d days, %d hours.\n', network, d, h-(d*24));
elseif (h > 0)
    fprintf('Total time for %s: %d hours, %d minutes.\n', network, h, minutes-(h*60));
else
    fprintf('Total time for %s: %d minutes, %g seconds.\n', network, minutes, s);
end

end
