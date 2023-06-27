function [ fva ] = FVA( FBAO, varargin)
% Performs flux variability analysis on all reactions of the given
% FBAObject. If the FBAObject has split reversible reactions (non-
% negative lower bounds), then a variable with name 'rev_split'
% specifies the forward directions of reversible reactions. This
% ensures that the reversed direction is set to zero when maximizing
% to avoid futile cycles with maximum flux equalling the upper bound.
% The returned array has size of the number of specified reactions:
% (to-from+1), if from and to are specified, or length(reactions),
% if reactions is specified.

    
	evalc('run initTomlab');
    
    tolerance = 0.0001;
    
    S = FBAO.S;
    n = FBAO.n;
    lowerBound = FBAO.B(1,:)';
    upperBound = FBAO.B(2,:)';
    
    % default arguments
    minimize = 1;
    maximize = 1;
    from = 1;
    to = n;
	rev_split = zeros(n, 1);
    
    arg = 1;
	while (arg+1 <= nargin)
        if (strcmp(inputname(arg+1), 'solver'))
            solver = varargin{arg};        
        elseif (isnumeric(varargin{arg}) || islogical(varargin{arg}))
            if (strcmp(inputname(arg+1), 'from'))
                from = varargin{arg};
            elseif (strcmp(inputname(arg+1), 'to'))
                to = varargin{arg};
			elseif (strcmp(inputname(arg+1), 'rev_split'))
				rev_split = varargin{arg};
			elseif (strcmp(inputname(arg+1), 'rev'))
				rev = varargin{arg};
            end
        elseif (strcmp(varargin{arg}, 'max'))
            minimize = 0;
        elseif (strcmp(varargin{arg}, 'min'))
            maximize = 0;
        elseif (iscell(varargin{arg}))
            reactions = varargin{arg};
        elseif (ischar(varargin{arg}))
            reactions = {varargin{arg}};            
        end
        arg = arg+1;
	end
    
    % determine the indices of the reactions 
    indices = from:to;
    if (exist('reactions', 'var'))
        indices = FBAO.nameToIndex(reactions);
    end
	rev_split = rev_split(indices);

    if (~exist('solver', 'var'))
        solver = 'cplex';
    end
    
    options = struct;
	options.name = 'FVA';
	options.solver = solver;    
    
    if (minimize && maximize)
        fva = zeros(2, to-from+1);
    else
        fva = zeros(1, to-from+1);
    end
            
    % decision variables
    v = tom('v', n, 1);
    
    % define the basic problem once, specify it later
    c = zeros(n, 1);
    objective = c'*v;
    constraints = { ...
        S*v == 0;
        lowerBound <= v <= upperBound;
    };

    Prob = sym2prob('lp', objective, constraints, [], options);
    Prob.optParam.IterPrint = 0; % Set to 1 to see iterations.
 	Prob.MIP.cpxControl.LPMETHOD = 1;
    Prob.MIP.cpxControl.SCAIND=1;
    Prob.Solver.Alg = 2; % Depth First, then Breadth search
	
	blocked = 0;
	blockedFluxes = 0;
    num = 0;
    charCount = 0;    
	for i = indices
        % print i as progress output
        num = num+1;
        if (charCount > 0)
            for k=1:charCount
                fprintf('\b');
            end
            charCount = 0;
        end
        if (minimize && maximize)
            message = sprintf('Minimizing/maximizing flux %d of %d.', num, length(indices));
        elseif (minimize)
            message = sprintf('Minimizing flux %d of %d.', num, length(indices));
        elseif(maximize)
            message = sprintf('Maximizing flux %d of %d.', num, length(indices));
        end
        fprintf(message);
        charCount = charCount+length(message);
        
        if (minimize)
            % set the objective: minimize v_i
            c = zeros(n, 1);
            c(i) = 1;
            Prob = modify_c(Prob, c);
            tomrun = tomRun(solver, Prob, 0);
            
            % check the solution
            if (tomrun.ExitFlag == 4)
				v_i_min = 0;
            elseif (tomrun.ExitFlag == 0)
				solution = getSolution(tomrun);
				v_i_min = solution.v(i);
				
				% fix numerical imprecisions
                if (v_i_min < lowerBound(i))
                    if (v_i_min < lowerBound(i) - tolerance)
						save(strcat('Error013_', datestr(now, 'yyyy-mm-dd_HHMM'), '.mat'));
						error('Maximum larger than upper bound for flux %d (reaction %d): %d, ub=%d.', i, rev_split(i), solution.v(i), upperBound(i));
                    end
                    v_i_min = lowerBound(i);
                end
            else
                save(strcat('Error012_', datestr(now, 'yyyy-mm-dd_HHMM'), '.mat'));
				error('Error while maximizing flux %d of reaction %d. ExitFlag=%d. %s.\n', i, unsplitIndex(i, rev), tomrun.ExitFlag, tomrun.ExitText);
            end
            
			% set the value in fva(1,i)
            fva(1, i) = v_i_min;
        end
        
        if (maximize)
            % set the objective: maximize v_i
            c = zeros(n, 1);
            c(i) = -1;
            Prob = modify_c(Prob, c);
			
			if (find(rev_split))
                % avoid futile cycles through both directions of a 
                % reversible reaction
                % reset lower and upper bounds
				x_L = lowerBound;
				x_U = upperBound;
%                 x_L(1) = 1;
                Prob = modify_x_L(Prob, x_L);
                Prob = modify_x_U(Prob, x_U);
                % set reversed of i to 0 to avoid reversible futile cycle
                if (rev_split(i))
					x_L(i+1) = 0;
					x_U(i+1) = 0;
                    Prob = modify_x_L(Prob, x_L);
                    Prob = modify_x_U(Prob, x_U);
                elseif (i > 1 && rev_split(i-1))
					x_L(i-1) = 0;
					x_U(i-1) = 0;
                    Prob = modify_x_L(Prob, x_L);
                    Prob = modify_x_U(Prob, x_U);
                end
			end
			
            tomrun = tomRun(solver, Prob, 0);
            
            % check the solution
            if (tomrun.ExitFlag == 4)
				v_i_max = 0;
			elseif (tomrun.ExitFlag == 0)
				solution = getSolution(tomrun);
				v_i_max = solution.v(i);
				
				% fix numerical imprecisions
				if (v_i_max < 0)
					% tolerate a small error and set to zero
					if (v_i_max > -tolerance)
						v_i_max = 0;
					else
						save(strcat('Error008_', datestr(now, 'yyyy-mm-dd_HHMM'), '.mat'));
						error('Negative maximum value for flux %d (reaction %d): %d.', i, rev_split(i), i, solution.v(i));
					end
				elseif (v_i_max > upperBound(i))
					if (v_i_max > upperBound(i) + tolerance)
						save(strcat('Error013_', datestr(now, 'yyyy-mm-dd_HHMM'), '.mat'));
						error('Maximum larger than upper bound for flux %d (reaction %d): %d, ub=%d.', i, rev_split(i), solution.v(i), upperBound(i));
					end
					v_i_max = upperBound(i);
				end
			else
                save(strcat('Error012_', datestr(now, 'yyyy-mm-dd_HHMM'), '.mat'));
				error('Error while maximizing flux %d of reaction %d. ExitFlag=%d. %s.\n', i, unsplitIndex(i, rev), tomrun.ExitFlag, tomrun.ExitText);
            end
			
			% count the blocked reactions
			if (v_i_max == 0)
				% if the reaction is irreversible or both directions are blocked count the
				% reactions as blocked
                irrev = (~rev_split(i) && (i == 1 || ~rev_split(i-1)));
				if (irrev || (~rev_split(i) && i > 1 && (minimize && fva(2, i-1) == 0 || ~minimize && fva(1, i-1) == 0)))
					blocked = blocked + 1;
                else
					blockedFluxes = blockedFluxes + 1;
				end
			end
			
			% set the value in fva(1,i) (max) or fva(2,i) (min/max)
            if (minimize)
				fva(2, i) = v_i_max;
				if (fva(1, i) > fva(2, i) + tolerance)
					save(strcat('Error011_', datestr(now, 'yyyy-mm-dd_HHMM'), '.mat'));
					error('\nmin v(%d)=%d, max v(%d)=%d.', i, fva(1, i), i, fva(2, i));
				end
            else
				fva(1, i) = v_i_max;
            end
        end
	end
    fprintf('\n%d blocked reactions, %d reversible reactions are irreversible.\n', blocked, blockedFluxes);
end

