function [ w ] = randomFluxConfiguration( name, m, mode, varargin )
% Generates m random binary vectors w(i,:), i=1,...m,
% such that Sv=0 for some v with w(i,:)=0 <=> v=0.
% w is a m x n matrix of binary vectors.
% Mode 1: deactivate a random set of imports, maximize a random set of
% export.
% Mode 3: start with an arbitrary flux configuration, generate the next by
% maximizing the manhattan distance to the average of all previous ones.
debug = false;

    % threshold for exchange fluxes
    epsilon = 1; % mode 1 is faster with epsilon = 0;
	if (mode == 1)
    % threshold for flux configurations: >0.001 doesn't maximize distance,
    % <1e-7 is slow for mode 3
%     epsilon2 = 1e-6; % allows violations of full coupling        
        epsilon2 = 1e-7; % works: 0.3, 0.1, 1e-2, 1e-4, 1e-6 (but small values allow violation of Sv=0)
	elseif (mode == 3)
        epsilon2 = 1e-7;
        if (strcmp(name, 'A_niger'))
            % A_niger doesn't work with 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-12, 1e-20
            epsilon2 = 1e-3;
        elseif (strcmp(name, 'H_sapiens_recon1') || strcmp(name, 'iUrothelialCancer1647'))
            % H_sapiens_recon1 and iUrothelialCancer1647 don't work with 1e-7
            epsilon2 = 1e-6;
        end
	end
    if (debug)
        fprintf('epsilon = %d. epsilon2 = %d.\n', epsilon, epsilon2);
    end
    
    evalc('run initTomlab');
    
    % re-seed the random number generator
    % this only affects different matlab sessions
    rng('shuffle');
    
	% parse the arguments
	arg = 1;
    while (arg+3 <= nargin)
        if (strcmp(inputname(arg+3), 'S'))
            S = varargin{arg};        
        elseif (strcmp(inputname(arg+3), 'rev'))
            rev = varargin{arg};
        elseif (strcmp(inputname(arg+3), 'blocked'))
            blocked = varargin{arg};
        elseif (strcmp(inputname(arg+3), 'Reactions'))
            Reactions = varargin{arg};
        elseif (strcmp(inputname(arg+3), 'coupling'))
            coupling = varargin{arg};
        elseif (strcmp(inputname(arg+2), 'solver'))
            solver = varargin{arg};            
        end
        arg = arg+1;
    end
    
    % load the data not provided as arguments
    if (~exist('blocked', 'var'))    
        blocked = load(strcat(name, '.blocked'));
    end
    if (~exist('S', 'var'))
        S = load(strcat(name, '.S'));
        S = S(:,~blocked);
    end
    if (~exist('rev', 'var'))
        rev = load(strcat(name, '.rev'));
        rev = rev(~blocked);
    end
    if (~exist('Reactions', 'var'))
%         Reactions = importdata(strcat(name, '.Reactions'));
        fid = fopen(strcat(name, '.Reactions'));
        Reactions = textscan(fid, '%s');
        fclose(fid);
        Reactions = Reactions{:};
        Reactions = Reactions(~blocked);
    end
    if (~exist('coupling', 'var'))
        coupling = load(strcat(name, '.coupling'));
    end    
    if (~exist('solver', 'var'))
        solver = 'cplex';
    end
    
    options = struct;
	options.name = name;
	options.solver = solver;
    n = size(S, 2);
    
    [coupledI,coupledJ,coupledV] = find(coupling);
    
    % set the default bounds
    lowerBound = zeros(n, 1);
    lowerBound(find(rev)) = -1000;
    upperBound = repmat(1000, n, 1);
    maxBound = max([abs(lowerBound); upperBound]);
    
    % initialize flux configuration to 0
    w = zeros(m, n);
        
    import = sum(abs(S), 1) == sum(S, 1);
    export = sum(abs(S), 1) == -sum(S, 1);
    numImport = length(find(import));
    numReversibleImport = length(find(import.*rev'));
    numExport = length(find(export));
    numReversibleExport = length(find(export.*rev'));
    exchange = import | export;
    numExchange = length(find(exchange));
    if (length(find(import)) + numReversibleExport < 2)
        fprintf(2, 'ERROR: %s has %d unblocked import/reversible export reaction (at least 2 required for sampling).\n', name, length(find(import)) + numReversibleExport);
        w = -1;
        return
    end
    if (length(find(export)) + numReversibleImport < 1)
        fprintf(2, 'ERROR: %s has %d unblocked export/reversible import reactions.\n', name, length(find(export)) + numReversibleImport);
        w = -1;
        return
    end
    
    tic
	
    if (mode == 1)
        fprintf('Sampling reaction activity patterns with %d import and %d export fluxes.\n', numImport+numReversibleExport, numExport+numReversibleImport);

        reject = 0;
        charCount = 0;
        i=1;
        while (i <= m)
            temp = -ones(1, n);
            % select random import reactions to block
            blockedImport = import .* (randi(2, 1, n)-1);
            blockedExport = export.*rev' .* (randi(2, 1, n)-1);
            % select random export reactions to maximize and their coefficients
            objImport = import.*rev' .* (rand(1, n));
            objExport = export .* (rand(1, n));
            % at least one import reaction must be left
            if (length(find(blockedImport)) == numImport && length(find(blockedExport.*rev')) == numReversibleExport)
                continue;
            end
            % at least one export reaction must be chosen
            if (length(find(objImport)) == 0 && length(find(objExport)) == 0)
                continue;
            end
            
            % reset lower and upper bounds
            lowerBound = zeros(n, 1);
            lowerBound(find(rev)) = -1000;
            upperBound = repmat(1000, n, 1);
            % set the upper bound of chosen import reactions to zero
            upperBound(find(blockedImport)) = 0;
            % set the lower bound of chosen reversible export reactions to zero
            lowerBound(find(blockedExport)) = 0;
            maxBound = max([abs(lowerBound); upperBound]);

            % determine a v with exchange reactions active as in x
            v = tom('v', n, 1);
            x = tom('x', numExchange, 1, 'integer');
            v_forward_ex = tom('v_forward_ex', numExchange, 1);
            v_reversed_ex = tom('v_reversed_ex', numExchange, 1);
            x_forward_ex = tom('x_forward_ex', numExchange, 1, 'integer');
            x_reversed_ex = tom('x_reversed_ex', numExchange, 1, 'integer');

            % maximize export through the chosen objective reactions
            obj = (objImport-objExport);
            objective = obj*v;
            constraints = { ...
                S*v == 0;
                lowerBound <= v <= upperBound;
                v(exchange) == v_forward_ex - v_reversed_ex;
                0 <= v_forward_ex <= x_forward_ex*maxBound;
                0 <= v_reversed_ex <= x_reversed_ex*maxBound;
                x_forward_ex + x_reversed_ex <= 1;
                epsilon*x <= v_forward_ex + v_reversed_ex <= maxBound*x;
                sum(x) >= 1;
                0 <= x_forward_ex <= 1;
                0 <= x_reversed_ex <= 1;
            };
            Prob = sym2prob('mip', objective, constraints, [], options);
            Prob.MIP.cpxControl.LPMETHOD = 1;
            Prob.MIP.cpxControl.SCAIND = 0;
            tomrun = tomRun(solver, Prob, 0);
            valid = 0;
            if (tomrun.ExitFlag == 0)
                solution = getSolution(tomrun);
                % check if the blocked exchange reactions are really zero
                if (length(find(solution.v(find(blockedImport)) > epsilon2)) ~= 0 || length(find(solution.v(find(blockedExport)) < -epsilon2)) ~= 0)
                    % possible numerical error. try different solver algorithm and scaling.
                    Prob.MIP.cpxControl.LPMETHOD = 0;
                    Prob.MIP.cpxControl.SCAIND = -1;
                    tomrun = tomRun(solver, Prob, 0);
                    if (tomrun.ExitFlag == 0)
                        solution = getSolution(tomrun);
                        if (length(find(solution.v(find(blockedImport)) > epsilon2)) == 0 && length(find(solution.v(find(blockedExport)) < -epsilon2)) == 0)
                            valid = 1;
                        else
                            % possible numerical error. try different solver algorithm and scaling.
                            Prob.MIP.cpxControl.LPMETHOD = 0;
                            Prob.MIP.cpxControl.SCAIND = 1;
                            tomrun = tomRun(solver, Prob, 0);
                            if (tomrun.ExitFlag == 0)
                                solution = getSolution(tomrun);
                                if (length(find(solution.v(find(blockedImport)) > epsilon2)) == 0 && length(find(solution.v(find(blockedExport)) < -epsilon2)) == 0)
                                    valid = 1;
                                end
                            end
                        end
                    end
                    Prob.MIP.cpxControl.LPMETHOD = 1;
                    Prob.MIP.cpxControl.SCAIND = 0;
                else
                    valid = 1;
                end
            end
            
            if (valid)
                temp = double(abs(solution.v) >= epsilon2);
                % check if there is a violation of full, partial,
                % directional, or anti-coupling
                for index=1:length(coupledI)
                    if ((coupledV(index) == 1 || coupledV(index) == 2) && (temp(coupledI(index)) ~= temp(coupledJ(index))))
                        valid = 0;
                        break;
                    elseif (coupledV(index) == 3 && temp(coupledI(index)) == 1 && temp(coupledJ(index)) == 0)
                        valid = 0;
                        break;
                    elseif (coupledV(index) == 5 && temp(coupledI(index)) == 0 && temp(coupledJ(index)) == 0)
                        valid = 0;
                        break;
                    end
                end
            end
            
            if (valid)
                % additional validation: check if w(i,:) has a valid v with x=0 <=> v=0 and Sv=0
                v = tom('v', n, 1);
                x = tom('x', numExchange, 1, 'integer');
                v_forward = tom('v_forward', n, 1);
                v_reversed = tom('v_reversed', n, 1);
                x_forward = tom('x_forward', n, 1, 'integer');
                x_reversed = tom('x_reversed', n, 1, 'integer');
                constraints = {...
                    S*v == 0;
                    lowerBound <= v <= upperBound;
                    v == v_forward - v_reversed;
                    epsilon2*x_forward <= v_forward <= x_forward*maxBound;
                    epsilon2*x_reversed <= v_reversed <= x_reversed*maxBound;
                    x_forward + x_reversed <= 1;
                    epsilon*x <= v_forward(exchange) + v_reversed(exchange) <= maxBound*x;
                    sum(x) >= 1;
                    temp == x_forward + x_reversed;
                    0 <= x_forward <= 1;
                    0 <= x_reversed <= 1;
                    0 <= v_forward;
                    0 <= v_reversed;
                    0 <= x <= 1;
                };
                Prob = sym2prob('mip', [], constraints, [], options);
                Prob.MIP.cpxControl.EPRHS = 1e-9; % min: 1e-9, max: 0.1, default: 1e-6
                tomrun = tomRun(solver, Prob, 0);
                if (tomrun.ExitFlag ~= 0)
                    % try different solver algorithm and scaling
                    Prob.MIP.cpxControl.LPMETHOD = 1;
                    Prob.MIP.cpxControl.SCAIND = -1;
                    tomrun = tomRun(solver, Prob, 0);
                    if (tomrun.ExitFlag ~= 0)
						Prob.MIP.cpxControl.LPMETHOD = 1;
						Prob.MIP.cpxControl.SCAIND = 1;
                        tomrun = tomRun(solver, Prob, 0);
                        if (tomrun.ExitFlag ~= 0)
                            Prob.MIP.cpxControl.LPMETHOD = 0;
                            Prob.MIP.cpxControl.SCAIND = -1;
                            tomrun = tomRun(solver, Prob, 0);
                            if (tomrun.ExitFlag ~= 0)
								Prob.MIP.cpxControl.LPMETHOD = 0;
								Prob.MIP.cpxControl.SCAIND = 0;
                                tomrun = tomRun(solver, Prob, 0);
								if (tomrun.ExitFlag ~= 0)
									Prob.MIP.cpxControl.LPMETHOD = 0;
									Prob.MIP.cpxControl.SCAIND = 1;
									tomrun = tomRun(solver, Prob, 0);
									if (tomrun.ExitFlag ~= 0)
										% give up
										valid = 0;
									end
								end
                            end
                        end
                    end
                    Prob.MIP.cpxControl.LPMETHOD = 1;
                    Prob.MIP.cpxControl.SCAIND = 0;
                end
            end

            if (charCount > 0)
                for k=1:charCount
                    fprintf('\b');
                end
                charCount = 0;
            end        
            if (valid)
                % temp is valid
                message = sprintf('Rejected: %d. Found: %d.', reject, i);
                fprintf(message);
                charCount = charCount + length(message);
                w(i,:) = temp;
                i = i + 1;
            else
                % repeat iteration if no solution
                reject = reject + 1;
                message = sprintf('Rejected: %d. Found: %d.', reject, i-1);
                fprintf(message);
                charCount = charCount + length(message);
            end
        end
    elseif (mode == 3)
		charCount = 0;
        for i=1:m
            % maximize the distance between the new w and the previous average
%             v = tom('v', n, 1);
%             x = tom('x', numExchange, 1, 'integer');
%             new = tom('new', n, 1, 'integer');
%             v_forward = tom('v_forward', n, 1);
%             v_reversed = tom('v_reversed', n, 1);        
%             x_forward = tom('x_forward', n, 1, 'integer');
%             x_reversed = tom('x_reversed', n, 1, 'integer');
			v = tom('v', n, 1);
			x = tom('x', numExchange, 1, 'integer');
            v_forward = tom('v_forward', n, 1);
            v_reversed = tom('v_reversed', n, 1);        
			x_forward = tom('x_forward', n, 1, 'integer');
			x_reversed = tom('x_reversed', n, 1, 'integer');
            
            if (i == 1)
                % start with a random average
                avg = rand(1, n, 1);
            else
                % generate the new average from all previous w
                avg = mean(w(1:(i-1), :), 1);
            end
            % transform the values to [-1, 1]
            avg = avg*2-1;
                v = tom('v', n, 1);
                x = tom('x', numExchange, 1, 'integer');
                v_forward_ex = tom('v_forward', n, 1);
                v_reversed_ex = tom('v_reversed', n, 1);        
                x_forward_ex = tom('x_forward', n, 1, 'integer');
                x_reversed_ex = tom('x_reversed', n, 1, 'integer');
                % maximize the difference between avg and new in [-1, 1]:
                % original objective maximizing distance to avg
                objective = avg*((x_forward + x_reversed)*2-1);
                % epsilon2 disallows values 0 < abs(v) < epsilon2
                constraints = { ...
                    S*v == 0;
                    lowerBound <= v <= upperBound;
                    v == v_forward - v_reversed;
                    0 <= v_forward <= x_forward*maxBound;
                    0 <= v_reversed <= x_reversed*maxBound;
                    x_forward + x_reversed <= 1;
                    epsilon*x <= v_forward(exchange) + v_reversed(exchange) <= maxBound*x;
                    sum(x) >= 1;
                    0 <= x_forward <= 1;
                    0 <= x_reversed <= 1;
                };
            Prob = sym2prob('mip', objective, constraints, [], options);
			Prob.MIP.cpxControl.LPMETHOD = 1;
			% mode 3 runs more stable without scaling
			Prob.MIP.cpxControl.SCAIND = -1;
            % 1e-9 decreases possible epsilon2 to 1e-8 without violating full coupling
            Prob.MIP.cpxControl.EPRHS = 1e-9; % min: 1e-9, max: 0.1, default: 1e-6
            tomrun = tomRun(solver, Prob, 0);
            if (tomrun.ExitFlag == 0)
                solution = getSolution(tomrun);
                temp = double(abs(solution.v) >= epsilon2);
                
                % check if there is a violation of full, partial,
                % directional, or anti-coupling
                for index=1:length(coupledI)
                    if ((coupledV(index) == 1 || coupledV(index) == 2) && (temp(coupledI(index)) ~= temp(coupledJ(index))))
                        save(strcat('Error021_', datestr(now, 'yyyy-mm-dd_HHMM'), '.mat'));
                        error('Sample %d violates coupling: coupling(%d,%d)=%d, w(%d,%d)=%d, w(%d,%d)=%d.', i, coupledI(index), coupledJ(index), coupledV(index), i, coupledI(index), temp(coupledI(index)), i, coupledJ(index), temp(coupledJ(index)));
                    elseif (coupledV(index) == 3 && temp(coupledI(index)) == 1 && temp(coupledJ(index)) == 0)
                        save(strcat('Error022_', datestr(now, 'yyyy-mm-dd_HHMM'), '.mat'));
                        error('Sample %d violates coupling: coupling(%d,%d)=%d, w(%d,%d)=%d, w(%d,%d)=%d.', i, coupledI(index), coupledJ(index), coupledV(index), i, coupledI(index), temp(coupledI(index)), i, coupledJ(index), temp(coupledJ(index)));
                    elseif (coupledV(index) == 5 && temp(coupledI(index)) == 0 && temp(coupledJ(index)) == 0)
                        save(strcat('Error023_', datestr(now, 'yyyy-mm-dd_HHMM'), '.mat'));
                        error('Sample %d violates coupling: coupling(%d,%d)=%d, w(%d,%d)=%d, w(%d,%d)=%d.', i, coupledI(index), coupledJ(index), coupledV(index), i, coupledI(index), temp(coupledI(index)), i, coupledJ(index), temp(coupledJ(index)));
                    end
                end
                
                % additional validation: check if w(i,:) has a valid v with x=0 <=> v=0 and Sv=0
                v = tom('v', n, 1);
                x = tom('x', numExchange, 1, 'integer');
                v_forward = tom('v_forward', n, 1);
                v_reversed = tom('v_reversed', n, 1);
                x_forward = tom('x_forward', n, 1, 'integer');
                x_reversed = tom('x_reversed', n, 1, 'integer');
                constraints = {...
                    S*v == 0;
                    lowerBound <= v <= upperBound;
                    v == v_forward - v_reversed;
                    epsilon2*x_forward <= v_forward <= x_forward*maxBound;
                    epsilon2*x_reversed <= v_reversed <= x_reversed*maxBound;
                    x_forward + x_reversed <= 1;
                    epsilon*x <= v_forward(exchange) + v_reversed(exchange) <= maxBound*x;
                    sum(x) >= 1;
                    temp == x_forward + x_reversed;
                    0 <= x_forward <= 1;
                    0 <= x_reversed <= 1;
                    0 <= v_forward;
                    0 <= v_reversed;
                    0 <= x <= 1;
                };
                Prob = sym2prob('mip', [], constraints, [], options);
                % decreases possible epsilon2 to 1e-8 without violating full coupling
                Prob.MIP.cpxControl.EPRHS = 1e-9; % min: 1e-9, max: 0.1, default: 1e-6
                if (strcmp(name, 'A_niger'))
                    % A_niger doesn't work with EPRHS = 1e-9
                    Prob.MIP.cpxControl.EPRHS = 1e-3;
                end
                tomrun = tomRun(solver, Prob, 0);
                if (tomrun.ExitFlag ~= 0)
                    % try different solver algorithm and scaling
                    Prob.MIP.cpxControl.LPMETHOD = 1;
                    Prob.MIP.cpxControl.SCAIND = 0;
                    tomrun = tomRun(solver, Prob, 0);
                    if (tomrun.ExitFlag ~= 0)
                        Prob.MIP.cpxControl.LPMETHOD = 1;
                        Prob.MIP.cpxControl.SCAIND = 1;
                        tomrun = tomRun(solver, Prob, 0);
                        if (tomrun.ExitFlag ~= 0)
                            Prob.MIP.cpxControl.LPMETHOD = 0;
                            Prob.MIP.cpxControl.SCAIND = -1;
                            tomrun = tomRun(solver, Prob, 0);
                            if (tomrun.ExitFlag ~= 0)
                                Prob.MIP.cpxControl.LPMETHOD = 0;
                                Prob.MIP.cpxControl.SCAIND = 0;
                                tomrun = tomRun(solver, Prob, 0);
                                if (tomrun.ExitFlag ~= 0)
                                    Prob.MIP.cpxControl.LPMETHOD = 0;
                                    Prob.MIP.cpxControl.SCAIND = 1;
                                    tomrun = tomRun(solver, Prob, 0);
                                    if (tomrun.ExitFlag ~= 0)
                                        % give up
                                        save(strcat('Error019_', datestr(now, 'yyyy-mm-dd_HHMM'), '.mat'));
                                        error('Could not find valid w (sampling mode %d, sample %d). tomRun error. ExitFlag: %d. %s\n', mode, i, tomrun.ExitFlag, tomrun.ExitText);                    
                                    end
                                end
                            end
                        end
                    end
                    Prob.MIP.cpxControl.LPMETHOD = 1;
                    Prob.MIP.cpxControl.SCAIND = -1;
                end
                
                % progress output
                if (charCount > 0)
                    for k=1:charCount
                        fprintf('\b \b');
                    end
                    charCount = 0;
                end
                
                % store the flux configuration
                w(i, :) = temp;
                
                alwaysZero = find(sum(w(1:i,:), 1) == 0);
                message = sprintf('Sample: %d. Distance to average: %g. Reactions always zero: %d.', i, sum(abs((avg+1)/2-(temp'))), length(alwaysZero));
                fprintf(message);
                charCount = charCount + length(message);
				
                % error checking
%                 if (max(abs(solution.v-(solution.v_forward-solution.v_reversed)) > 1e-4))
% 					violating = find(abs(solution.v-(solution.v_forward-solution.v_reversed)) == max(abs(solution.v-(solution.v_forward-solution.v_reversed))));
% 					j = violating(1);
% 					fprintf('v unequals v_forward-v_reversed: v(%d)=%g, v_forward(%d)=%g, v_reversed(%d)=%g.\n', j, solution.v(j), j, v_forward(j), j, v_reversed(j));
% 					charCount = 0;
%                 end
%                 if (length(find((abs(solution.x_forward + solution.x_reversed) >= epsilon2) ~= (abs(solution.v) >= epsilon2))) > 0)
% 					violating = find((solution.x_forward + solution.x_reversed) ~= (abs(solution.v) >= epsilon2));
% 					j = violating(1);
% 					fprintf('v violates x_forward-x_reversed: v(%d)=%g, x_forward(%d)=%g, x_reversed(%d)=%g.\n', j, solution.v(j), j, solution.x_forward(j), j, solution.x_reversed(j));
% 					charCount = 0;
%                 end
            else
                save(strcat('Error017_', datestr(now, 'yyyy-mm-dd_HHMM'), '.mat'));
				error('tomRun error. ExitFlag: %d. %s\n', tomrun.ExitFlag, tomrun.ExitText);
            end
        end
    end
    
    if (debug)
        % determine the average manhattan distance
        dist = 0;
        for i=1:(m-1)
            for j=(i+1):m
                dist = dist + sum(abs(w(i,:)-w(j,:)));
            end
        end
    end
    
    if (mode == 1)
        dlmwrite(strcat(name, '.samples'), w, '\t');
    else
        dlmwrite(strcat(name, '.samples', num2str(mode)), w, '\t');
    end
    
	t = toc;
    
    if (debug)
        fprintf('\nSampled %d flux configurations with an average of %g%% active reactions, average manhattan distance: %g.\n', m, (length(find(w))*100/numel(w)), dist/(m*(m-1)/2));
        alwaysZero = find(sum(w, 1) == 0);
        alwaysOne = find(sum(w, 1) == m);
        fprintf('%d reactions (%g%%) are inactive in all samples. %d reactions (%g%%) are active in all samples. Time: %g s.\n', length(alwaysZero), length(alwaysZero)*100/n, length(alwaysOne), length(alwaysOne)*100/n, t);
    else
        fprintf(' Sampling time: %g s.\n', t);
    end
end

