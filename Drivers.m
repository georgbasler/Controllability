function [ d, time ] = Drivers( name, w, specified, varargin )
% Calculates the drivers for the given flux target state w.
    
    % parse the arguments
    check = 0;
    arg = 1;
    while (arg+3 <= nargin)
        if (strcmp(inputname(arg+3), 'rev'))
            rev = varargin{arg};
        elseif (strcmp(inputname(arg+3), 'S'))
            S = varargin{arg};            
        elseif (strcmp(inputname(arg+3), 'blocked'))
            blocked = varargin{arg};
        elseif (strcmp(inputname(arg+3), 'coupling'))
            coupling = varargin{arg};
        elseif (strcmp(inputname(arg+3), 'solver'))
            solver = varargin{arg};
        elseif (ischar(varargin{arg}))
            if (strcmp(varargin{arg}, 'check'))
                check = 1;
            end
        end            
        arg = arg+1;
    end
    
    if (~exist('solver', 'var'))
        solver = 'cplex';
    end
    
    if (check)
        if (~exist('S', 'var'))
            S = loadData(name);
        end
        n = size(S, 2);
        
        epsilon = 1;
        epsilon2 = 1e-7; % smaller than 1e-8 gives valid solution for violating full coupling
        lowerBound = zeros(n, 1);
        lowerBound(find(rev)) = -1000;
        upperBound = repmat(1000, n, 1);
        maxBound = max(upperBound);
        exchange = (sum(abs(S), 1) == sum(S, 1) | sum(abs(S), 1) == -sum(S, 1));
        numExchange = length(find(exchange));
        
        x = tom('x', numExchange, 1, 'integer');
        v = tom('v', n, 1);
        v_forward = tom('v_forward', n, 1);
        v_reversed = tom('v_reversed', n, 1);
        x_forward = tom('x_forward', n, 1, 'integer');
        x_reversed = tom('x_reversed', n, 1, 'integer');
        v_forward_ex = tom('v_forward_ex', numExchange, 1);
        v_reversed_ex = tom('v_reversed_ex', numExchange, 1);   
        
%         specified = find(~isnan(w));
        
        constraints = {...
                S*v == 0;
                lowerBound <= v <= upperBound;
                v == v_forward - v_reversed;
                epsilon2*x_forward <= v_forward <= x_forward*maxBound;
                epsilon2*x_reversed <= v_reversed <= x_reversed*maxBound;
                x_forward + x_reversed <= 1;
                epsilon*x <= v_forward(exchange) + v_reversed(exchange) <= maxBound*x;
                sum(x) >= 1;
                w(specified)' == x_forward(specified) + x_reversed(specified);
                0 <= x_forward <= 1;
                0 <= x_reversed <= 1;
                0 <= v_forward;
                0 <= v_reversed;
                0 <= x <= 1;
        };
        options = struct;
        options.name = name;
        options.solver = solver;    
        Prob = sym2prob('mip', [], constraints, [], options);
%         Prob.MIP.cpxControl.LPMETHOD = 1;
%         Prob.MIP.cpxControl.SCAIND = -1;
%         Prob.MIP.cpxControl.EPINT = 1e-9; % min: 1e-9, max: 0.5, default: 1e-5
%         Prob.MIP.cpxControl.EPMRK = 0.99999; % min: 0.0001, max: 0.99999, default: 0.1
         Prob.MIP.cpxControl.EPRHS = 1e-9; % min: 1e-9, max: 0.1, default: 1e-6
%         Prob.MIP.cpxControl.NETEPRHS = 1e-11; % min: 1e-1, max: 1e-11, default: 1e-6
        
%         Prob = mipAssign(c, S, b_L, b_U, x_L, x_U, x_0, name, setupFile, nProblem, IntVars); %..., VarWeight, KNAPSACK, fIP, xIP, f_Low, x_min, x_max, f_opt, x_opt);
    %     Prob.MIP.cpxControl.LPMETHOD = 4;
        tomrun = tomRun('cplex', Prob, 1);
        fprintf('tomrun.ExitFlag = %d\n', tomrun.ExitFlag);
        if (tomrun.ExitFlag == 0)
            fprintf('Valid Sv=0 found for w.\n');
             solution = getSolution(tomrun);
             fprintf('w(specified)=%d,%d\n', w(specified));
             fprintf('solution.v(specified)=%d,%d\n', solution.v(specified)');
             fprintf('solution.v_forward(specified)=%d,%d\n', solution.v_forward(specified)');
             fprintf('solution.v_reversed(specified)=%d,%d\n', solution.v_reversed(specified)');
             fprintf('solution.x_forward(specified)=%d,%d\n', solution.x_forward(specified)');
             fprintf('solution.x_reversed(specified)=%d,%d\n', solution.x_reversed(specified)');
        elseif (tomrun.ExitFlag == 4)
            error('No Sv=0 found for w!');
        else
            save(strcat('Error020_', datestr(now, 'yyyy-mm-dd_HHMM'), '.mat'));
            error('tomRun error while checking for a solution of w with w=abs(sign(v)) and Sv=0.\nExitFlag: %d. %s\n', tomrun.ExitFlag, tomrun.ExitText);
        end
    end
    
    % calculate the activators and deactivators for w
    if (exist('rev', 'var') && exist('blocked', 'var') && exist('coupling', 'var'))
        [A, controllersTime] = Controllers(name, w, specified, rev, blocked, coupling);
    else
        [A, controllersTime] = Controllers(name, w, specified);
    end
    
    t = tic;
    n = size(A, 1);
%     coupling = load(strcat(name, '.coupling'));

    % mipAssign is 2-10 times faster than sym2prob
    % A: matrix A in Ax=b
    c = ones(n, 1); % objective: minimize sum of non-zero entries in x (min_x c'*x)
    x_L = zeros(n, 1); % lower bound of x: 0
    x_U = ones(n, 1); % upper bound of x: 1
    b_L = ones(n, 1); % lower bound of b: 1. (Ax >= 1)
    b_U = repmat(n, n, 1); % upper bound of b: n. (Ax <= n => 1 <= Ax <= n)
    x_0 = [];
    setupFile = [];
    nProblem = 1;
    IntVars = 1:n;
    
    % print w and control matrix
%     w
%     A(1:5,1:5)
    
    % OPTIONAL: this excludes exchange reactions as drivers
%     if (~exist('S', 'var'))
%         S = loadData(name);
%     end
%     exchange = (sum(abs(S), 1) == sum(S, 1) | sum(abs(S), 1) == -sum(S, 1));
%     A(:,exchange) = 0;
%     b_L(exchange, 1) = 0;
    % provide an initial solution (allows to check whether a
    % solution is returned valid). Note: this must be defined as logical.
    % requires to exclude exchange reactions (above)
%     drivers = [0 1 0 0 1 0 0 0 0];
%     x_L(logical(drivers), 1) = 1;
    
    Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, x_0, name, setupFile, nProblem, IntVars); %..., VarWeight, KNAPSACK, fIP, xIP, f_Low, x_min, x_max, f_opt, x_opt);
%     Prob.MIP.cpxControl.LPMETHOD = 4;
    tomrun = tomRun('cplex', Prob, 0);
    % 'mipSolve' seems to be a bit slower than 'cplex'
% 	solver = 'mipSolve';
%     [~, tomrun] = evalc('tomRun(solver, Prob, 0)');
    
    if (tomrun.ExitFlag == 0)
        d = tomrun.x_k;
%         fprintf('Found minimal set of %d drivers.\n', length(find(d)));
    else
        save(strcat('Error016_', datestr(now, 'yyyy-mm-dd_HHMM'), '.mat'));
        error('tomRun error while calculating drivers for %s.\nExitFlag: %d. %s\n', name, tomrun.ExitFlag, tomrun.ExitText);
    end
    
    time = toc(t) + controllersTime;
    
end

