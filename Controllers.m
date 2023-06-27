function [ c, time ] = Controllers( name, sample, specified, varargin )
% Calculates the activators and deactivators for the given
% reaction activity pattern w. The n x n matrix c contains a row for each
% reaction in w, and a column indicating its activators, if w(i)=1,
% or inactivators, if w(i)=0.
% NOTE: here it is NOT checked whether w is feasible! If this is not the case,
% then the returned drivers are meaningless.

    % parse the arguments
    arg = 1;
    while (arg+3 <= nargin)
        if (strcmp(inputname(arg+3), 'rev'))
            rev = varargin{arg};
        elseif (strcmp(inputname(arg+3), 'blocked'))
            blocked = varargin{arg};
        elseif (strcmp(inputname(arg+3), 'coupling'))
            coupling = varargin{arg};
        end
        arg = arg+1;
    end
    
    if (~exist('blocked', 'var'))
        blocked = load(strcat(name, '.blocked'));
    end
    if (~exist('rev', 'var'))
        rev = load(strcat(name, '.rev'));
        rev = rev(~blocked);
    end
    if (~exist('coupling', 'var'))
        % use the full inhibitor coupling
        coupling = load(strcat(name, '.coupling'));
    end
    
    if (length(sample) ~= size(coupling, 1))
        error('Unmatching sizes of w (%d) and coupling (%d).', length(w), size(coupling, 1));
    end
    
    t = tic;
    
    n = size(coupling, 1); 
    c = zeros(n, n);
    % NOTE: here w is the sample containing NaNs for unspecified,
    % sample is the full sample
    w = nan(1, n);
    w(specified) = sample(specified);
    
    essential = find(diag(coupling)==7);
    % find(w==1) skips NaN
    active = setdiff(find(w==1), essential);
    inactive = find(w==0);
    % essential and unspecified reactions are set to 1,
    % indicating they can be 'controlled' by any reaction
    c(essential, :) = 1;
    c(find(isnan(w)), :) = 1;
    
    if (~isempty(intersect(essential, inactive)))
        error('Inactive essential reactions in w.');
    end
    
%     % remove inhibitor coupling to yield the other types of multiple couplings
%     [i,j] = find(coupling > 7);
%     for k=1:length(i)
%         value = coupling(i(k), j(k));
%         bits = bitget(value, 7:-1:1);
%         inhibitorType = bin2dec(num2str([bits(1:4) 0 0 0]));
%         coupling(i(k), j(k)) = value - inhibitorType;
%     end
    
    % collect the activators of active reactions:
    % 1. active fully, partially, and directionally coupled
    % 2. inactive anti-coupled.
    for i=1:length(active)
        index = active(i);
        % NaNs are accepted by using w~=0 and w~=1
        % but now we take unspecified drivers only with their status in the
        % full sample, as this guarantess steady-state compatibility and
        % that an unspecified driver can have only 1 status
        activators = (sample'== 1 & coupling(:, index) == 1);
        activators = (activators | sample'==1 & coupling(:, index) == 2);
        activators = (activators | sample'==1 & coupling(:, index) == 3);
        activators = (activators | sample'==0 & coupling(:, index) == 5);
        if (size(activators, 2) ~= 1)
            error('Wrong activators dimension.');
        end
        c(index,:) = activators';
    end
    % collect the inactivators of inactive reactions:
    % 1. inactive fully, partially, and  reversed directionally
    % 2. active inhibitor coupled.
    for i=1:length(inactive)
        index = inactive(i);
        % NaNs are accepted by using w~=1
        inactivators = (sample'==0 & coupling(:, index) == 1);
        inactivators = (inactivators | sample'==0 & coupling(:, index) == 2);
        inactivators = (inactivators | sample'==0 & coupling(:, index) == 4);
        if (size(inactivators, 2) ~= 1)
            error('Wrong inactivators dimension.');
        end
        % add the active inhibitors
        % NaNs are accepted by using w~=0 (but see comment above)
        inhibitors = find(sample'==1 & coupling(:, index) > 7);
        for k=1:length(inhibitors)
            bits = bitget(coupling(inhibitors(k), index), 7:-1:1);
            % both directions affected by the same source direction (includes both to both)
            if (~rev(index) || ((bits(1) && bits(2)) || (bits(3) && bits(4))))
                inactivators(inhibitors(k), 1) = 1;
            end
        end
        c(index,:) = inactivators';
    end
    
    time = toc(t);
end

