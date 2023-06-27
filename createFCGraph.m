function [edges, profile] = createFCGraph ( name, varargin )
% Calculates the flux coupling profile and flux coupling graph.

	% parse the arguments
	arg = 1;
    while (arg+1 <= nargin)
        if (strcmp(inputname(arg+1), 'rev'))
            rev = varargin{arg};
        elseif (strcmp(inputname(arg+1), 'coupling'))
            coupling = varargin{arg};
        elseif (strcmp(inputname(arg+1), 'blocked'))
            blocked = varargin{arg};
        elseif (strcmp(inputname(arg+1), 'Reactions'))
            Reactions = varargin{arg};            
        elseif (strcmp(inputname(arg+1), 'ylimit'))
            ylimit = varargin{arg};            
        end        
        arg = arg+1;
    end
    
    % load the coupling data from inhibitors table
    if (~exist('coupling', 'var'))
        coupling = load(strcat(name, '.coupling'));
    end
    if (~exist('blocked', 'var'))
        blocked = load(strcat(name, '.blocked'));
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
    
    % check coupling table and reaction names for consistent lengths
    if (size(coupling, 1) ~= length(Reactions) || size(coupling, 2) ~= length(Reactions))
        error('Size of coupling table (%dx%d) does not match number of reaction names (%d).', size(coupling, 1), size(coupling, 2), length(Reactions));
    end
    % check coupling table and reversibility vector for consitent lengths
    if (size(coupling, 1) ~= length(rev) || size(coupling, 2) ~= length(rev))
        error('Size of coupling table (%dx%d) does not match size of reversibility vector (%d).', size(coupling, 1), size(coupling, 2), length(rev));
    end
    % check coupling table and blocked reactions for consitent lengths
    if (size(coupling, 1) ~= length(find(~blocked)) || size(coupling, 2) ~= length(find(~blocked)))
        error('Size of coupling table (%dx%d) does not match number of unblocked reactions in blocked vector (%d).', size(coupling, 1), size(coupling, 2), length(find(~blocked)));
    end
    
    fprintf('Creating FC graph for %s...', name);
    
    % determine and remove essential reactions by subtracting the 7s from
    % the diagonal
	essential = find(diag(coupling)==7);
    for i=1:length(essential)
        coupling(essential(i), essential(i)) = coupling(essential(i), essential(i)) -7;
    end
    
    % determine uncoupled pairs here for later error checking
    % note that removing 4's here would neglect those with additional
    % inhibitor coupling
    uncoupled = length(find(coupling==0)) + length(find(coupling==4)) - length(essential);
    multiples = 0;

	% set the diagonal to zero to remove implicit self-coupling
	coupling(find(eye(length(coupling)))) = 0;
    
	% get the non-zero values v and their indices (i,j)
	[i,j,v] = find(coupling);
	table = [i,j,v];
    % actual size of edges may be smaller, but difficult to determine
    edges = zeros(size(table, 1) + length(find(v>7)), 3);
    
    % replace summed couplings by multiple edges
    offset = 0;
    for k=1:size(table, 1)
        value = v(k);
        % only inhibitor coupling (8, 16, 32, 64) occurs together with other types
        if (value > 7)
            bits = bitget(value, 7:-1:1);
            inhibitor = bin2dec(num2str([bits(1:4) 0 0 0]));
            % create an edge for inhibitors of irreversible reactions and
            % those affecting both directions of reversible reactions
            value = value - inhibitor;
            % both directions affected by either source direction: ((bits(4) || bits(2)) && (bits(3) || bits(1)))
            % here: both directions affected by the same source direction (includes both to both)
            if (~rev(j(k)) || ((bits(1) && bits(2)) || (bits(3) && bits(4))))
                edges(k+offset, :) = [i(k) j(k) 6];
                offset = offset + 1;
                % count multiple edges
                if (value > 0 && value ~= 4)
                    multiples = multiples + 1;
                end
            elseif (value == 0 || value == 4)
                % count the new uncoupled pairs
                uncoupled = uncoupled + 1;
            end
        end
        if (value < 0 || value > 6)
            error('Invalid coupling value: %d at coupling(%d,%d), table(%d)=[%d,%d,%d], edges(%d)=[%d,%d,%d].', value, i(k), j(k), k, i(k), j(k), v(k), k+offset, edges(k+offset, 1), edges(k+offset, 2), edges(k+offset, 3));
        end
        % create an edge for the other coupling types
        % skip reversed directional coupling (v_i=0 => v_j=0)
        % here we may get empty entries in edges from value=value-inhibitor=0 and value==4
        if (value ~= 0 && value ~= 4)
            edges(k+offset, :) = [i(k) j(k) value];
        end
    end
    
    % remove empty entries (from removing inhibitors and v==4)
    edges(edges(:,1)==0, :) = [];
    edges(edges(:,2)==0, :) = [];
    edges(edges(:,3)==0, :) = [];
    
   	n = size(coupling, 1);
	pairs = n*(n-1); % number of ordered pairs involving different reactions (including essential reactions)
    couplingTypes = edges(:,3);
    if (length(find(couplingTypes>6)) > 0)
        error('edges contains %d illegal values (probably left over from inhibitor coupling).', length(find(couplingTypes>6)));
    end
    
    fprintf('done. Creating profile...');
    
    % frequency and probability of essential reactions
    freqEssential = length(essential);
    probEssential = freqEssential/n;
    % frequencies of coupling types
	fully = length(find(couplingTypes==1));
	partial = length(find(couplingTypes==2));
	directional = length(find(couplingTypes==3));
	anti = length(find(couplingTypes==5));
    inhibitor = length(find(couplingTypes==6));
    % probabilities of coupling types    
	probUncoupled = uncoupled / pairs;
    probMultiples = multiples / pairs;
	probFully = fully / pairs;
	probPartial = partial / pairs;
	probDirectional = directional / pairs;
	probAnti = anti / pairs;
    probInhibitor = inhibitor / pairs;
    
	% uncoupled pairs: the remaining ones minus multiple edges
	if (pairs ~= (uncoupled + fully + partial + directional + anti + inhibitor - multiples))
        error('Uncoupled+coupled should equal number of pairs.\npairs=%d, uncoupled=%d, coupled=%d (fully=%d, partial=%d, directional=%d, anti=%d, inhibitor=%d), multiples=%d.\n', pairs, uncoupled, fully+partial+directional+anti+inhibitor, fully, partial, directional, anti, inhibitor, multiples);
	end
	if (abs(probUncoupled + probFully + probPartial + probDirectional + probAnti + probInhibitor - probMultiples - 1) > 0.0001)
        fprintf('pairs=%d, uncoupled=%d, coupled=%d (fully=%d, partial=%d, directional=%d, anti=%d, inhibitor=%d), multiples=%d.\n', pairs, uncoupled, fully+partial+directional+anti+inhibitor, fully, partial, directional, anti, inhibitor, multiples);
		error('%d + %d + %d + %d + %d + %d - %d = %d (should be 1)', probUncoupled, probFully, probPartial, probDirectional, probAnti, probInhibitor, probMultiples, probUncoupled + probFully + probPartial + probDirectional + probAnti + probInhibitor - probMultiples);
	end
	
    % write the profile as struct to file
	frequencies = struct('Uncoupled', uncoupled, 'Full', fully, 'Partial', partial, 'Directional', directional, 'Anti', anti, 'Inhibitor', inhibitor, 'Essential', freqEssential);
	probabilities = struct('Uncoupled', probUncoupled, 'Full', probFully, 'Partial', probPartial, 'Directional', probDirectional, 'Anti', probAnti, 'Inhibitor', probInhibitor, 'Essential', probEssential);
    if (exist(strcat(name, '.profile.mat'), 'file') == 2)
        delete(strcat(name, '.profile.mat'));
    end
    profile = struct();
    profile.frequencies = frequencies;
    profile.probabilities = probabilities;
    save(strcat(name, '.profile.mat'), 'profile');
    % write the profiles as plain text to file
    if (exist(strcat(name, '.profile'), 'file') == 2)
        delete(strcat(name, '.profile'));
	end
	fText = evalc('frequencies');
    pText = evalc('probabilities');
	fid3 = fopen(strcat(name, '.profile'), 'w', 'n', 'UTF-8'); % encoding required to ensure proper line breaks
    fprintf(fid3, fText);
	fprintf(fid3, pText);
	fclose(fid3);
	
    fprintf('done. Plotting...');
    
    % plot the normalized distribution of coupling types (flux coupling profile) and save as PNG file
    total = (probFully+probPartial+probDirectional+probAnti+probInhibitor);
    graph = plot(linspace(0.5, 5.5, 5), [probFully/total probPartial/total probDirectional/total probAnti/total probInhibitor/total], 'b-o');
    set(gca,'XTick', linspace(0.5, 5.5, 5));
	set(gca,'XTickLabel', {'Full' 'Partial' 'Directional' 'Anti' 'Inhibitor'});
    ylabel('Frequency');
	title(strrep(name, '_', ' '));
    if (exist('ylimit', 'var'))
        ylim(ylimit);
    end
    saveas(graph, strcat(name, '.png'));
    
    fprintf('done. Writing graphs...');
    
	% print the graph to file
    file = strcat(name, '.fcgraph');
	fid = fopen(file, 'w', 'n', 'UTF-8'); % encoding required to ensure proper line breaks
	if (fid == -1)
		error('Cannot open file: %s', file);
    end
	for k=1:size(edges,1)
		fprintf(fid, '%s\t%s\t%d\n', Reactions{edges(k,1)}, Reactions{edges(k,2)}, edges(k,3));
	end
	fclose(fid);
    % print a separate graph with split source reactions
    file2 = strcat(name, '.split.fcgraph');
	fid2 = fopen(file2, 'w', 'n', 'UTF-8'); % encoding required to ensure proper line breaks
	if (fid2 == -1)
		error('Cannot open file: %s', file2);
    end
	for k=1:size(edges,1)
        i = edges(k,1);
        j = edges(k,2);
        value = coupling(i,j);
        if (rev(i) && value > 7)
            bits = bitget(value, 7:-1:1);
            if (bits(1) || bits(2))
                % i is forward: add an edge with [+]
                fprintf(fid2, '%s[+]\t%s\t%d\n', Reactions{edges(k,1)}, Reactions{edges(k,2)}, edges(k,3));
            end
            if (bits(3) || bits(4))
                % i is reversed: add an edge with [-]
                fprintf(fid2, '%s[-]\t%s\t%d\n', Reactions{edges(k,1)}, Reactions{edges(k,2)}, edges(k,3));
            end
        else
            % irreversible or not inhibitor coupled
            fprintf(fid2, '%s\t%s\t%d\n', Reactions{edges(k,1)}, Reactions{edges(k,2)}, edges(k,3));
        end
    end
    fclose(fid2);
    
    fprintf('done.\n');
end
