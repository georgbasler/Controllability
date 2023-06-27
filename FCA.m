function [S, rev, Reactions, fctable, blocked] = FCA(name, varargin)
% Calculate blocked reactions, full, partial, and directional coupling for
% the given network with 'name', stoichiometric matrix 'S', 'Reactions',
% and 'Metabolites'.

	% parse the arguments
	arg = 1;
	while (arg+1 <= nargin)
        if (strcmp(inputname(arg+1), 'S'))
            S = varargin{arg};
		elseif (strcmp(inputname(arg+1), 'rev'))
            rev = varargin{arg};
		elseif (strcmp(inputname(arg+1), 'Reactions'))
            Reactions = varargin{arg};
		elseif (strcmp(inputname(arg+1), 'Metabolites'))
            Metabolites = varargin{arg};
		end
        arg = arg+1;
    end
    
	% load the data
	if (~exist('S', 'var'))
		S = load(strcat(name, '.S'));
	end
	m = size(S, 1);
	n = size(S, 2);	
	if (~exist('rev', 'var'))
		rev = load(strcat(name, '.rev'));
    end
	if (~exist('Reactions', 'var'))
        if (exist(strcat(name, '.Reactions'), 'file') == 2)
% 			Reactions = importdata(strcat(name, '.Reactions'));
            fid = fopen(strcat(name, '.Reactions'));
            Reactions = textscan(fid, '%s');
            fclose(fid);
            Reactions = Reactions{:};
        else
			Reactions = strsplit(sprintf('R%d ', (1:n)));
        end
	end
	if (~exist('Metabolites', 'var'))
        if (exist(strcat(name, '.Metabolites'), 'file') == 2)
% 			Metabolites = importdata(strcat(name, '.Metabolites'));
            fid = fopen(strcat(name, '.Metabolites'));
            Metabolites = textscan(fid, '%s');
            fclose(fid);
            Metabolites = Metabolites{:};            
        else
			Metabolites = strsplit(sprintf('M%d ', (1:m)));
        end
	end
	
    fprintf('Running FCA on %s.', name);
    solver = 'glpk';
	network = struct('stoichiometricMatrix', S, 'reversibilityVector', rev, 'Reactions', char(Reactions), 'Metabolites', char(Metabolites));
	%     [fctable, blocked] = F2C2(solver, network);
    % suppress output
	[~, fctable, blocked] = evalc('F2C2(solver, network)');
    if (length(find(fctable < 0)) > 0)
        fprintf('Warning: fctable contains %d negative elements.\n', length(fctable < 0));
        fctable = abs(fctable);
    end

    % save flux coupling table and blocked reactions to files
	dlmwrite(strcat(name, '.fctable'), fctable, '\t');	
	dlmwrite(strcat(name, '.blocked'), blocked, '\t');

    % return the modified S, rev and Reactions for further calculations
    S = S(:, ~blocked);
    rev = rev(~blocked);
    Reactions = Reactions(~blocked);
    
    fprintf(' %d blocked reactions (%g%%).\n', length(find(blocked)), length(find(blocked))*100/length(blocked));
end

