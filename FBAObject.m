classdef FBAObject < handle

%  	properties (Constant)
%  	end
    
	% set by the constructor, get by anyone
	% SetAccess=immutable seems not to be supported by Matlab R2009a
	properties (SetAccess=private, GetAccess=public)
		name
		S
		M
		T
		C
		metabolites
		reactions
		m
		n
	end
	
	% set by the class (and subclasses), get by anyone
	properties (SetAccess=protected, GetAccess=public)
		B
	end

%  	properties (Dependent)
%  		<variables...>
%  	end
%  	properties
%  		<variables...>
%  	end
	
	methods
		
        % constructor arguments: name, S, [metabolites, reactions, M, T, B, C]
		function FBAObject = FBAObject(name, S, varargin)
            % constructor: create a new FBAObject from the provided data
            FBAObject.name = name;
            if (ischar(S))
                S = importdata(S);
            end
            FBAObject.S = S;
            FBAObject.m = size(S, 1);
			FBAObject.n = size(S, 2);
            arg = 1;
            
            clear metabolites; clear reactions; clear M; clear T; clear B; clear C; clear rev;
            
			% parse optional input arguments
            while (arg+2 <= nargin)
                if (ischar(varargin{arg}))
                    % extract variable name and value from the argument given as file name
                    [~,~,ext] = fileparts(varargin{arg});
                    variable = ext(2:length(ext));
%                     value = importdata(varargin{arg});
                    fid = fopen(varargin{arg});
                    value = textscan(fid, '%s');
                    fclose(fid);
                    value = value{:};
                else
                    % extract variable name and value from the argument given as variable
                    variable = inputname(arg+2);
                    value = varargin{arg};
                end
                % temporarily store in another struct to allow other
                % arguments than the fields of FBAObject (e.g. rev)
                data.(variable) = value;
                arg = arg+1;
            end
            
            % check and complete input data
            % metabolites
            if (isfield(data, 'metabolites'))
                if (size(data.metabolites, 1) == 1)
                    data.metabolites = data.metabolites';
                end
                if (size(data.metabolites, 1) ~= FBAObject.m)
                    error('Size of metabolites (%d) must correspond to rows in S (%d).', size(data.metabolites, 1), FBAObject.m);
                end
                FBAObject.metabolites = data.metabolites;
            else
                FBAObject.metabolites = strsplit(sprintf('M%d ', (1:FBAObject.m)));
            end
            % reactions
            if (isfield(data, 'reactions'))
                if (size(data.reactions, 2) == 1)
                    data.reactions = data.reactions';
                end
                if (size(data.reactions, 2) ~= FBAObject.n)
                    error('Size of reactions (%d) must correspond to columns in S (%d).', size(data.reactions, 2), FBAObject.n);
                end
                FBAObject.reactions = data.reactions;
            else
                % no reaction names: generate a string sequence
                FBAObject.reactions = strsplit(sprintf('R%d ', (1:FBAObject.n)));
            end
            % M
            if (isfield(data, 'M'))
                if (size(data.M,1) ~= FBAObject.m)
                    error('Size of M (%d) must correspond to number of metabolites (%d).', size(data.M,1), FBAObject.m);
                end
                FBAObject.M = data.M;
            end
            % T
            if (isfield(data, 'T'))
                if (size(data.T,1) ~= FBAObject.m)
                    error('Size of T (%d) must correspond to number of metabolites (%d).', size(data.T,1), FBAObject.m);
                end
                FBAObject.T = data.T;
            end
            % rev
            if (isfield(data, 'rev') && ~isfield(data, 'B'))
                if (size(data.rev,1) ~= FBAObject.n)
                    error('Size of rev (%d) must correspond to number of reactions (%d).', size(data.rev,1), FBAObject.n);
                end
                % create B from rev
                data.B = repmat([0; 1000], 1, FBAObject.n);
                data.B(1,find(data.rev)) = -1000;
                FBAObject.B = data.B;
            end
            % B
            if (isfield(data, 'B'))
                if (~(size(data.B,1) == 2 && size(data.B,2) == FBAObject.n) || (size(data.B,1) == FBAObject.n && size(data.B,2) == 2))
                    error('B (%dx%d) must be an array of size 2x%d).', size(data.B,1), size(data.B,2), FBAObject.n);
                end
                FBAObject.B = data.B;
            else
                error('Either flux bounds B or reversibilities rev must be given.');
            end
            % C
            if (isfield(data, 'C'))
                if (size(data.C,1) ~= FBAObject.m)
                    error('Size of C (%d) must correspond to the number of metabolites (%d).', size(data.C,1), FBAObject.m);
                end
                FBAObject.C = data.C;
            else 
                FBAObject.C = ones(FBAObject.m, 1);
            end
		end
		
		% returns the names of the given indices: metabolite names for type==1 or type=='metabolites', reaction names for type==2 or type=='reactions'
		function names = indexToName(obj, indices, type)
			if (type == 1)
				% metabolites
				names = obj.metabolites(indices);
			elseif (type == 2)
				% reactions
				names = obj.reactions(indices);
			end
		end
		
		% returns the indices of the given metabolite or reaction names
		% the type is inferred from the dimension of the name vector: column vector returns metabolites, row vector returns reactions.
		% for a single value (1x1), both metabolite and reaction names are searched
		function indices = nameToIndex(obj, names)
			if (size(names, 1) == 1)
				% reactions
				indices = arrayfun(@(x) find(strcmp(x,obj.reactions)), names);
			end
			if (size(names, 2) == 1 || isempty(indices))
				% metabolites
				indices = arrayfun(@(x) find(strcmp(x,obj.metabolites)), names);
			end
		end
	
	% adds the reactions with the specified bounds. only one set of bounds for all reactions is
	% supported.
	function addReactions(obj, reactions, bounds)
		if (~isstruct(reactions))
			error('Argument reactions must be a struct of reaction names and stoichiometric coefficients.');
		end
		if (size(obj.S, 1) ~= size(reactions.data, 1))
			error('Number of metabolites (%d) does not correspond to rows in S (%d).', size(reactions, 1), size(obj.S, 2));
		end
		if (~strcmp(class(bounds),'double') || size(bounds, 1) ~= 2 || size(bounds, 2) ~= 1)
			error('Argument bounds must be a double array of size (2x1).');
		end

		% update reactions, S, B, and n
		numReactions = size(reactions.data, 2);
		obj.reactions = [obj.reactions reactions.textdata];
		obj.S = [obj.S reactions.data];
		obj.B = [obj.B repmat(bounds, 1, numReactions)];
		obj.n = (obj.n + numReactions);
	end
	
	function removeReaction(obj, reaction)
		if (~ischar(reaction))
			error('Argument reaction must be a char.');
		end
		index = find(strcmp(obj.reactions, reaction));
		if (~isempty(index))
			obj.reactions(index) = [];
			obj.S(:, index) = [];
			obj.B(:, index) = [];
			obj.n = obj.n - 1;
		end
	end
	
%  	function addMetabolites(obj, metabolites, masses, energy)
%  		if (~isstruct(metabolites))
%  			error('Argument metabolites must be a struct of metabolite names and stoichiometric coefficients.');
%  		end
%  		if (size(obj.S, 2) ~= size(metabolites.data, 2))
%  			error('Number of reactions (%d) does not correspond to columns in S (%d).', size(metabolites, 2), size(obj.S, 1));
% 		end
% 		if (size(masses, 1) ~= size(metabolites.data, 1) || size(masses, 2) ~= size(obj.M, 2))
% 			error('Dimensions of masses must correspond to the number (%d) of metabolites and elements in M (%d).', size(metabolites.data, 1), size(obj.M, 2));
% 		end
%  
%  		% update metabolites, S, M, T, C, and m
%  		numMetabolites = size(metabolites.data, 1);
%  		obj.metabolites = [obj.metabolites; metabolites.textdata'];
%  		obj.S = [obj.S; metabolites.data];
% 		obj.M = [obj.M; masses];
% 		obj.T = [obj.T; sprintf(... T should be double. M and T should be optional.
%  		
% 	end
	
	function removeMetabolite(obj, metabolite)
		if (~ischar(metabolite))
			error('Argument metabolite must be a char.');
		end
		index = find(strcmp(obj.metabolites, metabolite));
		if (~isempty(index))
			obj.metabolites(index) = [];
			obj.S(index, :) = [];
			obj.M(index, :) = [];
			obj.T(index, :) = [];
			obj.C(index, :) = [];
			obj.m = obj.m - 1;
		end
	end
	
		
        function printEquation(obj, reaction)
            if (~ischar(reaction) && ~iscell(reaction))
				error('Argument reaction must be char or char array.');
            end
            
            rIndex = find(strcmp(obj.reactions, reaction));
			if (isempty(rIndex))
                error('No such reaction: %s\n', reaction);
		end
			
			obj.printEquations(rIndex);
		end
			
		function printEquations(obj, reactions)
            
			for i=1:size(reactions)
				
				rIndex = reactions(i);
			
				substrates = find(obj.S(:,rIndex) < 0);
				products = find(obj.S(:,rIndex) > 0);

				for j=1:size(substrates),
					fprintf('[%d] %s', -obj.S(substrates(j), rIndex), obj.metabolites{substrates(j)});
					if (j < size(substrates, 1))
						fprintf(' + ');
					end
				end
				if (obj.B(1, rIndex) < 0)
					fprintf(' <==> ');
				else
					fprintf(' --> ');
				end
				for j=1:size(products),
					fprintf('[%d] %s', obj.S(products(j), rIndex), obj.metabolites{products(j)});
					if (j < size(products, 1))
						fprintf(' + ');
					end
				end
				fprintf('\n');
			end
        end
		
		function setBounds(obj, reactions, bounds)
		% set the upper and lower flux bounds of the specified reaction(s)
			if (~ischar(reactions) && ~iscellstr(reactions))
  				error('Argument reactions to setBounds must be a char or char array of reaction names.');
  			elseif (~strcmp(class(bounds),'double') || size(bounds, 1) ~= 2 || size(bounds, 2) ~= 1)
				error('Argument bounds to setBounds must be a double array of size (2x1).');
  			elseif (ischar(reactions))
				reactions = {reactions};
			elseif (size(reactions,1) > 1)
				reactions = reactions';
			end
  			% get the column indices in B of the reactions specified in r
  			indices = arrayfun(@(x) find(strcmp(x,obj.reactions)), reactions);
  			% set the new values by replicating the values in b
  			obj.B(:,indices) = repmat(bounds, 1, size(indices,2));
		end
		
		function bounds = getBounds(obj, reactions)
		% get the upper and lower flux bounds of the specified reaction(s)
			if (~ischar(reactions) && ~iscellstr(reactions))
  				error('Argument to getBounds must be a char or char array of reaction names.');
			elseif (ischar(reactions))
				reactions = {reactions};
			end
			% get the columns in B specified by the indices of reactions in r
  			bounds = obj.B(:,arrayfun(@(x) find(strcmp(x,obj.reactions)), reactions));
		end
		
		function setAllBounds(obj, bounds)
		% set the upper and lower flux bounds of all reactions
			if (~strcmp(class(bounds),'double') || size(bounds, 1) ~= 2 || size(bounds, 2) ~= 1)
				error('Argument to setAllBounds must be a double array of size (2x1).');
			end
			obj.B = repmat(bounds, 1, size(obj.B,2));
		end
		
		function reactions = getImportReactions(obj)
			% get the names of all reactions with no substrates (only positive coefficients)
			indices = (sum(abs(obj.S),1) == sum(obj.S,1));
			reactions = obj.reactions(indices);
		end
		
		function reactions = getExportReactions(obj)
			% get the names of all reactions with no products (only negative coefficients)
			indices = (sum(-abs(obj.S),1) == sum(obj.S,1));
			reactions = obj.reactions(indices);
        end
		
        function deltaG = getDeltaGr(obj, reaction)
            
            % get the deltaG of the specified reaction(s)
			if (~ischar(reaction))
  				error('Argument to getDeltaGr must be a char.');
			end
                
			reaction = {reaction};
            
            rIndex = find(strcmp(obj.reactions, reaction));
            if (isempty(rIndex))
                error('No such reaction: %s\n', cell2str(reaction));
            end
            
            substrates = obj.S(:,rIndex);
            substrates(substrates>0) = 0;
            products = obj.S(:,rIndex);
            products(products<0) = 0;
            % set NaNs in T to 0 to allow vector multiplication
            T = obj.T;
        	T(isnan(T)) = 0;
    
            deltaG = T' * (products-substrates);
        end
        
        function mass = getMass(obj, metabolite)
        	mass = obj.M(find(strcmp(obj.metabolites, metabolite)),:);
        end
        
%  		function [r1,r2] = method(<obj>, <args...>)
%  		function prop = get.Property(<obj>)	% implicitly invoked by accessing the property
	
	% return the exchange reactions with non-zero flux in v
	function [reactions,fluxes] = getExchangeFluxes(obj, v)
		
		if (size(v,1) == size(obj.reactions,2) && size(v,2) == size(obj.reactions,1))
			v = v';
		end
		
		if (size(v) ~= size(obj.reactions))
			error('Dimension of v (%d,%d) does not correspond to reactions (%d,%d).', size(v,1), size(v,2), size(obj.reactions, 1), size(obj.reactions, 2));
		end
		
		% determine the indices of import reactions, non-zero fluxes, and their common non-zero indices
		exchangeIndices = (sum(abs(obj.S),1) == sum(obj.S,1));
		exchangeIndices = exchangeIndices | (sum(-abs(obj.S),1) == sum(obj.S,1));
		fluxIndices = (v ~= 0);
		exchangeFluxIndices = fluxIndices & exchangeIndices;
		
		reactions = obj.reactions(exchangeFluxIndices);
		fluxes = v(exchangeFluxIndices);
	end
	
		% return the import reactions with non-zero flux in v
	function [reactions,fluxes] = getFluxes(obj, v)
		
		if (size(v,1) == size(obj.reactions,2) && size(v,2) == size(obj.reactions,1))
			v = v';
		end
		
		if (size(v) ~= size(obj.reactions))
			error('Dimension of v (%d,%d) does not correspond to reactions (%d,%d).', size(v,1), size(v,2), size(obj.reactions, 1), size(obj.reactions, 2));
		end
		
		% determine the indices of import reactions, non-zero fluxes, and their common non-zero indices
		fluxIndices = (v ~= 0);
		
		reactions = obj.reactions(fluxIndices);
		fluxes = v(fluxIndices);
	end
	
	end
end
