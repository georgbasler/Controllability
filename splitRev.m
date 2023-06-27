function [S_split, rev_split] = splitRev(S, rev, factor)
% Creates a new matrix from S by adding columns for each non-zero entry in
% rev.
% factor is the value by which the previous (unsplit) column is multiplied
% to yield the new column, e.g. -1: negative (stoichiometric matrix), 0:
% zero, and 1: copied.

    if (size(S, 2) ~= length(rev))
        error('Number of columns of S (%d) must match length of rev (%d). Make sure to remove blocked reactions.', size(S, 2), length(rev));
    end
    
	m = size(S,1);
	n = size(S,2);
	S_split = zeros(m, n + length(find(rev)));
    rev_split = zeros(length(rev) + length(find(rev)), 1);
	offset = 0;

	for i=1:n
		% copy the forward reaction to the corresponding column
		S_split(:,i+offset) = S(:,i);
		if (rev(i) ~= 0)
            rev_split(i+offset) = 1;
			offset = offset+1;
            if (factor)
                % the new column is the old column multiplied by factor
                S_split(:,i+offset) = S(:,i)*factor;
            end
		end
	end
end
