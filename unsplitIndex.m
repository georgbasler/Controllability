function [ unsplitIndex ] = unsplitIndex( splitIndex, rev )
% calculates the index of splitIndex in the corresponding unsplit matrix by rev
	unsplitIndex = 0;
	i = 0;
	while (i < splitIndex)
		unsplitIndex = unsplitIndex + 1;
		if (rev(unsplitIndex))
			i = i + 1;
		end
		i = i + 1;
	end
end
