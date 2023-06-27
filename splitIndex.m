function [ splitIndex ] = splitIndex( unsplitIndex, rev )
% calculates the index of splitIndex in the corresponding split matrix by rev
    if (unsplitIndex == 1)
        splitIndex = 1;
    else
        splitIndex = unsplitIndex + length(find(rev(1:unsplitIndex-1, 1)));
    end
end
