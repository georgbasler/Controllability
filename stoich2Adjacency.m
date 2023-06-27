function [ A, substrateMatrix, productMatrix] = stoich2Adjacency( S, rev )
% Converts the given stoichiometric matrix S into the corresponding adjacency
% matrix A of the bipartite graph:
% 
%     [ 0 0  S S S ]
%     [ 0 0  S S S ]
% A = [ P'P' 0 0 0 ]
%     [ P'P' 0 0 0 ]
%     [ P'P' 0 0 0 ]

m = size(S, 1);
n = size(S, 2);

% create the adjacency matrix of substrate-reaction relations
substrateMatrix = zeros(m, n);
substrateMatrix(S<0) = 1;
% create the adjacency matrix of reaction-product relations
productMatrix = zeros(m, n);
productMatrix(S>0) = 1;

% switch the substrates and products of reversible reactions
S(:,logical(rev)) = S(:,logical(rev))*-1;
% add the switched substrates and products
substrateMatrix(S<0) = 1;
productMatrix(S>0) = 1;

A = [zeros(m, m) substrateMatrix; productMatrix' zeros(n, n)];

end

