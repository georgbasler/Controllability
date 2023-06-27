function [ S, rev, blocked, Reactions, Metabolites ] = loadData( name )
% Loads the stoichiometric matrix S, rev, blocked, Reactions, and
% Metabolites.

    blocked = load(strcat(name, '.blocked'));
    S = load(strcat(name, '.S'));
    S = S(:,~blocked);
    rev = load(strcat(name, '.rev'));
    rev = rev(~blocked);
%     Reactions = importdata(strcat(name, '.Reactions'));
    fid = fopen(strcat(name, '.Reactions'));
    Reactions = textscan(fid, '%s');
    fclose(fid);
    Reactions = Reactions{:};
    Reactions = Reactions(~blocked);
%     Metabolites = importdata(strcat(name, '.Metabolites'));
    fid = fopen(strcat(name, '.Metabolites'));
    Metabolites = textscan(fid, '%s');
    fclose(fid);
    Metabolites = Metabolites{:};  
end

