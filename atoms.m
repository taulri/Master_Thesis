function [AtomCoef] = atoms(Aindex, Acoef, Ach, Aloc)
% Unpack the index and coefficients of all the atoms

indexes = [];
coefficients = [];
ch = [];
loc = [];

for ind = 1:size(Aindex,2)
    indexes = [indexes,Aindex{1,ind}];
    coefficients = [coefficients, Acoef{1,ind}'];
    for i = 1:size(Aindex{1,ind},2)
        ch = [ch, Ach(ind)];
        loc = [loc, Aloc(ind)];
    end       
end

AtomCoef = [indexes; coefficients; ch; loc];

end

    



