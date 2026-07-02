% The function that allocates the values of the expression pattern to the
% entity.

% Input: 
% Ent - the vector of length N, the list of entities.
% coefs - the vector of polynome coefficients (usually of powers 0,1,2; but there are no restrictions in principle)


function values=pattern_allocation(Ent,coefs)
if size(Ent,1)>size(Ent,2)
    Ent=Ent';
end
values(1:size(Ent,1))=coefs(1);
    for c=2:length(coefs)
        values=values+coefs(c)*(Ent.^(c-1));
    end

end