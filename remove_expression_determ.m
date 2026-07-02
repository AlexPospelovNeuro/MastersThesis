% A function that removes a selected expression: no random choise, not probabilistic in any way

% INPUTS:
% - A properly structured Net2 genotype
% - An index of the expression to remove 


% OUTPUTS:
% - A properly structured Net2 genotype

function NewNet2=remove_expression_determ(Net2,exp)
if ~isfield(Net2,'Expression_patterns')
    error('There are no expressions in the genotype')
end
if exp>size(Net2.Expression_patterns,2)
     error('The target expression does not exist')
end
NewNet2=Net2;  
NewNet2.Expression_patterns(exp)=[];
%legality(NewNet2)
end