% A function that removes a selected link from a selected expression: no random choise, not probabilistic in any way

% INPUTS:
% - A properly structured Net2 genotype
% - An index of the expression to remove link from
% - An index of link to remove. This is a correct Index (2-based, as the first position is reserved for the expressing neuron)

% OUTPUTS:
% - A properly structured Net2 genotype



function NewNet2=remove_link_determ(Net2,exp,pat) 
if ~isfield(Net2,'Expression_patterns')
    error('There are no expressions in the genotype')
end
if exp>size(Net2.Expression_patterns,2)
     error('The target expression does not exist')
end
if pat>size(Net2.Expression_patterns{exp},2)-1
    error('The target pattern does not exist')
end
NewNet2=Net2;  

NewNet2.Expression_patterns{exp}(pat+1)=[];

% if size(NewNet2.Expression_patterns{exp},2)==1 % If no more links left in the expression pattern
%     NewNet2.Expression_patterns(exp)=[];
% end


%legality(Net2)
end
