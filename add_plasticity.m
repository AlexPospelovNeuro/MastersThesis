% A function to create a new plasticity pattern

% INPUTS:
% - A properly structured Net2 genotype
% - A probability p of the plasticity creation
% OUTPUTS:
% - A properly structured Net2 genotype


function [NewNet2,transition_tracker]=add_plasticity(Net2,p,transition_tracker)
NewNet2=Net2;

if rand(1)<p
    Plast_types=[1 2 3 1000];
    if rand<p
        Plast_types=[1 2 3 1000];
        if isfield(NewNet2,'Plasticity')
            Index=size(NewNet2.Plasticity,2)+1;
        else
            Index=1;
        end
        type=Plast_types(randperm(size(Plast_types,2),1)); % Choosing the random type for the plasticity out of the existing set;
        NewNet2.Plasticity{Index}{1}=type;
        if (type==1)||(type==1000)
            NewNet2.Plasticity{Index}{2}=NaN; % No observation frame for these types of the plasticity
        else
            NewNet2.Plasticity{Index}{2}(1,1)=randi([1 10]); % Beginning of the observation frame
            NewNet2.Plasticity{Index}{2}(1,2)=randi([NewNet2.Plasticity{Index}{2}(1,1) 20]); % End of the observation frame
            NewNet2.Plasticity{Index}{2}(2,1)=abs(randn); % SD for the beginning
            NewNet2.Plasticity{Index}{2}(2,2)=abs(randn); % SD for the end
        end
        NewNet2.Plasticity{Index}{3}(1)=randi([size(NewNet2.Cells.Input,2)+1 size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2)+size(NewNet2.Cells.Ion,2)+size(NewNet2.Cells.Mod,2)]); % Postsynapse, can be any class except the input
        if type ~= 1000 % if it is not relative refracterity
            NewNet2.Plasticity{Index}{3}(2)=randi([1 size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2)+size(NewNet2.Cells.Ion,2)]); % The presynapse for the plasticity cannot be modulatory 
        end
        NewNet2.Plasticity{Index}{4}{1}(1,1)=abs(randn)*0.1; % Mean single trigger effect
        NewNet2.Plasticity{Index}{4}{1}(2,1)=abs(randn)*0.01; % SD single trigger effect
        NewNet2.Plasticity{Index}{4}{2}(1,1)=randi([1 10]); % Mean max triggers number
        NewNet2.Plasticity{Index}{4}{2}(2,1)=abs(randn)*1; % SD max triggers number
        NewNet2.Plasticity{Index}{4}{3}(1,1)=randi([1 10]); % Mean delay
        NewNet2.Plasticity{Index}{4}{3}(2,1)=abs(randn)*1; % SD delay
        NewNet2.Plasticity{Index}{4}{3}(1,2)=randi([1 10]); % Mean in
        NewNet2.Plasticity{Index}{4}{3}(2,2)=abs(randn)*1; % SD in
        NewNet2.Plasticity{Index}{4}{3}(1,3)=randi([1 10]); % Mean plateau
        NewNet2.Plasticity{Index}{4}{3}(2,3)=abs(randn)*1; % SD plateau
        NewNet2.Plasticity{Index}{4}{3}(1,4)=randi([1 10]); % Mean out
        NewNet2.Plasticity{Index}{4}{3}(2,4)=abs(randn)*1; % SD out
    end
    if isstruct(transition_tracker)
        if isempty(transition_tracker.plasticities)
            transition_tracker.plasticities(1,1)=NaN;
            transition_tracker.plasticities(2,1)=1;
        else
            transition_tracker.plasticities(1,end+1)=NaN;
            transition_tracker.plasticities(2,end)=max(transition_tracker.plasticities(2,:))+1;
        end
    end
end
%legality(NewNet2)
end