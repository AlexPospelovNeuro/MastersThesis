% Function that adds a link (pattern) to an existing expression of the Net2 genotype

% INPUTS:
% - A properly structured Net2 genotype
% - An index of the expression to add link to
% - A probability p (with which the link will be added )

% OUTPUTS:
% - A properly structured Net2 genotype

function [NewNet2,flag]=add_link(Net2,exp,p)
flag=0;
%Failflag=0;
if ~isfield(Net2,'Expression_patterns')
    error('There are no expressions in the genotype')
end
if exp>size(Net2.Expression_patterns,2)
     error('The target expression does not exist')
end
NewNet2=Net2;  
if rand(1)<p
     flag=1;
     new_pat=size(NewNet2.Expression_patterns{exp},2)+1;
            Nonmod_classes=1:size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2)+size(NewNet2.Cells.Ion,2);
            Mod_classes=size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2)+size(NewNet2.Cells.Ion,2)+1:size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2)+size(NewNet2.Cells.Ion,2)+size(NewNet2.Cells.Mod,2);
            Features_set=[1 2 3 4 5 6 7 8 9 10 11];
            Feature=Features_set(randperm(size(Features_set,2),1));
           
            

            if Feature==1
                NewNet2.Expression_patterns{exp}{new_pat}{1}(1)=Feature;
                if ismember(NewNet2.Expression_patterns{exp}{1},Nonmod_classes) % If the expressing subclass is not modulatory 
                    item=randi([1,4]); % Then the effect on the connectionc can be both presynaptic and postsynaptic
                else % If the expressing subclass is modulatory
                    item_set=[1,3]; % It can affect only postsynaptic side of the connections
                    item=item_set(randperm(size(item_set,2),1));
                end
                NewNet2.Expression_patterns{exp}{new_pat}{1}(2)=item;
                if ismember(item,[1 3]) % If the connection pattern affects postsynaptic side of the synapse
                    target_class=randi([1, size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2)+size(NewNet2.Cells.Ion,2)]); % any nonmodulatory class can be the source
                elseif ismember(item,[2 4]) % If the connection affect presynaptic side
                    target_class=randi([size(NewNet2.Cells.Input,2)+1, size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2)+size(NewNet2.Cells.Ion,2)+size(NewNet2.Cells.Mod,2)]); % Any noninput subclass can be its target
                else
                    error('Unknown item for a newly generated connection-related expression pattern')
                end
                NewNet2.Expression_patterns{exp}{new_pat}{1}(3)=target_class;
                NewNet2.Expression_patterns{exp}{new_pat}{2}=randn(1,3);   
            elseif Feature==2
                NewNet2.Expression_patterns{exp}{new_pat}{1}(1)=Feature;
                if ismember(NewNet2.Expression_patterns{exp}{1},Nonmod_classes) % If the expressing subclass is not modulatory 
                    item=randi([1,2]); % It can affect the basic powers from both sides
                else
                    item=1; % modulatory can only have postsynaptic effect
                end
                NewNet2.Expression_patterns{exp}{new_pat}{1}(2)=item;
                if item==1
                    target_class=randi([1, size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2)+size(NewNet2.Cells.Ion,2)]); % any nonmodulatory class can be the source
                elseif item==2
                    target_class=randi([size(NewNet2.Cells.Input,2)+1, size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2)+size(NewNet2.Cells.Ion,2)+size(NewNet2.Cells.Mod,2)]); % Any noninput subclass can be its target
                else
                    error('Unknown item for a newly generated basic powers-related expression pattern')
                end
                NewNet2.Expression_patterns{exp}{new_pat}{1}(3)=target_class;
                NewNet2.Expression_patterns{exp}{new_pat}{2}=randn(1,3);
            elseif Feature==3
                NewNet2.Expression_patterns{exp}{new_pat}{1}(1)=Feature;
                if ismember(NewNet2.Expression_patterns{exp}{1},Nonmod_classes) % If the expressing subclass is not modulatory 
                    item=randi([1,2]); % It can affect the delays from both sides
                else
                    item=1; % modulatory can only have postsynaptic effect
                end
                NewNet2.Expression_patterns{exp}{new_pat}{1}(2)=item;
                if item==1
                    target_class=randi([1, size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2)+size(NewNet2.Cells.Ion,2)]); % any nonmodulatory class can be the source
                elseif item==2
                    target_class=randi([size(NewNet2.Cells.Input,2)+1, size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2)+size(NewNet2.Cells.Ion,2)+size(NewNet2.Cells.Mod,2)]); % Any noninput subclass can be its target
                else
                    error('Unknown item for a newly generated delays-related expression pattern')
                end
                NewNet2.Expression_patterns{exp}{new_pat}{1}(3)=target_class;
                NewNet2.Expression_patterns{exp}{new_pat}{2}=randn(1,3);
            elseif Feature==4
                NewNet2.Expression_patterns{exp}{new_pat}{1}(1)=Feature;
                NewNet2.Expression_patterns{exp}{new_pat}{2}=randn(1,3);
                % Here we have to do nothing, it works just as it is
            elseif Feature==5
                NewNet2.Expression_patterns{exp}{new_pat}{1}(1)=Feature;
                NewNet2.Expression_patterns{exp}{new_pat}{2}=randn(1,3);
                % Here we have to do nothing, it works just as it is
            elseif Feature==10
                NewNet2.Expression_patterns{exp}{new_pat}{1}(1)=Feature;
                NewNet2.Expression_patterns{exp}{new_pat}{2}=randn(1,3);    
            elseif Feature==11
                NewNet2.Expression_patterns{exp}{new_pat}{1}(1)=Feature;
                NewNet2.Expression_patterns{exp}{new_pat}{2}=randn(1,3);     
            elseif (Feature==6)&&(isfield(Net2,'Plasticity')) % If the feature is a plasticity feature and there is plasticity at all in the genotype
                Possible_plast_targets=[];
                for plast=1:size(NewNet2.Plasticity,2) % check if there is plasticity in the expressing type of the neurons, list the options
                    if (NewNet2.Plasticity{plast}{1}(1)==1000) && (NewNet2.Plasticity{plast}{3}(1)==NewNet2.Expression_patterns{exp}{1}) % If it is the relative refracterity and the expressing neuron fit
                        Possible_plast_targets=[Possible_plast_targets' [plast NewNet2.Plasticity{plast}{1}(1) 1]']';
                    end
                    if (ismember(NewNet2.Plasticity{plast}{1}(1),[1 2 3])) && (NewNet2.Plasticity{plast}{3}(1)==NewNet2.Expression_patterns{exp}{1}) % If it is not the relative refracterity and the expressing neuron fit and is postsynaptic
                        Possible_plast_targets=[Possible_plast_targets' [plast NewNet2.Plasticity{plast}{1}(1) 1]']';
                    end
                    if (ismember(NewNet2.Plasticity{plast}{1}(1),[1 2 3])) && (NewNet2.Plasticity{plast}{3}(2)==NewNet2.Expression_patterns{exp}{1}) % If it is not the relative refracterity and the expressing neuron fit and is presynaptic
                        if ismember(NewNet2.Expression_patterns{exp}{1},Nonmod_classes) % If the expressing subclass is not modulatory: modulatory subclass cannot be a source of the plasticity effect
                            Possible_plast_targets=[Possible_plast_targets' [plast NewNet2.Plasticity{plast}{1}(1) 2]']';
                        end
                    end
                end
                if ~isempty(Possible_plast_targets) %If there are possible plasticity targets, assign the link, otherwise nothing happens
                    NewNet2.Expression_patterns{exp}{new_pat}{1}(1)=Feature;
                    if size(Possible_plast_targets,1)>0 % if there are any valid plasticity targets
                        plast_target=randi([1 size(Possible_plast_targets,1)]); % The actual target pattern
                        if ismember(Possible_plast_targets(plast_target,2),[1 1000]) % If the target plastcity type does not have the observation frame
                            plast_featurelist=1:6;
                        elseif ismember(Possible_plast_targets(plast_target,2),[2 3]) % If the target plastcity type does have the observation frame
                            plast_featurelist=1:8;
                        else
                            error('Unknown plasticity type is targeted by the expression pattern')
                        end
                        target_feature=plast_featurelist(randperm(size(plast_featurelist,2),1));
                        NewNet2.Expression_patterns{exp}{new_pat}{1}(2)=target_feature;
                        NewNet2.Expression_patterns{exp}{new_pat}{1}(3)=Possible_plast_targets(plast_target,3);
                        NewNet2.Expression_patterns{exp}{new_pat}{1}(4)=Possible_plast_targets(plast_target,1);
                    end
                    NewNet2.Expression_patterns{exp}{new_pat}{2}=randn(1,3);
%                 else
%                     Failflag=1;
                end
            elseif Feature==7 % Since only the postsynaptic expression can affect the PSP shape pattern, no check for the modulatory/nonmodulatory issue
                NewNet2.Expression_patterns{exp}{new_pat}{1}(1)=Feature;
                NewNet2.Expression_patterns{exp}{new_pat}{1}(2)=randi([1,3]);
                NewNet2.Expression_patterns{exp}{new_pat}{1}(3)=randi([1,size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2)+size(NewNet2.Cells.Ion,2)]);
                NewNet2.Expression_patterns{exp}{new_pat}{2}=randn(1,3);
            elseif Feature==8
                NewNet2.Expression_patterns{exp}{new_pat}{1}(1)=Feature;
                if ismember(NewNet2.Expression_patterns{exp}{1},Nonmod_classes) % If the expressing subclass is not modulatory, it can be either of the members of the targeted pair 
                    item_set=[1,3]; % Since the expressing neuron is not modulatory, it can affect the modulatory connections only on the postsynaptic side
                    item=item_set(randperm(size(item_set,2),1));
                else
                    item_set=[1,2,3,4]; % The modulatory neuron can affect the feature as both presynapse and postsynapse
                    item=item_set(randperm(size(item_set,2),1));
                end
                NewNet2.Expression_patterns{exp}{new_pat}{1}(2)=item;
                if ismember(NewNet2.Expression_patterns{exp}{1},Nonmod_classes) % If the expressing subclass is not modulatory,
                    if rand(1)<0.5 % With p=0.5, the expressing neuron is presynaptic in the modulated pair
                        Postsyn=randi([size(NewNet2.Cells.Input,2)+1, size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2)+size(NewNet2.Cells.Ion,2)+size(NewNet2.Cells.Mod,2)]); % Then the postsynapse can be any noninput subclass
                        NewNet2.Expression_patterns{exp}{new_pat}{1}(3)=Postsyn;
                        NewNet2.Expression_patterns{exp}{new_pat}{1}(4)=NewNet2.Expression_patterns{exp}{1}; % Presynapse is the expressing subclass
                    else % The expressing neuron is postsynaptic in the modulated pair
                        Presyn=randi([1, size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2)+size(NewNet2.Cells.Ion,2)]); % Then the presynapse has to be nonmodulatory
                        NewNet2.Expression_patterns{exp}{new_pat}{1}(3)=NewNet2.Expression_patterns{exp}{1}; % Postsynapse is the expressing subclass
                        NewNet2.Expression_patterns{exp}{new_pat}{1}(4)=Presyn;
                    end
                    NewNet2.Expression_patterns{exp}{new_pat}{1}(5)=Mod_classes(randperm(size(Mod_classes,2),1)); % Add the modulatory input, may be any modulatory subclass
                else % If the expressing neuron is modulatory, it can be either a postsynapse in the modulated pair or the modulatory input
                    if ismember(NewNet2.Expression_patterns{exp}{new_pat}{1}(2),[2 4]) % If the link affects the presynaptic feature, it is certainly a modulatory input,
                        Postsyn=randi([size(NewNet2.Cells.Input,2)+1, size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2)+size(NewNet2.Cells.Ion,2)+size(NewNet2.Cells.Mod,2)]); % Then the postsynapse can be any noninput subclass
                        NewNet2.Expression_patterns{exp}{new_pat}{1}(3)=Postsyn;
                        Presyn=randi([1, size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2)+size(NewNet2.Cells.Ion,2)]); % Then the presynapse has to be nonmodulatory
                        NewNet2.Expression_patterns{exp}{new_pat}{1}(4)=Presyn;
                        NewNet2.Expression_patterns{exp}{new_pat}{1}(5)=NewNet2.Expression_patterns{exp}{1}; % The modulatory input of the triplet has to be the expressing subclass
                    else % If the link affects the postsynaptic feature yet the expressing subclass is modulatory, it has to be a postsynapse in the modulated pair
                        NewNet2.Expression_patterns{exp}{new_pat}{1}(3)=NewNet2.Expression_patterns{exp}{1}; % The postsynapse is the expressing subclass
                        Presyn=randi([1, size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2)+size(NewNet2.Cells.Ion,2)]); % Then the presynapse has to be nonmodulatory
                        NewNet2.Expression_patterns{exp}{new_pat}{1}(4)=Presyn;
                        NewNet2.Expression_patterns{exp}{new_pat}{1}(5)=Mod_classes(randperm(size(Mod_classes,2),1)); % Add the modulatory input, may be any modulatory subclass
                    end
                end    
                NewNet2.Expression_patterns{exp}{new_pat}{2}=randn(1,3);
            elseif Feature==9
                NewNet2.Expression_patterns{exp}{new_pat}{1}(1)=Feature;
                item=randi([1,6]);
                NewNet2.Expression_patterns{exp}{new_pat}{1}(2)=item;
                if ismember(NewNet2.Expression_patterns{exp}{1},Nonmod_classes) % If the expressing subclass is not modulatory, it can be either of the members of the targeted pair 
                    loc_set=[1]; % Since the expressing neuron is not modulatory, it can affect the modulatory connections only on the postsynaptic side
                    loc=loc_set(randperm(size(loc_set,2),1));
                else
                    loc_set=[1,2]; % The modulatory neuron can affect the feature as both presynapse and postsynapse
                    loc=loc_set(randperm(size(loc_set,2),1));
                end
                if ismember(NewNet2.Expression_patterns{exp}{1},Nonmod_classes) % If the expressing subclass is not modulatory 
                    if rand(1)<0.5 % With p=0.5, the expressing neuron is presynaptic in the modulated pair
                        Postsyn=randi([size(NewNet2.Cells.Input,2)+1, size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2)+size(NewNet2.Cells.Ion,2)+size(NewNet2.Cells.Mod,2)]); % Then the postsynapse can be any noninput subclass
                        NewNet2.Expression_patterns{exp}{new_pat}{1}(3)=Postsyn;
                        NewNet2.Expression_patterns{exp}{new_pat}{1}(4)=NewNet2.Expression_patterns{exp}{1};
                    else % The expressing neuron is postsynaptic in the modulated pair
                        Presyn=randi([1, size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2)+size(NewNet2.Cells.Ion,2)]); % Then the presynapse has to be nonmodulatory
                        NewNet2.Expression_patterns{exp}{new_pat}{1}(3)=NewNet2.Expression_patterns{exp}{1};
                        NewNet2.Expression_patterns{exp}{new_pat}{1}(4)=Presyn;
                    end
                    NewNet2.Expression_patterns{exp}{new_pat}{1}(5)=Mod_classes(randperm(size(Mod_classes,2),1)); % Add the modulatory input index
                    NewNet2.Expression_patterns{exp}{new_pat}{1}(6)=1; % The effect is certainly postsynaptic since the expressing neuron is the target
                else % If the expressing neuron is modulatory, it can be either a postsynapse in the modulated pair or the modulatory input
                   if loc==2 % If the link affects the presynaptic feature, it is certainly a modulatory input,
                        Postsyn=randi([size(NewNet2.Cells.Input,2)+1, size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2)+size(NewNet2.Cells.Ion,2)+size(NewNet2.Cells.Mod,2)]); % Then the postsynapse can be any noninput subclass
                        NewNet2.Expression_patterns{exp}{new_pat}{1}(3)=Postsyn;
                        Presyn=randi([1, size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2)+size(NewNet2.Cells.Ion,2)]); % Then the presynapse has to be nonmodulatory
                        NewNet2.Expression_patterns{exp}{new_pat}{1}(4)=Presyn;
                        NewNet2.Expression_patterns{exp}{new_pat}{1}(5)=NewNet2.Expression_patterns{exp}{1}; % The modulatory input of the triplet has to be the expressing subclass
                    else % If the link affects the postsynaptic feature yet the expressing subclass is modulatory, it has to be a postsynapse in the modulated pair
                        NewNet2.Expression_patterns{exp}{new_pat}{1}(3)=NewNet2.Expression_patterns{exp}{1}; % The postsynapse is the expressing subclass
                        Presyn=randi([1, size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2)+size(NewNet2.Cells.Ion,2)]); % Then the presynapse has to be nonmodulatory
                        NewNet2.Expression_patterns{exp}{new_pat}{1}(4)=Presyn;
                        NewNet2.Expression_patterns{exp}{new_pat}{1}(5)=Mod_classes(randperm(size(Mod_classes,2),1)); % Add the modulatory input, may be any modulatory subclass
                    end
                    NewNet2.Expression_patterns{exp}{new_pat}{1}(6)=loc;
                end
                NewNet2.Expression_patterns{exp}{new_pat}{2}=randn(1,3);
            end
%             if Failflag==0
%                 NewNet2.Expression_patterns{exp}{new_pat}{2}=randn(1,3);
%             end
    
%legality(NewNet2)    
end






end