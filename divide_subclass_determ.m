% A function that divides the selected subclass 
% INPUTS:
% - A properly structured Net2 genotype
% - An index of the subclass to divide
% OUTPUTS:
% - A properly structured Net2 genotype



function [NewNet2,transition_tracker]=divide_subclass_determ(Net2,cl,transition_tracker)

% report.subclasses=[cl;cl+1];
% report.plast=[];
% report.expressions=[];
% report.links={};
%report



NewNet2=Net2;
if cl<=(size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2))
    error('division of a reserved input and output neurons is forbidden')
end
if cl>(size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2)+size(NewNet2.Cells.Ion,2)+size(NewNet2.Cells.Mod,2))
    error('division of a non-existing subclass')
end
if cl>(size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2)+size(NewNet2.Cells.Ion,2))
    type=2; % modulatory class
    Index=cl-size(NewNet2.Cells.Input,2)-size(NewNet2.Cells.Output,2)-size(NewNet2.Cells.Ion,2); % Index inside the class
else
    type=1; % ionotropic class
    Index=cl-size(NewNet2.Cells.Input,2)-size(NewNet2.Cells.Output,2); % Index inside the class
end

transition_tracker.subclasses(2,transition_tracker.subclasses(2,:)>cl)=transition_tracker.subclasses(2,transition_tracker.subclasses(2,:)>cl)+1;
transition_tracker.subclasses(1,end+1)=cl;
transition_tracker.subclasses(2,end)=cl+1;



if type==1
    N_orig=NewNet2.Cells.Ion(:,Index);
    if NewNet2.Cells.Ion(1,Index)<2
        N_child1=1;
        N_child2=1;
    else
        %N_child1=randi(NewNet2.Cells.Ion(1,Index));
        N_child1=rand(1)*(NewNet2.Cells.Ion(1,Index)-1)+1;
        N_child2=NewNet2.Cells.Ion(1,Index)-N_child1;
    end
    NewNet2.Cells.Ion=[NewNet2.Cells.Ion(:,1:Index-1) [N_child1 NewNet2.Cells.Ion(2,Index)]' [N_child2 NewNet2.Cells.Ion(2,Index)]' NewNet2.Cells.Ion(:,Index+1:end)]; % removing the cell numbers
elseif type==2
    N_orig=NewNet2.Cells.Mod(:,Index);
    if NewNet2.Cells.Mod(1,Index)<2
        N_child1=1;
        N_child2=1;
    else
        %N_child1=randi(NewNet2.Cells.Mod(1,Index));
        N_child1=rand(1)*(NewNet2.Cells.Mod(1,Index)-1)+1;
        N_child2=NewNet2.Cells.Mod(1,Index)-N_child1;
    end
    NewNet2.Cells.Mod=[NewNet2.Cells.Mod(:,1:Index-1) [N_child1 NewNet2.Cells.Mod(2,Index)]' [N_child2 NewNet2.Cells.Mod(2,Index)]' NewNet2.Cells.Mod(:,Index+1:end)]; % removing the cell numbers
else
    error('Unknown subclass to divide')
end

Fract1=N_child1/(N_child1+N_child2);
Fract2=N_child2/(N_child1+N_child2);
% Here is the trick: the Net2.Connections(:,:,1) contains the post_mu, or the average number of the input connections per postsynaptic neuron in the ordered pair of subclasses Post <- Pre. When the subclass is being
% divided, for both child subclasses, this number remains the same in any ordered pair where the divided subclass is postsynaptic. If the divided subclass is presynaptic, the children subclasses should
% "split" their outputs proportional to their "part" in the split. The modulatory subclass can be only postsynaptic in with respect to the ionotropic connection, so the "split" is required only for the ionotropic
% dividing subclass here

if type==1
    Newconnect(1:size(NewNet2.Connections,1)+1,1:size(NewNet2.Connections,2)+1,1:size(NewNet2.Connections,3))=NaN;
    Newconnect(1:cl,1:cl,:)=NewNet2.Connections(1:cl,1:cl,:); 
    Newconnect(1:cl,cl+1:end,:)=NewNet2.Connections(1:cl,cl:end,:); 
    Newconnect(cl+1:end,1:cl,:)=NewNet2.Connections(cl:end,1:cl,:);
    Newconnect(cl+1:end,cl+1:end,:)=NewNet2.Connections(cl:end,cl:end,:);
    Newconnect(:,cl,1)=Newconnect(:,cl,1)*Fract1; % The splitting of the output connections
    Newconnect(:,cl+1,1)=Newconnect(:,cl+1,1)*Fract2;
elseif type==2
    Newconnect(1:size(NewNet2.Connections,1)+1,1:size(NewNet2.Connections,2),1:size(NewNet2.Connections,3))=NaN;
    Newconnect(1:cl,:,:)=NewNet2.Connections(1:cl,:,:); 
    Newconnect(cl+1:end,:,:)=NewNet2.Connections(cl:end,:,:);
end
NewNet2.Connections=Newconnect;

% The modulatory class can be presynaptic for the modulatory input, and it is a subject for "splitting" if so. 
if type==1
    Newconnect_Mod(1:size(NewNet2.Connections_Mod,1)+1,1:size(NewNet2.Connections_Mod,2)+1,1:size(NewNet2.Connections_Mod,3),1:size(NewNet2.Connections_Mod,4))=NaN;
    Newconnect_Mod(1:cl,1:cl,:,:)=NewNet2.Connections_Mod(1:cl,1:cl,:,:); 
    Newconnect_Mod(1:cl,cl+1:end,:,:)=NewNet2.Connections_Mod(1:cl,cl:end,:,:);
    Newconnect_Mod(cl+1:end,1:cl,:,:)=NewNet2.Connections_Mod(cl:end,1:cl,:,:);
    Newconnect_Mod(cl+1:end,cl+1:end,:,:)=NewNet2.Connections_Mod(cl:end,cl:end,:,:);
elseif type==2
    Newconnect_Mod(1:size(NewNet2.Connections_Mod,1)+1,1:size(NewNet2.Connections_Mod,2),1:size(NewNet2.Connections_Mod,3)+1,1:size(NewNet2.Connections_Mod,4))=NaN;
    Newconnect_Mod(1:cl,:,1:Index,:)=NewNet2.Connections_Mod(1:cl,:,1:Index,:); 
    Newconnect_Mod(1:cl,:,Index+1:end,:)=NewNet2.Connections_Mod(1:cl,:,Index:end,:); 
    Newconnect_Mod(cl+1:end,:,1:Index,:)=NewNet2.Connections_Mod(cl:end,:,1:Index,:);
    Newconnect_Mod(cl+1:end,:,Index+1:end,:)=NewNet2.Connections_Mod(cl:end,:,Index:end,:);
    Newconnect_Mod(:,:,Index,1)=Newconnect_Mod(:,:,Index,1)*Fract1; % The splitting of the output connections
    Newconnect_Mod(:,:,Index+1,1)=Newconnect_Mod(:,:,Index+1,1)*Fract2; 
end
NewNet2.Connections_Mod=Newconnect_Mod;

if type==1
    BasicPowers(1:size(NewNet2.BasicPowers,1)+1,1:size(NewNet2.BasicPowers,2)+1,1:size(NewNet2.BasicPowers,3))=NaN;
    BasicPowers(1:cl,1:cl,:)=NewNet2.BasicPowers(1:cl,1:cl,:); 
    BasicPowers(1:cl,cl+1:end,:)=NewNet2.BasicPowers(1:cl,cl:end,:);
    BasicPowers(cl+1:end,1:cl,:)=NewNet2.BasicPowers(cl:end,1:cl,:);
    BasicPowers(cl+1:end,cl+1:end,:)=NewNet2.BasicPowers(cl:end,cl:end,:);
elseif type==2
    BasicPowers(1:size(NewNet2.BasicPowers,1)+1,1:size(NewNet2.BasicPowers,2),1:size(NewNet2.BasicPowers,3))=NaN;
    BasicPowers(1:cl,:,:)=NewNet2.BasicPowers(1:cl,:,:); 
    BasicPowers(cl+1:end,:,:)=NewNet2.BasicPowers(cl:end,:,:);
end
NewNet2.BasicPowers=BasicPowers;

if type==1
    Delays(1:size(NewNet2.Delays,1)+1,1:size(NewNet2.Delays,2)+1,1:size(NewNet2.Delays,3))=NaN;
    Delays(1:cl,1:cl,:)=NewNet2.Delays(1:cl,1:cl,:); 
    Delays(1:cl,cl+1:end,:)=NewNet2.Delays(1:cl,cl:end,:);
    Delays(cl+1:end,1:cl,:)=NewNet2.Delays(cl:end,1:cl,:);
    Delays(cl+1:end,cl+1:end,:)=NewNet2.Delays(cl:end,cl:end,:);
elseif type==2
    Delays(1:size(NewNet2.Delays,1)+1,1:size(NewNet2.Delays,2),1:size(NewNet2.Delays,3))=NaN;
    Delays(1:cl,:,:)=NewNet2.Delays(1:cl,:,:); 
    Delays(cl+1:end,:,:)=NewNet2.Delays(cl:end,:,:);
end
NewNet2.Delays=Delays;

if type==1
    PSPshape(1:size(NewNet2.PSPshape,1)+1,1:size(NewNet2.PSPshape,2)+1)={NaN};
    PSPshape(1:cl,1:cl)=NewNet2.PSPshape(1:cl,1:cl); 
    PSPshape(1:cl,cl+1:end)=NewNet2.PSPshape(1:cl,cl:end);
    PSPshape(cl+1:end,1:cl)=NewNet2.PSPshape(cl:end,1:cl);
    PSPshape(cl+1:end,cl+1:end)=NewNet2.PSPshape(cl:end,cl:end);
elseif type==2
    PSPshape(1:size(NewNet2.PSPshape,1)+1,1:size(NewNet2.PSPshape,2))={NaN};
    PSPshape(1:cl,:)=NewNet2.PSPshape(1:cl,:); 
    PSPshape(cl+1:end,:)=NewNet2.PSPshape(cl:end,:);
end
NewNet2.PSPshape=PSPshape;

NewNet2.Thresholds=[NewNet2.Thresholds(:,1:cl) NewNet2.Thresholds(:,cl:end)]; %  thresholds
NewNet2.AbsRefract=[NewNet2.AbsRefract(:,1:cl) NewNet2.AbsRefract(:,cl:end)]; %  absolute refracterity
NewNet2.ThreshNoise=[NewNet2.ThreshNoise(:,1:cl) NewNet2.ThreshNoise(:,cl:end)]; %  threshold noise
NewNet2.RecurCon=[NewNet2.RecurCon(:,1:cl) NewNet2.RecurCon(:,cl:end)]; %  recurrent connections

if type==1
    Mod(1:size(NewNet2.Mod,1)+1,1:size(NewNet2.Mod,2)+1,1:size(NewNet2.Mod,3))={NaN};
    Mod(1:cl,1:cl,:)=NewNet2.Mod(1:cl,1:cl,:); 
    Mod(1:cl,cl+1:end,:)=NewNet2.Mod(1:cl,cl:end,:);
    Mod(cl+1:end,1:cl,:)=NewNet2.Mod(cl:end,1:cl,:);
    Mod(cl+1:end,cl+1:end,:)=NewNet2.Mod(cl:end,cl:end,:);
elseif type==2
    Mod(1:size(NewNet2.Mod,1)+1,1:size(NewNet2.Mod,2),1:size(NewNet2.Mod,3)+1)={NaN};
    Mod(1:cl,:,1:Index)=NewNet2.Mod(1:cl,:,1:Index); 
    Mod(1:cl,:,Index+1:end)=NewNet2.Mod(1:cl,:,Index:end);
    Mod(cl+1:end,:,1:Index)=NewNet2.Mod(cl:end,:,1:Index);
    Mod(cl+1:end,:,Index+1:end)=NewNet2.Mod(cl:end,:,Index:end);
end
NewNet2.Mod=Mod;

if isfield(NewNet2,'Plasticity') % If there are plasticity patterns in the genotype
    
    Newplast=[]; % A variable to keep track of new plasticities created by the class division (for the sake of expression pattern references)
    N_plast=size(NewNet2.Plasticity,2);
    for plast=1:N_plast % Moving from low to high index this time, adding new plasticity patterns to the end
        affected=NewNet2.Plasticity{plast}{3}; % The indices of the classes that are related to the plasticity with index "plast"
        if size(affected,2)==1 % Plasticity is defined for a single subclass (type 1000 in the current version of the model)
            if affected==cl % If the plasticity is defined for the dividing subclass
                NextIndex=size(NewNet2.Plasticity,2)+1; % There will be a new plasticity with the new index, appended to the list of plasticities.
                NewNet2.Plasticity(NextIndex)=NewNet2.Plasticity(plast); % The new plasticity features are copy of the old one
                NewNet2.Plasticity{NextIndex}{3}(1)=NewNet2.Plasticity{NextIndex}{3}(1)+1; % A reference to the new subclass in the new plasticity.
                Newplast=[Newplast [plast; NextIndex]]; % Keep the indices of the dividing plastity in a from of pair of indices: [plast] (what was divided) and [NextIndex] (Index of the new plasicity)
            else
                if affected>cl 
                    NewNet2.Plasticity{plast}{3}(1)=NewNet2.Plasticity{plast}{3}(1)+1; % Increase index of the subclass reference in the plasticities if the Index is bigger than of the class that is being divided
                    
                end % If the index of the subclass that has the current plasticity is less than the index of the divided subclass, no changes have to be done. 
            end
        elseif size(affected,2)==2 % Plasticity is defined for a pair of subclasses
            if affected(1)>cl  % if the index of the subclass that has the plasticity is bigger than the index of one subclass being divided
                NewNet2.Plasticity{plast}{3}(1)=NewNet2.Plasticity{plast}{3}(1)+1; % If it is the first one in the pair, increase the index of the first one in the pair
            end
            if affected(2)>cl
                NewNet2.Plasticity{plast}{3}(2)=NewNet2.Plasticity{plast}{3}(2)+1; % If it is the second one in the pair, increase the second one in the pair
            end
            % If the index of one (or both) subclasses that form the pair is less than the index of the divided sublcass, do nothing. These actions are done independently for both indices in the pair. 
            if (affected(1)==cl)&&(affected(2)~=cl) % If the postsynaptic subclass was divided (and the presynaptic was not)
                NextIndex=size(NewNet2.Plasticity,2)+1; % An index of a new plasticity, to append in the end of the list
                NewNet2.Plasticity(NextIndex)=NewNet2.Plasticity(plast); % Copy the features of the plasticity, append
                NewNet2.Plasticity{NextIndex}{3}(1)=NewNet2.Plasticity{NextIndex}{3}(1)+1; % A reference to the new postsynaptic subclass. The index of the second subclass in the pair was increased previously if required. 
                Newplast=[Newplast [plast; NextIndex]]; 
            elseif (affected(1)~=cl)&&(affected(2)==cl) % If the presynaptic subclass was divided (and the postsynaptic was not)
                NextIndex=size(NewNet2.Plasticity,2)+1; % An index of a new plasticity, to append in the end of the list
                NewNet2.Plasticity(NextIndex)=NewNet2.Plasticity(plast); % Copy and append
                NewNet2.Plasticity{NextIndex}{3}(2)=NewNet2.Plasticity{NextIndex}{3}(2)+1; % A reference to the new presynaptic subclass. The index of the second subclass in the pair was increased previously if required. 
                Newplast=[Newplast [plast; NextIndex]];
            elseif (affected(1)==cl)&&(affected(2)==cl) % If presynaptic and postsynaptic subclass is the same and was divided. So, original plasticity is in [cl cl] pair of subclasses
                NextIndex=size(NewNet2.Plasticity,2)+1; % An index of a new plasticity, to append in the end of the list
                NewNet2.Plasticity(NextIndex)=NewNet2.Plasticity(plast); % Copy and append for the first time
                NewNet2.Plasticity{NextIndex}{3}(1)=NewNet2.Plasticity{NextIndex}{3}(1)+1; % A reference to the new postsynaptic subclass. First new plasticity is in [cl+1 cl] pair of subclasses
                NewNet2.Plasticity(NextIndex+1)=NewNet2.Plasticity(plast); % Copy and append for the second time
                NewNet2.Plasticity{NextIndex+1}{3}(2)=NewNet2.Plasticity{NextIndex+1}{3}(2)+1; % A reference to the new presynaptic subclass. Second new plasticity is in [cl cl+1] pair of subclasses
                NewNet2.Plasticity(NextIndex+2)=NewNet2.Plasticity(plast); % Copy and append for the third time
                NewNet2.Plasticity{NextIndex+2}{3}(1)=NewNet2.Plasticity{NextIndex+2}{3}(1)+1; % A reference to the new postsynaptic subclass
                NewNet2.Plasticity{NextIndex+2}{3}(2)=NewNet2.Plasticity{NextIndex+2}{3}(2)+1; % A reference to the new presynaptic subclass. Third new plasticity is in [cl+1 cl+1] pair of subclasses
                Newplast=[Newplast [plast; NextIndex] [plast; NextIndex+1] [plast; NextIndex+2]];
            end
        end
    end
    %report.plast=Newplast;
    for update=1:size(Newplast,2)
        if isstruct(transition_tracker)
            transition_tracker.plasticities(1,end+1)=transition_tracker.plasticities(1,transition_tracker.plasticities(2,:)==Newplast(1,update)); % The new item in the tracker has some "origin" in the original indexing, that now has a new index similarequal to the "old" plasticity in the Newplast
            transition_tracker.plasticities(2,end)=Newplast(2,update);   % The "current" index of the new item is equal to the "new" index from the Newplast
        end
    end
    
    %transition_tracker.plasticities(2,transition_tracker.plasticities(2,:)==Newplast(1,update))
end



% Routine for the expression division


if isfield(NewNet2,'Expression_patterns') % If there are any expressions
    %disp('Expressions division beginning')
    %plast_exp_summary(NewNet2)
    N_exp=size(NewNet2.Expression_patterns,2); % The initial number of the expression
    % STAGE ONE: division of the expressions in the divided subclass, modification of the indices in the expressions
    for exp=1:N_exp % for every expression
        if NewNet2.Expression_patterns{exp}{1}(1)>cl % If expression is in the subclass with the index higher then the divided one
            NewNet2.Expression_patterns{exp}{1}(1)=NewNet2.Expression_patterns{exp}{1}(1)+1; % Increase the index
        elseif NewNet2.Expression_patterns{exp}{1}(1)==cl % If expression is in the subclass which is being divided
            NextIndex=size(NewNet2.Expression_patterns,2)+1; % The index for the new expression 
            NewNet2.Expression_patterns(NextIndex)=NewNet2.Expression_patterns(exp); % Make a copy of the expression 
            NewNet2.Expression_patterns{NextIndex}{1}(1)=NewNet2.Expression_patterns{NextIndex}{1}(1)+1; % Assign the index of the new subclass to the new expression pattern
%            report.expressions=[report.expressions [exp; NextIndex]];
            if isstruct(transition_tracker)
                transition_tracker.expressions(1,end+1)=transition_tracker.expressions(1,transition_tracker.expressions(2,:)==exp);
                transition_tracker.expressions(2,end)=NextIndex;
                transition_tracker.links{end+1}=transition_tracker.links{transition_tracker.expressions(2,:)==exp}; % add new links branch since we have a new expression 
            end
            
        end % If the index of the expressing subclass is less than the index of the divided subclass, nothing should be done yet.
    end
    %disp('Expressions division end')
    %plast_exp_summary(NewNet2)
    % END OF STAGE 1. The expressions are divided and assigned to the correct subclasses, but the links within them are not. The division and modification of the links depend on types of the links. 
    % STAGE TWO: division and modification of the links
    N_exp=size(NewNet2.Expression_patterns,2); % The new number of the expressions
    %disp('link division beginning')
    expressions_counter=0; % the expressions
    for exp=1:N_exp % for every expression
        %link_duplication_counter=0; % For the report
        N_link=size(NewNet2.Expression_patterns{exp},2);
        for pat=2:N_link % for every pattern related to the expression
            if ismember(NewNet2.Expression_patterns{exp}{pat}{1}(1),[1 2 3 7]) % If the referred subclass of the link is in the third position, there is no presynaptic/postsynaptic specifics to the downstream process
                if NewNet2.Expression_patterns{exp}{pat}{1}(3)==cl % and it is the class that is being divided
                    NextIndex=size(NewNet2.Expression_patterns{exp},2)+1; % Reserve a new index for the new pattern
                    NewNet2.Expression_patterns{exp}(NextIndex)=NewNet2.Expression_patterns{exp}(pat); % Make a new pattern as the copy of the existing one
                    NewNet2.Expression_patterns{exp}{NextIndex}{1}(3)=NewNet2.Expression_patterns{exp}{NextIndex}{1}(3)+1; % Change the index of the targer class to the new one. 
                    
                    if isstruct(transition_tracker)
                        transition_tracker.links{transition_tracker.expressions(2,:)==exp}(1,end+1)=transition_tracker.links{transition_tracker.expressions(2,:)==exp}(1,transition_tracker.links{transition_tracker.expressions(2,:)==exp}(2,:)==pat); % getting the "original" duplicating link index and setting it as the orginal for the "new" link
                        transition_tracker.links{transition_tracker.expressions(2,:)==exp}(2,end)=NextIndex; % The new link index
                    end
%                     link_duplication_counter=link_duplication_counter+1;
%                     if link_duplication_counter==1
%                         expressions_counter=expressions_counter+1;
%                         report.links{expressions_counter}{1}=exp;
%                         report.links{expressions_counter}{2}=[pat; NextIndex];
%                     else
%                         report.links{expressions_counter}{2}=[report.links{expressions_counter}{2} [pat; NextIndex]];
%                     end
                elseif NewNet2.Expression_patterns{exp}{pat}{1}(3)>cl % If the referred class has a bigger index than the one that is being divided
                    NewNet2.Expression_patterns{exp}{pat}{1}(3)=NewNet2.Expression_patterns{exp}{pat}{1}(3)+1; % Increase it by one because  of the indexation shift
                end % If the link refers to a subclass with index less than the index of the divided one, do nothing
            elseif ismember(NewNet2.Expression_patterns{exp}{pat}{1}(1),[8 9]) % If the link is "modulatory" and has triplet of subclasses inside the addresation
                % A workflow note for this part. cl is the divided subclass, it is constant for the function call. After the division, the "children" have indices cl ("original") and cl+1 ("copy"). At the stage 1 of the expression duplication, the expressions themselves were duplicated, and the
                % sunclass indexation of the expressing subclass (NewNet2.Expression_patterns{exp}{1}) was updated. But the subclass within the links was not yet, it is done in parallel with links' duplication when necessary. Therefore, inside the links, all the class indices are "old"
                
                % Mandatory indexation shifts;
                if NewNet2.Expression_patterns{exp}{pat}{1}(3)>cl % If the postsynaptic neuron of the "target" pair has an index higher than the divided subclass
                    NewNet2.Expression_patterns{exp}{pat}{1}(3)=NewNet2.Expression_patterns{exp}{pat}{1}(3)+1; % Increase that index by 1
                end
                if NewNet2.Expression_patterns{exp}{pat}{1}(4)>cl % If the presynaptic neuron of the "target" pair has an index higher than the divided subclass
                    NewNet2.Expression_patterns{exp}{pat}{1}(4)=NewNet2.Expression_patterns{exp}{pat}{1}(4)+1; % Increase that index by 1
                end
                if NewNet2.Expression_patterns{exp}{pat}{1}(5)>cl % If the modulatory input of the "target" pair has an index higher than the divided subclass
                    NewNet2.Expression_patterns{exp}{pat}{1}(5)=NewNet2.Expression_patterns{exp}{pat}{1}(5)+1; % Increase that index by 1
                end
                
                % A workflow note: at this moment (specifically for the links of types 8 and 9), the address "cl+1" should not exist: the "old" subclass cl+1 was already shifted to be cl+2, and the new cl+1 has not yet been assigned. But cl+1 exists in NewNet2.Expression_patterns{exp}{1} (the
                % expressing subclass), and it is a "copy" addressation
                
                % The expression duplication routine
                if (NewNet2.Expression_patterns{exp}{pat}{1}(3)~=cl)&&(NewNet2.Expression_patterns{exp}{pat}{1}(4)~=cl)&&(NewNet2.Expression_patterns{exp}{pat}{1}(5)~=cl)   % If none of the classes of the triplet was duplicated
                    % do nothing
                elseif (NewNet2.Expression_patterns{exp}{pat}{1}(3)==cl)&&(NewNet2.Expression_patterns{exp}{pat}{1}(4)~=cl)&&(NewNet2.Expression_patterns{exp}{pat}{1}(5)~=cl) % If the postsynaptic subclass of the pair is being divided
                    if ismember(NewNet2.Expression_patterns{exp}{1},[NewNet2.Expression_patterns{exp}{pat}{1}(3) NewNet2.Expression_patterns{exp}{pat}{1}(3)+1]) % If the dividing class is the one that expresses the pattern
                        if NewNet2.Expression_patterns{exp}{1}==cl+1 % If we are operating in the "copy" of the divided subclass
                            NewNet2.Expression_patterns{exp}{pat}{1}(3)=cl+1; % Edit the triplet within the link to match the new expressing subclass
                        elseif NewNet2.Expression_patterns{exp}{1}==cl % If we are operating in the "original" of the divided subclass, do nothing
                        else error(['Impossible situation during the expression duplication, the situation is: cl=' num2str(cl) ' , e_cl=' num2str(NewNet2.Expression_patterns{exp}{1}) ', x1=' num2str(NewNet2.Expression_patterns{exp}{pat}{1}(3)) ', y1=' num2str(NewNet2.Expression_patterns{exp}{pat}{1}(4)) ', z1=' num2str(NewNet2.Expression_patterns{exp}{pat}{1}(5))])
                        end
                    elseif (ismember(NewNet2.Expression_patterns{exp}{1},[NewNet2.Expression_patterns{exp}{pat}{1}(4) NewNet2.Expression_patterns{exp}{pat}{1}(4)+1]))||(ismember(NewNet2.Expression_patterns{exp}{1},[NewNet2.Expression_patterns{exp}{pat}{1}(5) NewNet2.Expression_patterns{exp}{pat}{1}(5)+1])) % The expressing subclass is presynapse or modulatory input
                        % The expressing class was not divided, therefore no pattern duplication; but the link duplicates
                        NextIndex=size(NewNet2.Expression_patterns{exp},2)+1; % Get the index of last pattern of the expression to append
                        NewNet2.Expression_patterns{exp}(NextIndex)=NewNet2.Expression_patterns{exp}(pat); % Append the expression with the new pattern (link)
                        NewNet2.Expression_patterns{exp}{NextIndex}{1}(3)=NewNet2.Expression_patterns{exp}{NextIndex}{1}(3)+1; % Modify the copy of the link to fit the copy of the subclass
                        if isstruct(transition_tracker)
                            transition_tracker.links{transition_tracker.expressions(2,:)==exp}(1,end+1)=transition_tracker.links{transition_tracker.expressions(2,:)==exp}(1,transition_tracker.links{transition_tracker.expressions(2,:)==exp}(2,:)==pat); % getting the "original" duplicating link index and setting it as the orginal for the "new" link
                            transition_tracker.links{transition_tracker.expressions(2,:)==exp}(2,end)=NextIndex; % The new link index
                        end
%                         link_duplication_counter=link_duplication_counter+1;
%                         if link_duplication_counter==1
%                             expressions_counter=expressions_counter+1;
%                             report.links{expressions_counter}{1}=exp;
%                             report.links{expressions_counter}{2}=[pat; NextIndex];
%                         else
%                             report.links{expressions_counter}{2}=[report.links{expressions_counter}{2} [pat; NextIndex]];
%                         end
                    else error(['Impossible situation during the expression duplication, the situation is: cl=' num2str(cl) ' , e_cl=' num2str(NewNet2.Expression_patterns{exp}{1}) ', x1=' num2str(NewNet2.Expression_patterns{exp}{pat}{1}(3)) ', y1=' num2str(NewNet2.Expression_patterns{exp}{pat}{1}(4)) ', z1=' num2str(NewNet2.Expression_patterns{exp}{pat}{1}(5))])
                    end
                
                elseif (NewNet2.Expression_patterns{exp}{pat}{1}(3)~=cl)&&(NewNet2.Expression_patterns{exp}{pat}{1}(4)==cl)&&(NewNet2.Expression_patterns{exp}{pat}{1}(5)~=cl) % If the presynaptic subclass of the pair is being divided
                    if ismember(NewNet2.Expression_patterns{exp}{1},[NewNet2.Expression_patterns{exp}{pat}{1}(4) NewNet2.Expression_patterns{exp}{pat}{1}(4)+1]) % If the dividing class is the one that expresses the pattern
                        if NewNet2.Expression_patterns{exp}{1}==cl+1 % If we are operating in the "copy" of the divided subclass
                            NewNet2.Expression_patterns{exp}{pat}{1}(4)=cl+1; % Edit the triplet within the link to match the new expressing subclass
                        elseif NewNet2.Expression_patterns{exp}{1}==cl % If we are operating in the "original" of the divided subclass, do nothing
                        else error(['Impossible situation during the expression duplication, the situation is: cl=' num2str(cl) ' , e_cl=' num2str(NewNet2.Expression_patterns{exp}{1}) ', x1=' num2str(NewNet2.Expression_patterns{exp}{pat}{1}(3)) ', y1=' num2str(NewNet2.Expression_patterns{exp}{pat}{1}(4)) ', z1=' num2str(NewNet2.Expression_patterns{exp}{pat}{1}(5))])
                        end
                    elseif (ismember(NewNet2.Expression_patterns{exp}{1},[NewNet2.Expression_patterns{exp}{pat}{1}(3) NewNet2.Expression_patterns{exp}{pat}{1}(3)+1]))||(ismember(NewNet2.Expression_patterns{exp}{1},[NewNet2.Expression_patterns{exp}{pat}{1}(5) NewNet2.Expression_patterns{exp}{pat}{1}(5)+1])) % The expressing subclass is postsynapse or modulatory input
                        % The expressing class was not divided, therefore no pattern duplication; but the link duplicates
                        NextIndex=size(NewNet2.Expression_patterns{exp},2)+1; % Get the index of last pattern of the expression to append
                        NewNet2.Expression_patterns{exp}(NextIndex)=NewNet2.Expression_patterns{exp}(pat); % Append the expression with the new pattern (link)
                        NewNet2.Expression_patterns{exp}{NextIndex}{1}(4)=NewNet2.Expression_patterns{exp}{NextIndex}{1}(4)+1; % Modify the copy of the link to fit the copy of the subclass
                        if isstruct(transition_tracker)
                            transition_tracker.links{transition_tracker.expressions(2,:)==exp}(1,end+1)=transition_tracker.links{transition_tracker.expressions(2,:)==exp}(1,transition_tracker.links{transition_tracker.expressions(2,:)==exp}(2,:)==pat); % getting the "original" duplicating link index and setting it as the orginal for the "new" link
                            transition_tracker.links{transition_tracker.expressions(2,:)==exp}(2,end)=NextIndex; % The new link index
                        end
%                         link_duplication_counter=link_duplication_counter+1;
%                         if link_duplication_counter==1
%                             expressions_counter=expressions_counter+1;
%                             report.links{expressions_counter}{1}=exp;
%                             report.links{expressions_counter}{2}=[pat; NextIndex];
%                         else
%                             report.links{expressions_counter}{2}=[report.links{expressions_counter}{2} [pat; NextIndex]];
%                         end
                    else error(['Impossible situation during the expression duplication, the situation is: cl=' num2str(cl) ' , e_cl=' num2str(NewNet2.Expression_patterns{exp}{1}) ', x1=' num2str(NewNet2.Expression_patterns{exp}{pat}{1}(3)) ', y1=' num2str(NewNet2.Expression_patterns{exp}{pat}{1}(4)) ', z1=' num2str(NewNet2.Expression_patterns{exp}{pat}{1}(5))])
                    end
                    
                elseif (NewNet2.Expression_patterns{exp}{pat}{1}(3)~=cl)&&(NewNet2.Expression_patterns{exp}{pat}{1}(4)~=cl)&&(NewNet2.Expression_patterns{exp}{pat}{1}(5)==cl) % If the modulatory input is being divided    
                    if ismember(NewNet2.Expression_patterns{exp}{1},[NewNet2.Expression_patterns{exp}{pat}{1}(5) NewNet2.Expression_patterns{exp}{pat}{1}(5)+1]) % If the dividing class is the one that expresses the pattern (and therefore the expressing class is either "original" or "copy")
                        if NewNet2.Expression_patterns{exp}{1}==cl+1 % If we are operating in the "copy" of the divided subclass
                            NewNet2.Expression_patterns{exp}{pat}{1}(5)=cl+1; % Edit the triplet within the link to match the new expressing subclass
                        elseif NewNet2.Expression_patterns{exp}{1}==cl % If we are operating in the "original" of the divided subclass, do nothing
                        else error(['Impossible situation during the expression duplication, the situation is: cl=' num2str(cl) ' , e_cl=' num2str(NewNet2.Expression_patterns{exp}{1}) ', x1=' num2str(NewNet2.Expression_patterns{exp}{pat}{1}(3)) ', y1=' num2str(NewNet2.Expression_patterns{exp}{pat}{1}(4)) ', z1=' num2str(NewNet2.Expression_patterns{exp}{pat}{1}(5))])
                        end
                        
                    elseif (ismember(NewNet2.Expression_patterns{exp}{1},[NewNet2.Expression_patterns{exp}{pat}{1}(3) NewNet2.Expression_patterns{exp}{pat}{1}(3)+1]))||(ismember(NewNet2.Expression_patterns{exp}{1},[NewNet2.Expression_patterns{exp}{pat}{1}(4) NewNet2.Expression_patterns{exp}{pat}{1}(4)+1])) % The expressing subclass is target postsynapse 
                        % The expressing class was not divided, therefore no pattern duplication; but the link duplicates
                        NextIndex=size(NewNet2.Expression_patterns{exp},2)+1; % Get the index of last pattern of the expression to append
                        NewNet2.Expression_patterns{exp}(NextIndex)=NewNet2.Expression_patterns{exp}(pat); % Append the expression with the new pattern (link)
                        NewNet2.Expression_patterns{exp}{NextIndex}{1}(5)=NewNet2.Expression_patterns{exp}{NextIndex}{1}(5)+1; % Modify the copy of the link to fit the copy of the subclass
                        if isstruct(transition_tracker)
                            transition_tracker.links{transition_tracker.expressions(2,:)==exp}(1,end+1)=transition_tracker.links{transition_tracker.expressions(2,:)==exp}(1,transition_tracker.links{transition_tracker.expressions(2,:)==exp}(2,:)==pat); % getting the "original" duplicating link index and setting it as the orginal for the "new" link
                            transition_tracker.links{transition_tracker.expressions(2,:)==exp}(2,end)=NextIndex; % The new link index
                        end
%                         link_duplication_counter=link_duplication_counter+1;
%                         if link_duplication_counter==1
%                             expressions_counter=expressions_counter+1;
%                             report.links{expressions_counter}{1}=exp;
%                             report.links{expressions_counter}{2}=[pat; NextIndex];
%                         else
%                             report.links{expressions_counter}{2}=[report.links{expressions_counter}{2} [pat; NextIndex]];
%                         end
                    else error(['Impossible situation during the expression duplication, the situation is: cl=' num2str(cl) ' , e_cl=' num2str(NewNet2.Expression_patterns{exp}{1}) ', x1=' num2str(NewNet2.Expression_patterns{exp}{pat}{1}(3)) ', y1=' num2str(NewNet2.Expression_patterns{exp}{pat}{1}(4)) ', z1=' num2str(NewNet2.Expression_patterns{exp}{pat}{1}(5))])
                    end
                    
                elseif (NewNet2.Expression_patterns{exp}{pat}{1}(3)==cl)&&(NewNet2.Expression_patterns{exp}{pat}{1}(4)==cl)&&(NewNet2.Expression_patterns{exp}{pat}{1}(5)~=cl) % If the presynaptic and postsynaptic class is divided, recursive link. 
                    if ismember(NewNet2.Expression_patterns{exp}{1},[NewNet2.Expression_patterns{exp}{pat}{1}(3) NewNet2.Expression_patterns{exp}{pat}{1}(3)+1]) % If the dividing class is the one that expresses the pattern
                        % By this point, the expression itself has been duplicated and the the expressing class of the "copy" is set to cl+1
                        if NewNet2.Expression_patterns{exp}{1}==cl+1 % If we are operating in the "copy" of the divided subclass
                            NextIndex=size(NewNet2.Expression_patterns{exp},2)+1; % Get the index of last pattern of the expression to append
                            NewNet2.Expression_patterns{exp}(NextIndex)=NewNet2.Expression_patterns{exp}(pat); % Append the expression with the new pattern (link) 
                            NewNet2.Expression_patterns{exp}{NextIndex}{1}(4)=NewNet2.Expression_patterns{exp}{NextIndex}{1}(4)+1; % Change the adressation of the presynapse of te target
                            NewNet2.Expression_patterns{exp}{pat}{1}(3)=NewNet2.Expression_patterns{exp}{pat}{1}(3)+1; % For the "original", recursive link, just change the indices of both targets
                            NewNet2.Expression_patterns{exp}{pat}{1}(4)=NewNet2.Expression_patterns{exp}{pat}{1}(4)+1;
                            
                            if isstruct(transition_tracker)
                                transition_tracker.links{transition_tracker.expressions(2,:)==exp}(1,end+1)=transition_tracker.links{transition_tracker.expressions(2,:)==exp}(1,transition_tracker.links{transition_tracker.expressions(2,:)==exp}(2,:)==pat); % getting the "original" duplicating link index and setting it as the orginal for the "new" link
                                transition_tracker.links{transition_tracker.expressions(2,:)==exp}(2,end)=NextIndex; % The new link index
                            end
%                             link_duplication_counter=link_duplication_counter+1;
%                             if link_duplication_counter==1
%                                 expressions_counter=expressions_counter+1;
%                                 report.links{expressions_counter}{1}=exp;
%                                 report.links{expressions_counter}{2}=[pat; NextIndex ];
%                             else
%                                 report.links{expressions_counter}{2}=[report.links{expressions_counter}{2} [pat; NextIndex]];
%                             end
                        elseif NewNet2.Expression_patterns{exp}{1}==cl % If we work inside the "original" divided subclass
                            % The reference to a "recursive" connection does not change 
                            NextIndex=size(NewNet2.Expression_patterns{exp},2)+1; % Get the index of last pattern of the expression to append
                            NewNet2.Expression_patterns{exp}(NextIndex)=NewNet2.Expression_patterns{exp}(pat); % Append the expression with the new pattern (link) 
                            NewNet2.Expression_patterns{exp}{NextIndex}{1}(3)=NewNet2.Expression_patterns{exp}{NextIndex}{1}(3)+1; % Change the adressation of the postsynapse of te target
                            if isstruct(transition_tracker)
                                transition_tracker.links{transition_tracker.expressions(2,:)==exp}(1,end+1)=transition_tracker.links{transition_tracker.expressions(2,:)==exp}(1,transition_tracker.links{transition_tracker.expressions(2,:)==exp}(2,:)==pat); % getting the "original" duplicating link index and setting it as the orginal for the "new" link
                                transition_tracker.links{transition_tracker.expressions(2,:)==exp}(2,end)=NextIndex; % The new link index
                            end
%                             link_duplication_counter=link_duplication_counter+1;
%                             if link_duplication_counter==1
%                                 expressions_counter=expressions_counter+1;
%                                 report.links{expressions_counter}{1}=exp;
%                                 report.links{expressions_counter}{2}=[pat; NextIndex ];
%                             else
%                                 report.links{expressions_counter}{2}=[report.links{expressions_counter}{2} [pat; NextIndex]];
%                             end
                        else error(['Impossible situation during the expression duplication, the situation is: cl=' num2str(cl) ' , e_cl=' num2str(NewNet2.Expression_patterns{exp}{1}) ', x1=' num2str(NewNet2.Expression_patterns{exp}{pat}{1}(3)) ', y1=' num2str(NewNet2.Expression_patterns{exp}{pat}{1}(4)) ', z1=' num2str(NewNet2.Expression_patterns{exp}{pat}{1}(5))])
                        end
                    elseif ismember(NewNet2.Expression_patterns{exp}{1},[NewNet2.Expression_patterns{exp}{pat}{1}(5) NewNet2.Expression_patterns{exp}{pat}{1}(5)+1]) % If the modulatory input is the one that expresses the pattern
                        % Here the link quadruplication should happen
                        NextIndex=size(NewNet2.Expression_patterns{exp},2)+1; % Get the index of last pattern of the expression to append
                        NewNet2.Expression_patterns{exp}(NextIndex)=NewNet2.Expression_patterns{exp}(pat); % Append the expression with the new pattern (link) referring to the postsynaptic subclass only
                        NewNet2.Expression_patterns{exp}{NextIndex}{1}(3)=NewNet2.Expression_patterns{exp}{NextIndex}{1}(3)+1;  % Modify the copy of the link to fit the copy of the subclass
                        NewNet2.Expression_patterns{exp}(NextIndex+1)=NewNet2.Expression_patterns{exp}(pat); % Append the expression with the new pattern (link) referring to the presynaptic subclass only
                        NewNet2.Expression_patterns{exp}{NextIndex+1}{1}(4)=NewNet2.Expression_patterns{exp}{NextIndex+1}{1}(4)+1; % Modify the copy of the link to fit the copy of the subclass
                        NewNet2.Expression_patterns{exp}(NextIndex+2)=NewNet2.Expression_patterns{exp}(pat); % Append the expression with the new pattern (link) referring to postsynaptic and presynaptic subclass only
                        NewNet2.Expression_patterns{exp}{NextIndex+2}{1}(3)=NewNet2.Expression_patterns{exp}{NextIndex+2}{1}(3)+1;  % Modify the copy of the link to fit the copy of the subclass (postsynaptic)
                        NewNet2.Expression_patterns{exp}{NextIndex+2}{1}(4)=NewNet2.Expression_patterns{exp}{NextIndex+2}{1}(4)+1;  % Modify the copy of the link to fit the copy of the subclass (presynaptic)
                        if isstruct(transition_tracker)
                            transition_tracker.links{transition_tracker.expressions(2,:)==exp}(1,end+1:end+3)=transition_tracker.links{transition_tracker.expressions(2,:)==exp}(1,transition_tracker.links{transition_tracker.expressions(2,:)==exp}(2,:)==pat); % getting the "original" duplicating link index and setting it as the orginal for the "new" link
                            transition_tracker.links{transition_tracker.expressions(2,:)==exp}(2,end-2:end)=[NextIndex NextIndex+1 NextIndex+2]; % The new link index
                        end
%                        link_duplication_counter=link_duplication_counter+1;
%                         if link_duplication_counter==1
%                             expressions_counter=expressions_counter+1;
%                             report.links{expressions_counter}{1}=exp;
%                             report.links{expressions_counter}{2}=[pat pat pat; NextIndex NextIndex+1 NextIndex+2];
%                         else
%                             report.links{expressions_counter}{2}=[report.links{expressions_counter}{2} [pat pat pat; NextIndex NextIndex+1 NextIndex+2]];
%                         end
                    else error(['Impossible situation during the expression duplication, the situation is: cl=' num2str(cl) ' , e_cl=' num2str(NewNet2.Expression_patterns{exp}{1}) ', x1=' num2str(NewNet2.Expression_patterns{exp}{pat}{1}(3)) ', y1=' num2str(NewNet2.Expression_patterns{exp}{pat}{1}(4)) ', z1=' num2str(NewNet2.Expression_patterns{exp}{pat}{1}(5))])
                    end
                elseif (NewNet2.Expression_patterns{exp}{pat}{1}(3)==cl)&&(NewNet2.Expression_patterns{exp}{pat}{1}(4)~=cl)&&(NewNet2.Expression_patterns{exp}{pat}{1}(5)==cl) % If the postsynaptic target and the modulatory input is the same and it was divided.
                    if ismember(NewNet2.Expression_patterns{exp}{1},[NewNet2.Expression_patterns{exp}{pat}{1}(3) NewNet2.Expression_patterns{exp}{pat}{1}(3)+1]) % If the dividing class is the one that expresses the pattern
                        % Duplication of expression, then the link duplicates. This is the only place where presynaptic/postsynaptic matters directly, because the the dividing subclass can be both and it cannot be determined from the subclass itself. 
                        if NewNet2.Expression_patterns{exp}{1}==cl+1 % If we are operating in the "copy" of the divided subclass
                            NextIndex=size(NewNet2.Expression_patterns{exp},2)+1; % Get the index of last pattern of the expression to append
                            NewNet2.Expression_patterns{exp}(NextIndex)=NewNet2.Expression_patterns{exp}(pat); % Append the expression with the new pattern (link) 
                            if (ismember(NewNet2.Expression_patterns{exp}{pat}{1}(1),[8]) && ismember(NewNet2.Expression_patterns{exp}{pat}{1}(2),[1 3]))||(ismember(NewNet2.Expression_patterns{exp}{pat}{1}(1),[9]) && ismember(NewNet2.Expression_patterns{exp}{pat}{1}(6),[1])) % If postsynaptic
                                NewNet2.Expression_patterns{exp}{NextIndex}{1}(3)=NewNet2.Expression_patterns{exp}{NextIndex}{1}(3)+1; % Change the adressation of the postsynapse of the target
                            else
                                NewNet2.Expression_patterns{exp}{NextIndex}{1}(5)=NewNet2.Expression_patterns{exp}{NextIndex}{1}(5)+1; % Change the adressation of the modulatory input
                            end
                            NewNet2.Expression_patterns{exp}{pat}{1}(3)=NewNet2.Expression_patterns{exp}{pat}{1}(3)+1; % For the "original", recursive link, just change the indices of both targets
                            NewNet2.Expression_patterns{exp}{pat}{1}(5)=NewNet2.Expression_patterns{exp}{pat}{1}(5)+1;
                            if isstruct(transition_tracker)
                                transition_tracker.links{transition_tracker.expressions(2,:)==exp}(1,end+1)=transition_tracker.links{transition_tracker.expressions(2,:)==exp}(1,transition_tracker.links{transition_tracker.expressions(2,:)==exp}(2,:)==pat); % getting the "original" duplicating link index and setting it as the orginal for the "new" link
                                transition_tracker.links{transition_tracker.expressions(2,:)==exp}(2,end)=NextIndex; % The new link index
                            end
%                             link_duplication_counter=link_duplication_counter+1;
%                             if link_duplication_counter==1
%                                 expressions_counter=expressions_counter+1;
%                                 report.links{expressions_counter}{1}=exp;
%                                 report.links{expressions_counter}{2}=[pat; NextIndex ];
%                             else
%                                 report.links{expressions_counter}{2}=[report.links{expressions_counter}{2} [pat; NextIndex]];
%                             end
                        elseif NewNet2.Expression_patterns{exp}{1}==cl % If we work inside the "original" divided subclass
                            % The reference to a "recursive" connection does not change 
                            NextIndex=size(NewNet2.Expression_patterns{exp},2)+1; % Get the index of last pattern of the expression to append
                            NewNet2.Expression_patterns{exp}(NextIndex)=NewNet2.Expression_patterns{exp}(pat); % Append the expression with the new pattern (link) 
                            if (ismember(NewNet2.Expression_patterns{exp}{pat}{1}(1),[8]) && ismember(NewNet2.Expression_patterns{exp}{pat}{1}(2),[1 3]))||(ismember(NewNet2.Expression_patterns{exp}{pat}{1}(1),[9]) && ismember(NewNet2.Expression_patterns{exp}{pat}{1}(6),[1])) % If postsynaptic
                                NewNet2.Expression_patterns{exp}{NextIndex}{1}(5)=NewNet2.Expression_patterns{exp}{NextIndex}{1}(5)+1; % Change the adressation of the modulatory input
                            else
                                NewNet2.Expression_patterns{exp}{NextIndex}{1}(3)=NewNet2.Expression_patterns{exp}{NextIndex}{1}(3)+1; % Change the adressation of the postsynapse of the target
                            end
                            if isstruct(transition_tracker)
                                transition_tracker.links{transition_tracker.expressions(2,:)==exp}(1,end+1)=transition_tracker.links{transition_tracker.expressions(2,:)==exp}(1,transition_tracker.links{transition_tracker.expressions(2,:)==exp}(2,:)==pat); % getting the "original" duplicating link index and setting it as the orginal for the "new" link
                                transition_tracker.links{transition_tracker.expressions(2,:)==exp}(2,end)=NextIndex; % The new link index
                            end
%                             link_duplication_counter=link_duplication_counter+1;
%                             if link_duplication_counter==1
%                                 expressions_counter=expressions_counter+1;
%                                 report.links{expressions_counter}{1}=exp;
%                                 report.links{expressions_counter}{2}=[pat; NextIndex ];
%                             else
%                                 report.links{expressions_counter}{2}=[report.links{expressions_counter}{2} [pat; NextIndex]];
%                             end
                        else error(['Impossible situation during the expression duplication, the situation is: cl=' num2str(cl) ' , e_cl=' num2str(NewNet2.Expression_patterns{exp}{1}) ', x1=' num2str(NewNet2.Expression_patterns{exp}{pat}{1}(3)) ', y1=' num2str(NewNet2.Expression_patterns{exp}{pat}{1}(4)) ', z1=' num2str(NewNet2.Expression_patterns{exp}{pat}{1}(5))])
                        end
                    elseif ismember(NewNet2.Expression_patterns{exp}{1},[NewNet2.Expression_patterns{exp}{pat}{1}(4) NewNet2.Expression_patterns{exp}{pat}{1}(4)+1]) % If the target presynapse expresses 
                        % Here the link quadruplication should happen
                        NextIndex=size(NewNet2.Expression_patterns{exp},2)+1; % Get the index of last pattern of the expression to append
                        NewNet2.Expression_patterns{exp}(NextIndex)=NewNet2.Expression_patterns{exp}(pat); % Append the expression with the new pattern (link) referring to the postsynaptic subclass only
                        NewNet2.Expression_patterns{exp}{NextIndex}{1}(3)=NewNet2.Expression_patterns{exp}{NextIndex}{1}(3)+1;  % Modify the copy of the link to fit the copy of the subclass
                        NewNet2.Expression_patterns{exp}(NextIndex+1)=NewNet2.Expression_patterns{exp}(pat); % Append the expression with the new pattern (link) referring to the presynaptic subclass only
                        NewNet2.Expression_patterns{exp}{NextIndex+1}{1}(5)=NewNet2.Expression_patterns{exp}{NextIndex+1}{1}(5)+1; % Modify the copy of the link to fit the copy of the subclass
                        NewNet2.Expression_patterns{exp}(NextIndex+2)=NewNet2.Expression_patterns{exp}(pat); % Append the expression with the new pattern (link) referring to postsynaptic and presynaptic subclass only
                        NewNet2.Expression_patterns{exp}{NextIndex+2}{1}(3)=NewNet2.Expression_patterns{exp}{NextIndex+2}{1}(3)+1;  % Modify the copy of the link to fit the copy of the subclass (postsynaptic)
                        NewNet2.Expression_patterns{exp}{NextIndex+2}{1}(5)=NewNet2.Expression_patterns{exp}{NextIndex+2}{1}(5)+1;  % Modify the copy of the link to fit the copy of the subclass (modulation)
                        if isstruct(transition_tracker)
                            transition_tracker.links{transition_tracker.expressions(2,:)==exp}(1,end+1:end+3)=transition_tracker.links{transition_tracker.expressions(2,:)==exp}(1,transition_tracker.links{transition_tracker.expressions(2,:)==exp}(2,:)==pat); % getting the "original" duplicating link index and setting it as the orginal for the "new" link
                            transition_tracker.links{transition_tracker.expressions(2,:)==exp}(2,end-2:end)=[NextIndex NextIndex+1 NextIndex+2]; % The new link index
                        end
%                         link_duplication_counter=link_duplication_counter+1;
%                         if link_duplication_counter==1
%                             expressions_counter=expressions_counter+1;
%                             report.links{expressions_counter}{1}=exp;
%                             report.links{expressions_counter}{2}=[pat pat pat; NextIndex NextIndex+1 NextIndex+2];
%                         else
%                             report.links{expressions_counter}{2}=[report.links{expressions_counter}{2} [pat pat pat; NextIndex NextIndex+1 NextIndex+2]];
%                         end
                    else error(['Impossible situation during the expression duplication, the situation is: cl=' num2str(cl) ' , e_cl=' num2str(NewNet2.Expression_patterns{exp}{1}) ', x1=' num2str(NewNet2.Expression_patterns{exp}{pat}{1}(3)) ', y1=' num2str(NewNet2.Expression_patterns{exp}{pat}{1}(4)) ', z1=' num2str(NewNet2.Expression_patterns{exp}{pat}{1}(5))])
                    end
                end
                
%                 % If the index of one (or more) subclasses that form the triplet is less than the index of the divided sublcass, do nothing. These actions are done independently for both indices in the pair. 
%                 if (NewNet2.Expression_patterns{exp}{pat}{1}(3)==cl)&&(NewNet2.Expression_patterns{exp}{pat}{1}(4)~=cl)&&(NewNet2.Expression_patterns{exp}{pat}{1}(5)~=cl) % If the postsynaptic subclass of the pair is being divided
%                     NextIndex=size(NewNet2.Expression_patterns{exp},2)+1; % Get the index of last pattern of the expression to append
%                     NewNet2.Expression_patterns{exp}(NextIndex)=NewNet2.Expression_patterns{exp}(pat); % Append the expression with the new pattern (link)
%                     NewNet2.Expression_patterns{exp}{NextIndex}{1}(3)=NewNet2.Expression_patterns{exp}{NextIndex}{1}(3)+1; % Modify the copy of the link to fit the copy of the subclass
%                 elseif (NewNet2.Expression_patterns{exp}{pat}{1}(3)~=cl)&&(NewNet2.Expression_patterns{exp}{pat}{1}(4)==cl)&&(NewNet2.Expression_patterns{exp}{pat}{1}(5)~=cl) % If the presynaptic subclass of the pair is being divided
%                     NextIndex=size(NewNet2.Expression_patterns{exp},2)+1; % Get the index of last pattern of the expression to append
%                     NewNet2.Expression_patterns{exp}(NextIndex)=NewNet2.Expression_patterns{exp}(pat); % Append the expression with the new pattern (link)
%                     NewNet2.Expression_patterns{exp}{NextIndex}{1}(4)=NewNet2.Expression_patterns{exp}{NextIndex}{1}(4)+1; % Modify the copy of the link to fit the copy of the subclass
%                 elseif (NewNet2.Expression_patterns{exp}{pat}{1}(3)~=cl)&&(NewNet2.Expression_patterns{exp}{pat}{1}(4)~=cl)&&(NewNet2.Expression_patterns{exp}{pat}{1}(5)==cl) % If the modulating subclass of the pair is being divided
%                     % In this case, the expression duplicates, but no duplication of the link is required. However, it is convenient to edit the triplet indexation here
%                     if NewNet2.Expression_patterns{exp}{1}==cl+1 % If we are working in the "copy" of the divided subclass
%                         NewNet2.Expression_patterns{exp}{pat}{1}(5)=cl+1; % The index of the modulatory input of the triplet should match the expressing subclass
%                     end
%                 elseif (NewNet2.Expression_patterns{exp}{pat}{1}(3)==cl)&&(NewNet2.Expression_patterns{exp}{pat}{1}(4)==cl)&&(NewNet2.Expression_patterns{exp}{pat}{1}(5)~=cl)  % If the divided subclass is both post- and pre- synaptic
%                     if (ismember(NewNet2.Expression_patterns{exp}{pat}{1}(1),[8]) && ismember(NewNet2.Expression_patterns{exp}{pat}{1}(2),[1 3]))||(ismember(NewNet2.Expression_patterns{exp}{pat}{1}(1),[9]) && ismember(NewNet2.Expression_patterns{exp}{pat}{1}(6),[1])) % If postsynaptic
%                         if (NewNet2.Expression_patterns{exp}{1}==NewNet2.Expression_patterns{exp}{pat}{1}(3))||(NewNet2.Expression_patterns{exp}{1}==NewNet2.Expression_patterns{exp}{pat}{1}(3)+1) % If the expression is in the recurrent modulated "pair", the "original" or "copy"
%                             if NewNet2.Expression_patterns{exp}{1}==cl % If we work inside the "original" divided subclass
%                                 % The first link stays as it is (the indexation is already eddited)
%                                 NextIndex=size(NewNet2.Expression_patterns{exp},2)+1; % Get the index of last pattern of the expression to append
%                                 NewNet2.Expression_patterns{exp}(NextIndex)=NewNet2.Expression_patterns{exp}(pat); % Append the expression with the new pattern (link) 
%                                 NewNet2.Expression_patterns{exp}{NextIndex}{1}(4)=NewNet2.Expression_patterns{exp}{NextIndex}{1}(4)+1; % Change the adressation of the ionotropic input of te target
%                             elseif NewNet2.Expression_patterns{exp}{1}==cl+1 % If we work in the "copy" of the duvided subclass
%                                 NewNet2.Expression_patterns{exp}{pat}{1}(3)=NewNet2.Expression_patterns{exp}{pat}{1}(3)+1; % For the "original", recursive link, just change the indices of both targets
%                                 NewNet2.Expression_patterns{exp}{pat}{1}(4)=NewNet2.Expression_patterns{exp}{pat}{1}(4)+1;
%                                 NextIndex=size(NewNet2.Expression_patterns{exp},2)+1; % Get the index of last pattern of the expression to append
%                                 NewNet2.Expression_patterns{exp}(NextIndex)=NewNet2.Expression_patterns{exp}(pat); % Append the expression with the new pattern (link) 
%                                 NewNet2.Expression_patterns{exp}{NextIndex}{1}(3)=NewNet2.Expression_patterns{exp}{NextIndex}{1}(3)+1; % Change the adressation of the ionotropic input of te target
%                             else
%                                 error('A problem with modulation triplet adressation during the class division')
%                             end
%                            И СНОВА ЛАЖА В ЛОГИКЕ
%                         else % If it not within the modulated pair but is postsynaptic, then this cannot be
%                             error('A problem with modulation triplet adressation during the class division')
%                         end
%                       else %If presynaptic
%                         if (NewNet2.Expression_patterns{exp}{1}==NewNet2.Expression_patterns{exp}{pat}{1}(5)) % If the expression is in the modulatory input
%                             NextIndex=size(NewNet2.Expression_patterns{exp},2)+1; % Get the index of last pattern of the expression to append
%                             NewNet2.Expression_patterns{exp}(NextIndex)=NewNet2.Expression_patterns{exp}(pat); % Append the expression with the new pattern (link) referring to the postsynaptic subclass only
%                             NewNet2.Expression_patterns{exp}{NextIndex}{1}(3)=NewNet2.Expression_patterns{exp}{NextIndex}{1}(3)+1;  % Modify the copy of the link to fit the copy of the subclass
%                             NewNet2.Expression_patterns{exp}(NextIndex+1)=NewNet2.Expression_patterns{exp}(pat); % Append the expression with the new pattern (link) referring to the presynaptic subclass only
%                             NewNet2.Expression_patterns{exp}{NextIndex+1}{1}(4)=NewNet2.Expression_patterns{exp}{NextIndex+1}{1}(4)+1; % Modify the copy of the link to fit the copy of the subclass
%                             NewNet2.Expression_patterns{exp}(NextIndex+2)=NewNet2.Expression_patterns{exp}(pat); % Append the expression with the new pattern (link) referring to postsynaptic and presynaptic subclass only
%                             NewNet2.Expression_patterns{exp}{NextIndex+2}{1}(3)=NewNet2.Expression_patterns{exp}{NextIndex+2}{1}(3)+1;  % Modify the copy of the link to fit the copy of the subclass (postsynaptic)
%                             NewNet2.Expression_patterns{exp}{NextIndex+2}{1}(4)=NewNet2.Expression_patterns{exp}{NextIndex+2}{1}(4)+1;  % Modify the copy of the link to fit the copy of the subclass (presynaptic)
%                         else error('A problem with modulation triplet adressation during the class division')
%                         end
%                     end
                  
%                 elseif (NewNet2.Expression_patterns{exp}{pat}{1}(3)==cl)&&(NewNet2.Expression_patterns{exp}{pat}{1}(4)~=cl)&&(NewNet2.Expression_patterns{exp}{pat}{1}(5)==cl)  % If the divided subclass is both postsynaptic in the target and the modulatory input
%                     % Such situation means that the divided subclass modulates its own input. The division causes duplication of the expression pattern, but within each patter the link must be duplicated as well. Moreover, the duplication of the pattern changed the
%                     % addressation within the expression (NewNet2.Expression_patterns{exp}{1}), but not within the triplet, so we must fix it here
%                     NextIndex=size(NewNet2.Expression_patterns{exp},2)+1; % Get the index of last pattern of the expression to append
%                     NewNet2.Expression_patterns{exp}(NextIndex)=NewNet2.Expression_patterns{exp}(pat); % Append the expression with the new pattern (link)
%                     if NewNet2.Expression_patterns{exp}{1}==cl % If we are now working within the "original" divided subclass
%                         % The "original" link stays the same
%                         if (ismember(NewNet2.Expression_patterns{exp}{pat}{1}(1),[8]) && ismember(NewNet2.Expression_patterns{exp}{pat}{1}(2),[1 3]))||(ismember(NewNet2.Expression_patterns{exp}{pat}{1}(1),[9]) && ismember(NewNet2.Expression_patterns{exp}{pat}{1}(6),[1])) % If postsynaptic
%                             NewNet2.Expression_patterns{exp}{NextIndex}{1}(5)=cl+1;
%                         else
%                             NewNet2.Expression_patterns{exp}{NextIndex}{1}(3)=cl+1;
%                         end
%                     else % If we are now working within the "copy" of the divided subclass
%                         NewNet2.Expression_patterns{exp}{pat}{1}(3)=cl+1; % Edit the adressation within the triplet to address the "recursive" modulation in the "copy"
%                         NewNet2.Expression_patterns{exp}{pat}{1}(5)=cl+1;
%                         if (ismember(NewNet2.Expression_patterns{exp}{pat}{1}(1),[8]) && ismember(NewNet2.Expression_patterns{exp}{pat}{1}(2),[1 3]))||(ismember(NewNet2.Expression_patterns{exp}{pat}{1}(1),[9]) && ismember(NewNet2.Expression_patterns{exp}{pat}{1}(6),[1])) % If postsynaptic
%                             NewNet2.Expression_patterns{exp}{NextIndex}{1}(3)=cl+1;
%                         else
%                             NewNet2.Expression_patterns{exp}{NextIndex}{1}(5)=cl+1;
%                         end
%                                            
%                     end         
            elseif NewNet2.Expression_patterns{exp}{pat}{1}(1)==6 % A special subroutine for the plasticity references. It is a bit tricky because of a specific data structure for the plasticity
                if ~isempty(Newplast) % If there are any new plasticity patterns
                    Plast_targets=find(Newplast(1,:)==NewNet2.Expression_patterns{exp}{pat}{1}(4)); % Find the ones that were divided or quadruplicated because of the subclass division 
                    
                    if (size(Plast_targets,2)==1) % If the plasticity was divided and if the copy is presented in the currently processing subclass, there are basically four options
                        if (NewNet2.Expression_patterns{exp}{pat}{1}(3)==1) && (NewNet2.Plasticity{Newplast(1,Plast_targets)}{3}(1)==cl) % If the plasticity is postsynaptic and the divided subclass is postsynaptic in the pair
                            if NewNet2.Expression_patterns{exp}{1}~=cl % if we are not in the "original" copy of the divided subclass
                                NewNet2.Expression_patterns{exp}{pat}{1}(4)=Newplast(2,Plast_targets); % Assigning new plasticity to the link, no link division required
                            end
                        elseif (NewNet2.Expression_patterns{exp}{pat}{1}(3)==1) && (NewNet2.Plasticity{Newplast(1,Plast_targets)}{3}(2)==cl) % If the plasticity is postsynaptic but the divided subclass is presynaptic in the pair
                            NextIndex=size(NewNet2.Expression_patterns{exp},2)+1; % Get the index of last pattern of the expression to append
                            NewNet2.Expression_patterns{exp}(NextIndex)=NewNet2.Expression_patterns{exp}(pat); % Append the expression with the new pattern (link)
                            NewNet2.Expression_patterns{exp}{NextIndex}{1}(4)=Newplast(2,Plast_targets); % Modify the copy of the link to refer new link to a new plasticity 
                            if isstruct(transition_tracker)
                                transition_tracker.links{transition_tracker.expressions(2,:)==exp}(1,end+1)=transition_tracker.links{transition_tracker.expressions(2,:)==exp}(1,transition_tracker.links{transition_tracker.expressions(2,:)==exp}(2,:)==pat); % getting the "original" duplicating link index and setting it as the orginal for the "new" link
                                transition_tracker.links{transition_tracker.expressions(2,:)==exp}(2,end)=NextIndex; % The new link index
                            end
%                             link_duplication_counter=link_duplication_counter+1;
%                             if link_duplication_counter==1
%                                 expressions_counter=expressions_counter+1;
%                                 report.links{expressions_counter}{1}=exp;
%                                 report.links{expressions_counter}{2}=[pat; NextIndex ];
%                             else
%                                 report.links{expressions_counter}{2}=[report.links{expressions_counter}{2} [pat; NextIndex]];
%                             end
                        elseif (NewNet2.Expression_patterns{exp}{pat}{1}(3)==2) && (NewNet2.Plasticity{Newplast(1,Plast_targets)}{3}(2)==cl)  % If the plasticity is presynaptic and the divided subclass is presynaptic in the pair
                            if NewNet2.Expression_patterns{exp}{1}~=cl % if we are not in the "original" copy of the divided subclass
                                NewNet2.Expression_patterns{exp}{pat}{1}(4)=Newplast(2,Plast_targets); % Assigning new plasticity to the link, no link division required
                            end
                        elseif (NewNet2.Expression_patterns{exp}{pat}{1}(3)==2) && (NewNet2.Plasticity{Newplast(1,Plast_targets)}{3}(1)==cl)  % If the plasticity is presynaptic but the divided subclass is postesynaptic in the pair
                            NextIndex=size(NewNet2.Expression_patterns{exp},2)+1; % Get the index of last pattern of the expression to append
                            NewNet2.Expression_patterns{exp}(NextIndex)=NewNet2.Expression_patterns{exp}(pat); % Append the expression with the new pattern (link)
                            NewNet2.Expression_patterns{exp}{NextIndex}{1}(4)=Newplast(2,Plast_targets); % Modify the copy of the link to refer new link to a new plasticity
                            if isstruct(transition_tracker)
                                transition_tracker.links{transition_tracker.expressions(2,:)==exp}(1,end+1)=transition_tracker.links{transition_tracker.expressions(2,:)==exp}(1,transition_tracker.links{transition_tracker.expressions(2,:)==exp}(2,:)==pat); % getting the "original" duplicating link index and setting it as the orginal for the "new" link
                                transition_tracker.links{transition_tracker.expressions(2,:)==exp}(2,end)=NextIndex; % The new link index
                            end
%                             link_duplication_counter=link_duplication_counter+1;
%                             if link_duplication_counter==1
%                                 expressions_counter=expressions_counter+1;
%                                 report.links{expressions_counter}{1}=exp;
%                                 report.links{expressions_counter}{2}=[pat; NextIndex ];
%                             else
%                                 report.links{expressions_counter}{2}=[report.links{expressions_counter}{2} [pat; NextIndex]];
%                             end
                        else error('NEVER SUPPOSED TO HAPPEN')
                        end
                    elseif size(Plast_targets,2)==3 % If the plasticity was quadruplicated
                        % The specifics of this situation is: no matter if the plasticity was presynaptic or postsynaptic, if if was in the recurrent connection that connects the neurons of the same subclass, we now have two subclasses, in each of which there is a recurrent connection and a new external connection 
                        % We now have four plasticities instead of one, and two expressions instead of one. This creates a bit of a complex situation.
                        if NewNet2.Expression_patterns{exp}{1}==cl % If the expression we are working with is within the "original" subclass 
                            NextIndex=size(NewNet2.Expression_patterns{exp},2)+1; % Get the index of last pattern of the expression to append
                            NewNet2.Expression_patterns{exp}(NextIndex)=NewNet2.Expression_patterns{exp}(pat); % Append the expression with the new pattern (link)
                            if (NewNet2.Expression_patterns{exp}{pat}{1}(3)==1) % If the plasticity is postsynaptic
                                NewNet2.Expression_patterns{exp}{NextIndex}{1}(4)=Newplast(2,Plast_targets(2)); % The new link is appended by the "second" copy of the quadruplicated plasticity. This is an input from the "copy" subclass. 
                            else % If the plasticity is presynaptic
                                NewNet2.Expression_patterns{exp}{NextIndex}{1}(4)=Newplast(2,Plast_targets(1)); % The new link is appended by the "first" copy of the quadruplicated plasticity. This is an output to the "copy" sublass
                            end
                            if isstruct(transition_tracker)
                                transition_tracker.links{transition_tracker.expressions(2,:)==exp}(1,end+1)=transition_tracker.links{transition_tracker.expressions(2,:)==exp}(1,transition_tracker.links{transition_tracker.expressions(2,:)==exp}(2,:)==pat); % getting the "original" duplicating link index and setting it as the orginal for the "new" link
                                transition_tracker.links{transition_tracker.expressions(2,:)==exp}(2,end)=NextIndex; % The new link index
                            end
%                            link_duplication_counter=link_duplication_counter+1;
%                             if link_duplication_counter==1
%                                 expressions_counter=expressions_counter+1;
%                                 report.links{expressions_counter}{1}=exp;
%                                 report.links{expressions_counter}{2}=[pat; NextIndex ];
%                             else
%                                 report.links{expressions_counter}{2}=[report.links{expressions_counter}{2} [pat; NextIndex]];
%                             end
                        elseif NewNet2.Expression_patterns{exp}{1}==cl+1 % If the expression we are working with is within the "copy" of the subclass
                            NextIndex=size(NewNet2.Expression_patterns{exp},2)+1; % Get the index of last pattern of the expression to append
                            NewNet2.Expression_patterns{exp}(NextIndex)=NewNet2.Expression_patterns{exp}(pat); % Append the expression with the new pattern (link)
                            if (NewNet2.Expression_patterns{exp}{pat}{1}(3)==1) % If the plasticity is postsynaptic
                                NewNet2.Expression_patterns{exp}{NextIndex}{1}(4)=Newplast(2,Plast_targets(1)); % The new link is appended by the "first" copy of the quadruplicated plasticity. This is an input from the "original" subclass
                            else % If the plasticity is presynaptic
                                NewNet2.Expression_patterns{exp}{NextIndex}{1}(4)=Newplast(2,Plast_targets(2)); % The new link is appended by the "second" copy of the quadruplicated plasticity. This is an output to the "original" subclass
                            end
                            NewNet2.Expression_patterns{exp}{pat}{1}(4)=Newplast(2,Plast_targets(3)); % Alter the plasticity index in the original "link" with the "third" copy of the plasticity (in other words, edit the recurrent connection)
                            if isstruct(transition_tracker)
                                transition_tracker.links{transition_tracker.expressions(2,:)==exp}(1,end+1)=transition_tracker.links{transition_tracker.expressions(2,:)==exp}(1,transition_tracker.links{transition_tracker.expressions(2,:)==exp}(2,:)==pat); % getting the "original" duplicating link index and setting it as the orginal for the "new" link
                                transition_tracker.links{transition_tracker.expressions(2,:)==exp}(2,end)=NextIndex; % The new link index
                            end
%                             link_duplication_counter=link_duplication_counter+1;
%                             if link_duplication_counter==1
%                                 expressions_counter=expressions_counter+1;
%                                 report.links{expressions_counter}{1}=exp;
%                                 report.links{expressions_counter}{2}=[pat; NextIndex ];
%                             else
%                                 report.links{expressions_counter}{2}=[report.links{expressions_counter}{2} [pat; NextIndex]];
%                             end
                        else % If the expression we are working with is within some other subclass, basically we must do the following: copy every link that refers to the connection with divided subclass.
                            % NOTE: within the plasticities, we already have a new indexation for the subclasses. So we don't need to change the references to the plasticities here, because even if the indexation of the subclasses changed, the "old" plasticities refer to the correct subclasses.  
                            warning('WHY AM I HERE? CHECK LEGALITY OF THE SOURCE GENOTYPE')
                        end
                         
                        %plast_exp_summary(NewNet2)
                    elseif size(Plast_targets,2)==0 % The plasticity pattern targeted by the current link was not affected by the class division, do nothing
                        
                        
                    else % This should never happen
                        exp
                        pat
                        Newplast
                        Plast_targets
                        error('Impossible situation with the plasticity and expression patterns during the class division ')
                    end
                end
            end    
                
            end
        
            
        end
    end
end




%legality(NewNet2)

