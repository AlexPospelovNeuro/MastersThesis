
% The function that introduces the point mutations into the Net2 genotype

% INPUTS:
% - A properly structured Net2 genotype
% - A Point_mutation structure that contains probabilities and magnitudes for the mutations of the Net2 parameters
% OUTPUTS:
% - A new Ne2 genotype, structured similarly to the input (same set of the interneuron classes, plasticity patterns and expression patterns)


function NewNet2=point_mutation(Net2,Point_Mut)

%legality(Net2) % Initial legality test

NewNet2=Net2;
    % Ionotropic interneurons number
    for cl=1:size(NewNet2.Cells.Ion,2)
        if rand(1)<Point_Mut.Cells.Ion(1,1)
            NewNet2.Cells.Ion(1,cl)=max(1,NewNet2.Cells.Ion(1,cl)+randn*Point_Mut.Cells.Ion(1,2)); % At this version of the function, the mean number of neurons per subclass cannot become less than 1. 
        end
        if rand(1)<Point_Mut.Cells.Ion(2,1)
            NewNet2.Cells.Ion(2,cl)=max(0,NewNet2.Cells.Ion(2,cl)+randn*Point_Mut.Cells.Ion(2,2));
        end
    end
    % Modulatory interneurons number
    for cl=1:size(NewNet2.Cells.Mod,2)
        if rand(1)<Point_Mut.Cells.Mod(1,1)
            NewNet2.Cells.Mod(1,cl)=max(1,NewNet2.Cells.Mod(1,cl)+randn*Point_Mut.Cells.Mod(1,2));
        end
        if rand(1)<Point_Mut.Cells.Mod(2,1)
            NewNet2.Cells.Mod(2,cl)=max(0,NewNet2.Cells.Mod(2,cl)+randn*Point_Mut.Cells.Mod(2,2));
        end
    end
    % Thresholds, threshold noise, recurrent connections, absolute refracterity
    for cl=1:size(NewNet2.Thresholds,2) 
        if rand(1)<Point_Mut.Thresholds(1,1)
            NewNet2.Thresholds(1,cl)=NewNet2.Thresholds(1,cl)+randn*Point_Mut.Thresholds(1,2);
        end
        if rand(1)<Point_Mut.Thresholds(2,1)
            NewNet2.Thresholds(2,cl)=max(0,NewNet2.Thresholds(2,cl)+randn*Point_Mut.Thresholds(2,2));
        end
        if rand(1)<Point_Mut.ThreshNoise(1,1)
            NewNet2.ThreshNoise(1,cl)=NewNet2.ThreshNoise(1,cl)+randn*Point_Mut.ThreshNoise(1,2);
        end
        if rand(1)<Point_Mut.ThreshNoise(2,1)
            NewNet2.ThreshNoise(2,cl)=max(0,NewNet2.ThreshNoise(2,cl)+randn*Point_Mut.ThreshNoise(2,2));
        end
        if rand(1)<Point_Mut.AbsRefract(1,1)
            NewNet2.AbsRefract(1,cl)=NewNet2.AbsRefract(1,cl)+randn*Point_Mut.AbsRefract(1,2);
        end
        if rand(1)<Point_Mut.AbsRefract(2,1)
            NewNet2.AbsRefract(2,cl)=max(0,NewNet2.AbsRefract(2,cl)+randn*Point_Mut.AbsRefract(2,2));
        end
        if rand(1)<Point_Mut.RecurCon(1,1)
            NewNet2.RecurCon(1,cl)=NewNet2.RecurCon(1,cl)+randn*Point_Mut.RecurCon(1,2);
        end
        if rand(1)<Point_Mut.AbsRefract(2,1)
            NewNet2.RecurCon(2,cl)=max(0,NewNet2.RecurCon(2,cl)+randn*Point_Mut.RecurCon(2,2));
        end
    end
    % Ionotropic connections
    for post_cl=size(NewNet2.Cells.Input,2)+1:size(NewNet2.Connections,1)
        for pre_cl=1:size(NewNet2.Connections,2)
            if rand(1)<Point_Mut.Connections(1,1)
                NewNet2.Connections(post_cl,pre_cl,1)=NewNet2.Connections(post_cl,pre_cl,1)+randn*Point_Mut.Connections(1,2); % mean
            end
            if rand(1)<Point_Mut.Connections(2,1)
                NewNet2.Connections(post_cl,pre_cl,2)=max(0,NewNet2.Connections(post_cl,pre_cl,2)+randn*Point_Mut.Connections(2,2)); % SD in
            end
            if rand(1)<Point_Mut.Connections(3,1)
                NewNet2.Connections(post_cl,pre_cl,3)=max(0,NewNet2.Connections(post_cl,pre_cl,3)+randn*Point_Mut.Connections(3,2)); % SD out
            end
        end
    end
    % Basic powers
    for post_cl=1:size(NewNet2.BasicPowers,1)
        for pre_cl=1:size(NewNet2.BasicPowers,2)
            if rand(1)<Point_Mut.BasicPowers(1,1)
                NewNet2.BasicPowers(post_cl,pre_cl,1)=NewNet2.BasicPowers(post_cl,pre_cl,1)+randn*Point_Mut.BasicPowers(1,2); % mean
            end
            if rand(1)<Point_Mut.BasicPowers(2,1)
                NewNet2.BasicPowers(post_cl,pre_cl,2)=max(0,NewNet2.BasicPowers(post_cl,pre_cl,2)+randn*Point_Mut.BasicPowers(2,2)); % SD 
            end
            
        end
    end
    % Delays
    for post_cl=1:size(NewNet2.Delays,1)
        for pre_cl=1:size(NewNet2.Delays,2)
            if rand(1)<Point_Mut.Delays(1,1)
                NewNet2.Delays(post_cl,pre_cl,1)=NewNet2.Delays(post_cl,pre_cl,1)+randn*Point_Mut.Delays(1,2); % mean
            end
            if rand(1)<Point_Mut.Delays(2,1)
                NewNet2.Delays(post_cl,pre_cl,2)=max(0,NewNet2.Delays(post_cl,pre_cl,2)+randn*Point_Mut.Delays(2,2)); % SD 
            end
            
        end
    end
    % PSP shapes
    for post_cl=1:size(NewNet2.PSPshape,1)
        for pre_cl=1:size(NewNet2.PSPshape,2)
            if rand(1)<Point_Mut.PSPshape(1,1)
                NewNet2.PSPshape{post_cl,pre_cl}(1,1)=NewNet2.PSPshape{post_cl,pre_cl}(1,1)+randn*Point_Mut.PSPshape(1,2); % mean
            end
            if rand(1)<Point_Mut.PSPshape(2,1)
                NewNet2.PSPshape{post_cl,pre_cl}(2,1)=max(0,NewNet2.PSPshape{post_cl,pre_cl}(2,1)+randn*Point_Mut.PSPshape(2,2)); % SD 
            end
            if rand(1)<Point_Mut.PSPshape(1,1)
                NewNet2.PSPshape{post_cl,pre_cl}(1,2)=NewNet2.PSPshape{post_cl,pre_cl}(1,2)+randn*Point_Mut.PSPshape(1,2); % mean
            end
            if rand(1)<Point_Mut.PSPshape(2,1)
                NewNet2.PSPshape{post_cl,pre_cl}(2,2)=max(0,NewNet2.PSPshape{post_cl,pre_cl}(2,2)+randn*Point_Mut.PSPshape(2,2)); % SD 
            end
            if rand(1)<Point_Mut.PSPshape(1,1)
                NewNet2.PSPshape{post_cl,pre_cl}(1,3)=NewNet2.PSPshape{post_cl,pre_cl}(1,3)+randn*Point_Mut.PSPshape(1,2); % mean
            end
            if rand(1)<Point_Mut.PSPshape(2,1)
                NewNet2.PSPshape{post_cl,pre_cl}(2,3)=max(0,NewNet2.PSPshape{post_cl,pre_cl}(2,3)+randn*Point_Mut.PSPshape(2,2)); % SD 
            end
            
        end
    end
    
    % Modulation connections
    for post_cl=1:size(NewNet2.Connections_Mod,1)
        for pre_cl=1:size(NewNet2.Connections_Mod,2)
            for mod_cl=1:size(NewNet2.Connections_Mod,3)
                if rand(1)<Point_Mut.Connections_Mod(1,1)
                    NewNet2.Connections_Mod(post_cl,pre_cl,mod_cl,1)=NewNet2.Connections_Mod(post_cl,pre_cl,mod_cl,1)+randn*Point_Mut.Connections_Mod(1,2); % mean
                end
                if rand(1)<Point_Mut.Connections_Mod(2,1)
                    NewNet2.Connections_Mod(post_cl,pre_cl,mod_cl,2)=max(0,NewNet2.Connections_Mod(post_cl,pre_cl,mod_cl,2)+randn*Point_Mut.Connections_Mod(2,2)); % SD in
                end
                if rand(1)<Point_Mut.Connections(3,1)
                    NewNet2.Connections_Mod(post_cl,pre_cl,mod_cl,3)=max(0,NewNet2.Connections_Mod(post_cl,pre_cl,mod_cl,3)+randn*Point_Mut.Connections_Mod(3,2)); % SD out
                end
            end
        end
    end
    
    % Modulation parameters
    for post_cl=1:size(NewNet2.Mod,1)
        for pre_cl=1:size(NewNet2.Mod,2)
            for mod_cl=1:size(NewNet2.Mod,3)
                % Single trigger effect
                if rand(1)<Point_Mut.Mod{1}(1,1)
                    NewNet2.Mod{post_cl,pre_cl,mod_cl}{1}(1,1)=NewNet2.Mod{post_cl,pre_cl,mod_cl}{1}(1,1)+randn*Point_Mut.Mod{1}(1,2); % mean
                end
                if rand(1)<Point_Mut.Mod{1}(2,1)
                    NewNet2.Mod{post_cl,pre_cl,mod_cl}{1}(2,1)=max(0,NewNet2.Mod{post_cl,pre_cl,mod_cl}{1}(2,1)+randn*Point_Mut.Mod{1}(2,2)); % SD 
                end
                % Maximal number of triggers
                if rand(1)<Point_Mut.Mod{2}(1,1)
                    NewNet2.Mod{post_cl,pre_cl,mod_cl}{2}(1,1)=NewNet2.Mod{post_cl,pre_cl,mod_cl}{2}(1,1)+randn*Point_Mut.Mod{2}(1,2); % mean
                end
                if rand(1)<Point_Mut.Mod{2}(2,1)
                    NewNet2.Mod{post_cl,pre_cl,mod_cl}{2}(2,1)=max(0,NewNet2.Mod{post_cl,pre_cl,mod_cl}{2}(2,1)+randn*Point_Mut.Mod{2}(2,2)); % SD 
                end
                % temporal characteristics
                if rand(1)<Point_Mut.Mod{3}(1,1)
                    NewNet2.Mod{post_cl,pre_cl,mod_cl}{3}(1,1)=NewNet2.Mod{post_cl,pre_cl,mod_cl}{3}(1,1)+randn*Point_Mut.Mod{3}(1,2); % mean
                end
                if rand(1)<Point_Mut.Mod{3}(2,1)
                    NewNet2.Mod{post_cl,pre_cl,mod_cl}{3}(2,1)=max(0,NewNet2.Mod{post_cl,pre_cl,mod_cl}{3}(2,1)+randn*Point_Mut.Mod{3}(2,2)); % SD 
                end
                
                if rand(1)<Point_Mut.Mod{3}(3,1)
                    NewNet2.Mod{post_cl,pre_cl,mod_cl}{3}(1,2)=NewNet2.Mod{post_cl,pre_cl,mod_cl}{3}(1,2)+randn*Point_Mut.Mod{3}(3,2); % mean
                end
                if rand(1)<Point_Mut.Mod{3}(4,1)
                    NewNet2.Mod{post_cl,pre_cl,mod_cl}{3}(2,2)=max(0,NewNet2.Mod{post_cl,pre_cl,mod_cl}{3}(2,2)+randn*Point_Mut.Mod{3}(4,2)); % SD 
                end
                
                if rand(1)<Point_Mut.Mod{3}(5,1)
                    NewNet2.Mod{post_cl,pre_cl,mod_cl}{3}(1,3)=NewNet2.Mod{post_cl,pre_cl,mod_cl}{3}(1,3)+randn*Point_Mut.Mod{3}(5,2); % mean
                end
                if rand(1)<Point_Mut.Mod{3}(6,1)
                    NewNet2.Mod{post_cl,pre_cl,mod_cl}{3}(2,3)=max(0,NewNet2.Mod{post_cl,pre_cl,mod_cl}{3}(2,3)+randn*Point_Mut.Mod{3}(6,2)); % SD 
                end
                
                if rand(1)<Point_Mut.Mod{3}(7,1)
                    NewNet2.Mod{post_cl,pre_cl,mod_cl}{3}(1,4)=NewNet2.Mod{post_cl,pre_cl,mod_cl}{3}(1,4)+randn*Point_Mut.Mod{3}(7,2); % mean
                end
                if rand(1)<Point_Mut.Mod{3}(8,1)
                    NewNet2.Mod{post_cl,pre_cl,mod_cl}{3}(2,4)=max(0,NewNet2.Mod{post_cl,pre_cl,mod_cl}{3}(2,4)+randn*Point_Mut.Mod{3}(8,2)); % SD 
                end
            end
        end
    end 
    
    % Plasticity patterns
    if isfield(NewNet2,'Plasticity')
        for pl=1:size(NewNet2.Plasticity,2)
            if ~isnan(NewNet2.Plasticity{pl}{2}) %If there is the observation frame specified
                if rand(1)<Point_Mut.Plasticity{2}(1,1)
                    NewNet2.Plasticity{pl}{2}(1,1)=NewNet2.Plasticity{pl}{2}(1,1)+randn*Point_Mut.Plasticity{2}(1,2); % mean, Beginning of the observation frame
                end
                if rand(1)<Point_Mut.Plasticity{2}(2,1)
                    NewNet2.Plasticity{pl}{2}(2,1)=max(0,NewNet2.Plasticity{pl}{2}(2,1)+randn*Point_Mut.Plasticity{2}(2,2)); % SD 
                end
                
                if rand(1)<Point_Mut.Plasticity{2}(1,1)
                    NewNet2.Plasticity{pl}{2}(1,2)=NewNet2.Plasticity{pl}{2}(1,2)+randn*Point_Mut.Plasticity{2}(1,2); % mean End of the observation frame
                end
                if rand(1)<Point_Mut.Plasticity{2}(2,1)
                    NewNet2.Plasticity{pl}{2}(2,2)=max(0,NewNet2.Plasticity{pl}{2}(2,2)+randn*Point_Mut.Plasticity{2}(2,2)); % SD 
                end
            end
            
            if rand(1)<Point_Mut.Plasticity{4}{1}(1,1)
                NewNet2.Plasticity{pl}{4}{1}(1,1)=NewNet2.Plasticity{pl}{4}{1}(1,1)+randn*Point_Mut.Plasticity{4}{1}(1,2); % mean, single trigger effect
            end
            if rand(1)<Point_Mut.Plasticity{4}{1}(2,1)
                NewNet2.Plasticity{pl}{4}{1}(2,1)=max(0,NewNet2.Plasticity{pl}{4}{1}(2,1)+randn*Point_Mut.Plasticity{4}{1}(2,2)); % SD 
            end
            if rand(1)<Point_Mut.Plasticity{4}{2}(1,1)
                NewNet2.Plasticity{pl}{4}{2}(1,1)=NewNet2.Plasticity{pl}{4}{2}(1,1)+randn*Point_Mut.Plasticity{4}{2}(1,2); % mean, maximal number of triggers
            end
            if rand(1)<Point_Mut.Plasticity{4}{2}(2,1)
                NewNet2.Plasticity{pl}{4}{2}(2,1)=max(0,NewNet2.Plasticity{pl}{4}{2}(2,1)+randn*Point_Mut.Plasticity{4}{2}(2,2)); % SD 
            end
            
            % temporal characteristics
            if rand(1)<Point_Mut.Plasticity{4}{3}(1,1)
                NewNet2.Plasticity{pl}{4}{3}(1,1)=NewNet2.Plasticity{pl}{4}{3}(1,1)+randn*Point_Mut.Plasticity{4}{3}(1,2); % mean
            end
            if rand(1)<Point_Mut.Plasticity{4}{3}(2,1)
                NewNet2.Plasticity{pl}{4}{3}(2,1)=max(0,NewNet2.Plasticity{pl}{4}{3}(2,1)+randn*Point_Mut.Plasticity{4}{3}(2,2)); % SD 
            end
                
            if rand(1)<Point_Mut.Plasticity{4}{3}(3,1)
                NewNet2.Plasticity{pl}{4}{3}(1,2)=NewNet2.Plasticity{pl}{4}{3}(1,2)+randn*Point_Mut.Plasticity{4}{3}(3,2); % mean
            end
            if rand(1)<Point_Mut.Plasticity{4}{3}(4,1)
                NewNet2.Plasticity{pl}{4}{3}(2,2)=max(0,NewNet2.Plasticity{pl}{4}{3}(2,2)+randn*Point_Mut.Plasticity{4}{3}(4,2)); % SD 
            end
                
            if rand(1)<Point_Mut.Plasticity{4}{3}(5,1)
                NewNet2.Plasticity{pl}{4}{3}(1,3)=NewNet2.Plasticity{pl}{4}{3}(1,3)+randn*Point_Mut.Plasticity{4}{3}(5,2); % mean
            end
            if rand(1)<Point_Mut.Plasticity{4}{3}(6,1)
                NewNet2.Plasticity{pl}{4}{3}(2,3)=max(0,NewNet2.Plasticity{pl}{4}{3}(2,3)+randn*Point_Mut.Plasticity{4}{3}(6,2)); % SD 
            end
                
            if rand(1)<Point_Mut.Plasticity{4}{3}(7,1)
                NewNet2.Plasticity{pl}{4}{3}(1,4)=NewNet2.Plasticity{pl}{4}{3}(1,4)+randn*Point_Mut.Plasticity{4}{3}(7,2); % mean
            end
            if rand(1)<Point_Mut.Plasticity{4}{3}(8,1)
                NewNet2.Plasticity{pl}{4}{3}(2,4)=max(0,NewNet2.Plasticity{pl}{4}{3}(2,4)+randn*Point_Mut.Plasticity{4}{3}(8,2)); % SD 
            end
            
            
        end
    
    end
    
    
    
    % Expression patterns
    if isfield(NewNet2,'Expression_patterns')
        for exp=1:size(NewNet2.Expression_patterns,2) % for each expression
            for pat=2:size(NewNet2.Expression_patterns{exp},2) % and each pattern inside it
                if rand(1)<Point_Mut.Expression_patterns(1,1)
                    NewNet2.Expression_patterns{exp}{pat}{2}(1)=NewNet2.Expression_patterns{exp}{pat}{2}(1)+Point_Mut.Expression_patterns(1,2);
                end
                if rand(1)<Point_Mut.Expression_patterns(2,1)
                    NewNet2.Expression_patterns{exp}{pat}{2}(2)=NewNet2.Expression_patterns{exp}{pat}{2}(2)+Point_Mut.Expression_patterns(2,2);
                end
                if rand(1)<Point_Mut.Expression_patterns(3,1)
                    NewNet2.Expression_patterns{exp}{pat}{2}(3)=NewNet2.Expression_patterns{exp}{pat}{2}(3)+Point_Mut.Expression_patterns(3,2);
                end
            end
        end
    end




%legality(NewNet2) % Final legality test
end