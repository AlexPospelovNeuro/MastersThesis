% This function tests the legality genotype of the Net2 type of the circuit

% INPUT - Net2 structure.
% Output - None, warnings in the command window.

function legality(Net2)

Genotype_legal=1;

%%%%% SUBCLASSES STRUCTURE %%%%%
if ~isfield(Net2,"Cells")
    warning('Genotype is illegal: no Cells field, information about subclasses is missing')
    Genotype_legal=0;
end
if ~isfield(Net2.Cells,"Input")
    warning('Genotype is illegal: no Cells.Input field, information about Input subclass is missing')
    Genotype_legal=0;
end
if ~isfield(Net2.Cells,"Output")
    warning('Genotype is illegal: no Cells.Output field, information about Output subclass is missing')
    Genotype_legal=0;
end
if ~isfield(Net2.Cells,"Ion")
    warning('Genotype is illegal: no Cells.Ion field, information about Ion subclass is missing')
    Genotype_legal=0;
end
if ~isfield(Net2.Cells,"Mod")
    warning('Genotype is illegal: no Cells.Mod field, information about Mod subclass is missing')
    Genotype_legal=0;
end

if size(Net2.Cells.Input,1)~=2
    warning('Genotype is illegal: field Net2.Cells.Input has incorrect format; correct format: [[means];[sigmas]]')
    Genotype_legal=0;
end
if size(Net2.Cells.Output,1)~=2
    warning('Genotype is illegal: field Net2.Cells.Output has incorrect format; correct format: [[means];[sigmas]]')
    Genotype_legal=0;
end
if size(Net2.Cells.Ion,1)~=2
    warning('Genotype is illegal: field Net2.Cells.Ion has incorrect format; correct format: [[means];[sigmas]]')
    Genotype_legal=0;
end
if size(Net2.Cells.Mod,1)~=2
    warning('Genotype is illegal: field Net2.Cells.Mod has incorrect format; correct format: [[means];[sigmas]]')
    Genotype_legal=0;
end
n_subclasses=size(Net2.Cells.Input,2)+size(Net2.Cells.Output,2)+size(Net2.Cells.Ion,2)+size(Net2.Cells.Mod,2);
n_nonmod=size(Net2.Cells.Input,2)+size(Net2.Cells.Output,2)+size(Net2.Cells.Ion,2);
%%%%% CONNECTIONS %%%%%
if ~isfield(Net2,"Connections")
    warning('Genotype is illegal: no Connections field, information about connections is missing')
    Genotype_legal=0;
end
if size(Net2.Connections,1)~=n_subclasses
    warning('Genotype is illegal: number of rows in the Connection matrix does not match the number of subclasses')
    Genotype_legal=0;
end
if size(Net2.Connections,2)~=n_nonmod
    warning('Genotype is illegal: number of columns in the Connection matrix does not match the number of nonmodulatory subclasses')
    Genotype_legal=0;
end
if size(Net2.Connections,3)~=3
    warning('Genotype is illegal: incorrect format for the connection parameters; for each pair of neurons, there should be [mean postsynaptic_sigma presynaptic_sigma]')
    Genotype_legal=0;
end
for input_class=1:size(Net2.Cells.Input,2)
    if (sum(Net2.Connections(input_class,:,1))>0) || (sum(Net2.Connections(input_class,:,2))>0) || sum(Net2.Connections(input_class,:,3))>0
        warning(['Genotype is illegal: The input class ' num2str(input_class) ' has nonzero mean or sd of the input connections'])
        Genotype_legal=0;
    end
end
if ~isfield(Net2,"Connections_Mod")
    warning('Genotype is illegal: no Connections_Mod field, information about modulatory connections is missing')
    Genotype_legal=0;
end
if size(Net2.Connections_Mod,1)~=n_subclasses
    warning('Genotype is illegal: number of rows in the Connection_Mod matrix does not match the number of subclasses')
    Genotype_legal=0;
end
if size(Net2.Connections_Mod,2)~=n_nonmod
    warning('Genotype is illegal: number of columns in the Connection_Mod matrix does not match the number of nonmodulatory subclasses')
    Genotype_legal=0;
end
if size(Net2.Connections_Mod,3)~=n_subclasses-n_nonmod
    warning('Genotype is illegal: third dimention of the Connection_Mod matrix does not match the number of modulatory subclasses')
    Genotype_legal=0;
end
if size(Net2.Connections_Mod,4)~=3
    warning('Genotype is illegal: incorrect format for the modulatory connection parameters; for each triplet of neurons, there should be [mean postsynaptic_sigma presynaptic_sigma]')
    Genotype_legal=0;
end


%%%%% POWERS AND DELAYS %%%%%
if ~isfield(Net2,"BasicPowers")
    warning('Genotype is illegal: no BasicPowers field, information about basic powers is missing')
    Genotype_legal=0;
end
if size(Net2.BasicPowers,1)~=n_subclasses
    warning('Genotype is illegal: number of rows in the BasicPowers matrix does not match the number of subclasses')
    Genotype_legal=0;
end
if size(Net2.BasicPowers,2)~=n_nonmod
    warning('Genotype is illegal: number of columns in the BasicPowers matrix does not match the number of nonmodulatory subclasses')
    Genotype_legal=0;
end
if ~isfield(Net2,"Delays")
    warning('Genotype is illegal: no Delays field, information about delays is missing')
    Genotype_legal=0;
end
if size(Net2.Delays,1)~=n_subclasses
    warning('Genotype is illegal: number of rows in the Delays matrix does not match the number of subclasses')
    Genotype_legal=0;
end
if size(Net2.Delays,2)~=n_nonmod
    warning('Genotype is illegal: number of columns in the Delays matrix does not match the number of nonmodulatory subclasses')
    Genotype_legal=0;
end

%%%%% PSP SHAPES %%%%%
if ~isfield(Net2,"PSPshape")
    warning('Genotype is illegal: no PSPshape field, information about PSP shapes is missing')
    Genotype_legal=0;
end
if size(Net2.PSPshape,1)~=n_subclasses
    warning('Genotype is illegal: number of rows in the PSPshape cell array does not match the number of subclasses')
    Genotype_legal=0;
end
if size(Net2.PSPshape,2)~=n_nonmod
    warning('Genotype is illegal: number of columns in the PSPshape cell array does not match the number of nonmodulatory subclasses')
    Genotype_legal=0;
end
for post=1:n_subclasses
    for pre=1:n_nonmod % The PSP is nonapplicable for the modulatore "presynapse"
        if size(Net2.PSPshape{post,pre},1)~=2
            warning(['Genotype is illegal: Net2.PSPshape{' num2str(post) ',' num2str(pre) '} has incorrect format; correct format: [[means];[sigmas]]'])
            Genotype_legal=0;
        end
        if size(Net2.PSPshape{post,pre},2)~=3
            warning(['Genotype is illegal: Net2.PSPshape{' num2str(post) ',' num2str(pre) '} has incorrect format; correct format: [in plateau out]'])
            Genotype_legal=0;
        end
    end
end

%%%%% THRESHOLDS, NOISE, RECURRENT CONNECTIONS AND ABSOLUTE REFRACTERITY %%%%%
if ~isfield(Net2,"Thresholds")
    warning('Genotype is illegal: no Thresholds field, information about thresholds is missing')
    Genotype_legal=0;
end
if size(Net2.Thresholds,1)~=2
    warning('Genotype is illegal: field Net2.Thresholds has incorrect format; correct format: [[means];[sigmas]]')
    Genotype_legal=0;
end
if size(Net2.Thresholds,2)~=n_subclasses
    warning('Genotype is illegal: field Net2.Thresholds does not match the number of subclasses')
    Genotype_legal=0;
end
if ~isfield(Net2,"AbsRefract")
    warning('Genotype is illegal: no AbsRefract field, information about absolute refracterity is missing')
    Genotype_legal=0;
end
if size(Net2.AbsRefract,1)~=2
    warning('Genotype is illegal: field Net2.AbsRefract has incorrect format; correct format: [[means];[sigmas]]')
    Genotype_legal=0;
end
if size(Net2.AbsRefract,2)~=n_subclasses
    warning('Genotype is illegal: field Net2.AbsRefract does not match the number of subclasses')
    Genotype_legal=0;
end
if ~isfield(Net2,"ThreshNoise")
    warning('Genotype is illegal: no ThreshNoise field, information about threshold noise is missing')
    Genotype_legal=0;
end
if size(Net2.ThreshNoise,1)~=2
    warning('Genotype is illegal: field Net2.ThreshNoise has incorrect format; correct format: [[means];[sigmas]]')
    Genotype_legal=0;
end
if size(Net2.ThreshNoise,2)~=n_subclasses
    warning('Genotype is illegal: field Net2.ThreshNoise does not match the number of subclasses')
    Genotype_legal=0;
end
if ~isfield(Net2,"RecurCon")
    warning('Genotype is illegal: no RecurCon field, information about recurrent connections is missing')
    Genotype_legal=0;
end
if size(Net2.RecurCon,1)~=2
    warning('Genotype is illegal: field Net2.RecurCon has incorrect format; correct format: [[means];[sigmas]]')
    Genotype_legal=0;
end
if size(Net2.RecurCon,2)~=n_subclasses
    warning('Genotype is illegal: field Net2.RecurCon does not match the number of subclasses')
    Genotype_legal=0;
end


%%%%% MODULATION PATTERNS %%%%%

if ~isfield(Net2,"Mod")
    warning('Genotype is illegal: no MOD field, information about modulation patterns is missing')
    Genotype_legal=0;
end
if size(Net2.Mod,1)~=n_subclasses
    warning('Genotype is illegal: number of rows in the Mod cell array does not match the number of subclasses')
    Genotype_legal=0;
end
if size(Net2.Mod,2)~=n_nonmod
    warning('Genotype is illegal: number of columns in the Mod cell array does not match the number of nonmodulatory subclasses')
    Genotype_legal=0;
end
if size(Net2.Mod,3)~=n_subclasses-n_nonmod
    warning('Genotype is illegal: third dimention of the Mod cell array does not match the number of modulatory subclasses')
    Genotype_legal=0;
end
for post=1:n_subclasses
    for pre=1:n_nonmod
        for mod_c=1:n_subclasses-n_nonmod
            if size(Net2.Mod{post,pre,mod_c},1)~=1
                warning(['Genotype is illegal: Net2.Mod{' num2str(post) ',' num2str(pre) ',' num2str(mod_c) '} has incorrect format; should be 1x3 cell array'])
                Genotype_legal=0;
            end
            if size(Net2.Mod{post,pre,mod_c},2)~=3
                warning(['Genotype is illegal: Net2.Mod{' num2str(post) ',' num2str(pre) ',' num2str(mod_c) '} has incorrect format; should be 1x3 cell array'])
                Genotype_legal=0;
            end
            if (size(Net2.Mod{post,pre,mod_c}{1},1)~=2)&&(size(Net2.Mod{post,pre,mod_c}{1},2)~=1)
                warning(['Genotype is illegal: Net2.Mod{' num2str(post) ',' num2str(pre) ',' num2str(mod_c) '} has incorrect format; the first field should contain [[mean];[sigma]] vector for single effects'])
                Genotype_legal=0;
            end
            if (size(Net2.Mod{post,pre,mod_c}{2},1)~=2)&&(size(Net2.Mod{post,pre,mod_c}{2},2)~=1)
                warning(['Genotype is illegal: Net2.Mod{' num2str(post) ',' num2str(pre) ',' num2str(mod_c) '} has incorrect format; the second field should contain [[mean];[sigma]] vector for maximal number of effects'])
                Genotype_legal=0;
            end
            if (size(Net2.Mod{post,pre,mod_c}{3},1)~=2)&&(size(Net2.Mod{post,pre,mod_c}{3},2)~=4)
                warning(['Genotype is illegal: Net2.Mod{' num2str(post) ',' num2str(pre) ',' num2str(mod_c) '} has incorrect format; the second field should contain [[means];[sigmas]] vector for [delay in plateau out]'])
                Genotype_legal=0;
            end
        end
    end
end

%%%%% PLASTICITY PATTERNS

if ~isfield(Net2,"Plasticity")
    warning('No plasticity patterns detected, check if they should be in the genotype')
else
    if size(Net2.Plasticity,1)~=1
        warning('Genotype is illegal: the plasticity field should contain 1xPl cell array, but the number of rows is more than 1 now')
        Genotype_legal=0;
    end

    n_plast=size(Net2.Plasticity,2); % number of plasticity patterns
    legal_plasticity_codes=[1 2 3 1000]; % legal types of plasticity in the current Net2 genotype version
    for pl=1:n_plast
        if (size(Net2.Plasticity{pl}{1},1)~=1) && (size(Net2.Plasticity{pl}{1},2)~=1)
            warning(['Genotype is illegal: the plasticity pattern ' num2str(pl) ' has incorrect format for the type: should be a single integer '])
            Genotype_legal=0;
        end
        if ~ismember(Net2.Plasticity{pl}{1},legal_plasticity_codes) 
            warning(['Genotype is illegal: the plasticity pattern ' num2str(pl) ' belongs to an illegal type ' num2str(Net2.Plasticity{pl}{1})])
            Genotype_legal=0;
        end
        if ismember(Net2.Plasticity{pl}{1},[1 1000])
            if (size(Net2.Plasticity{pl}{2},1)~=1) && (size(Net2.Plasticity{pl}{2},2)~=1) 
                if (~isnan(Net2.Plasticity{pl}{2}))
                    warning(['Genotype is illegal: the plasticity pattern ' num2str(pl) ' belongs to a type with undefined observation window: the second value of the pattern should be [NaN] '])
                    Genotype_legal=0;
                end
            end
        end
        if ismember(Net2.Plasticity{pl}{1},[2 3])
            if (size(Net2.Plasticity{pl}{2},1)~=2) && (size(Net2.Plasticity{pl}{2},2)~=2)
                warning(['Genotype is illegal: the plasticity pattern ' num2str(pl) ' belongs to a type with defined observation window: the second value of the pattern should be 2x2 numeric matrix with parameters of observation window'])
                Genotype_legal=0;
            end
        end
        if ismember(Net2.Plasticity{pl}{1},[1000])
            if (size(Net2.Plasticity{pl}{3},1)~=1)
                warning(['Genotype is illegal: the plasticity pattern ' num2str(pl) ' belongs to a type assigned to a subclass: the third element should contain a single number referring to a cell class'])
                Genotype_legal=0;
            end
            if ~ismember(Net2.Plasticity{pl}{3}(1),1:n_subclasses) 
                warning(['Genotype is illegal: the plasticity pattern ' num2str(pl) ' refers to an undefined cell subclass'])
                Genotype_legal=0;
            end
        end
        if ismember(Net2.Plasticity{pl}{1},[1 2 3])
            if (size(Net2.Plasticity{pl}{3},1)~=1) && (size(Net2.Plasticity{pl}{3},2)~=2)
                warning(['Genotype is illegal: the plasticity pattern ' num2str(pl) ' belongs to a type assigned to a pair of subclasses: the third element should contain a 1x2 vector referring to [post pre] subclasses'])
                Genotype_legal=0;
            end
            if ~ismember(Net2.Plasticity{pl}{3}(1),1:n_subclasses) 
                warning(['Genotype is illegal: the plasticity pattern ' num2str(pl) ' refers to an undefined cell subclass of postsynaptic neurons'])
                Genotype_legal=0;
            end
            if ~ismember(Net2.Plasticity{pl}{3}(2),1:n_subclasses) 
                warning(['Genotype is illegal: the plasticity pattern ' num2str(pl) ' refers to an undefined cell subclass of presynaptic neurons'])
                Genotype_legal=0;
            end
            if ismember(Net2.Plasticity{pl}{3}(2),n_nonmod+1:n_subclasses) 
                warning(['Genotype is illegal: the plasticity pattern ' num2str(pl) ' refers to a modulatory subclass, which is illegal'])
                Genotype_legal=0;
            end
        end
        if (size(Net2.Plasticity{pl}{4}{1},1)~=2) && (size(Net2.Plasticity{pl}{4}{1},2)~=1)
            warning(['Genotype is illegal: the plasticity pattern ' num2str(pl) ' has an illegal format for the single trigger power: should be [mean; sigma]'])
            Genotype_legal=0;
        end
        if (size(Net2.Plasticity{pl}{4}{2},1)~=2) && (size(Net2.Plasticity{pl}{4}{2},2)~=1)
            warning(['Genotype is illegal: the plasticity pattern ' num2str(pl) ' has an illegal format for the maximal number of triggers: should be [mean; sigma]'])
            Genotype_legal=0;
        end
        if (size(Net2.Plasticity{pl}{4}{3},1)~=2) && (size(Net2.Plasticity{pl}{4}{3},2)~=4)
            warning(['Genotype is illegal: the plasticity pattern ' num2str(pl) ' has an illegal format for the timings: should be [means; sigmas] for [delay in plateau out]'])
            Genotype_legal=0;
        end
    end
end
% % 1000 - extrasynaptic input-independent activity-dependent (relative refracterity). Does not depend on the input activity, only on the spike of the neuron. Synapse-non specific, specific for a neuron class, not an ordered pair. Independent from an absolute threshold. Does not use the observation frame (check and remove)
% % 1 - synaptic input-dependent activity-independent (pair-pulse facilitation/inhibition). Does not depend on the postsynaptic neuron activity (but may depend on postsynaptic neuron parameters). No parameters of induction, triggered by the event of an input signal
% % 2 - synaptic input-independent activity-dependent (spike-timing plasticity). Triggered by the second spike of the neuron in a short time. Frame of induction - [beginning end], timing for the observation window for the first spike in a pair
% % 3 - synaptic input-dependent activity-dependent (hebbian plasticity). Similar to pair-pulse plasticity, but is triggered by a temporal coincidence of the synaptic event and the postsynaptic spike within a timeframe. Parameters of induction - [beginning end], timing for the observation window for the first event in the pair
% 
% 
% 
% 
% % The plasticity pattern for the synapse is legal only if the presynaptic neuron is non-modulatory
% Net2.Plasticity{1}={1,[NaN],[3 4],{[-0.3;0.1],[10;1],[[0 1 10 5];[0 0 2 1]]}}; % The coordinate in the outer cell array is the pattern's ID (immutable, reassignment is possible if the pattern in eliminated by mutation)
% Net2.Plasticity{2}={1000,[1 6;0 1],[4],{[0.4;0.1],[10;1],[[2 3 20 10];[0 0 2 1]]}};
% Net2.Plasticity{3}={2,[3 6;1 1],[4 2],{[0.4;0.1],[10;1],[[2 2 2 15];[1 0 2 1]]}};
% Net2.Plasticity{4}={3,[3 6;1 1],[4 4],{[0.4;0.1],[10;1],[[2 2 2 15];[1 0 2 1]]}};
% Net2.Plasticity{5}={3,[3 6;1 1],[6 3],{[0.4;0.1],[10;1],[[2 2 2 15];[1 0 2 1]]}};
% Net2.Plasticity{6}={2,[3 6;1 1],[4 4],{[0.4;0.1],[10;1],[[2 2 2 15];[1 0 2 1]]}};
% Net2.Plasticity{7}={1000,[1 6;0 1],[5],{[0.4;0.1],[10;1],[[2 3 20 10];[0 0 2 1]]}};


%%%%% EXPRESSION PATTERNS %%%%%   CONTINUE FROM HERE


if ~isfield(Net2,"Expression_patterns")
    warning('No expression patterns detected, check if they should be in the genotype')
else
    if size(Net2.Expression_patterns,1)~=1
        warning('Genotype is illegal: the Expression_patterns field should contain 1xPl cell array, but the number of rows is more than 1 now')
        Genotype_legal=0;
    end

    n_exp=size(Net2.Expression_patterns,2); % number of expression patterns
    for exp=1:n_exp
        if ~ismember(Net2.Expression_patterns{exp}{1}(1),1:n_subclasses) 
            warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' referes to an undefined neuron subclass'])
            Genotype_legal=0;
        end
        for pat=2:size(Net2.Expression_patterns{exp},2)
            %%%%% THE INTRINSIC LEGALITY OF THE PATTERNS %%%%%
            if ~ismember(Net2.Expression_patterns{exp}{pat}{1}(1),1:11) % OVERALL LEGALITY OF REFERENCE
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' has an unknown feature to influence'])
                Genotype_legal=0;
            end
        
            if (Net2.Expression_patterns{exp}{pat}{1}(1)==1) && (size(Net2.Expression_patterns{exp}{pat}{1},2)~=3) % CONNECTIONS AND AFFINITIES, IONOTROPIC
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' has wrong set of connections reference parameters: must be 3'])
                Genotype_legal=0;
            else
                if (Net2.Expression_patterns{exp}{pat}{1}(1)==1) && (~ismember(Net2.Expression_patterns{exp}{pat}{1}(2),1:4)) 
                    warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' refers to an unknown connection feature'])
                    Genotype_legal=0;
                end
            end

            if (Net2.Expression_patterns{exp}{pat}{1}(1)==4) && (size(Net2.Expression_patterns{exp}{pat}{1},2)~=1) % THRESHOLD 
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' has wrong set of threshold reference parameters: must be 1'])
                Genotype_legal=0;
            end
            if (Net2.Expression_patterns{exp}{pat}{1}(1)==5) && (size(Net2.Expression_patterns{exp}{pat}{1},2)~=1) % ABS REFRACT
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' has wrong set of absolute refracterity reference parameters: must be 1'])
                Genotype_legal=0;
            end
            if (Net2.Expression_patterns{exp}{pat}{1}(1)==10) && (size(Net2.Expression_patterns{exp}{pat}{1},2)~=1) % RECURRENT CONNECTIONS
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' has wrong set of recurrent connections reference parameters: must be 1'])
                Genotype_legal=0;
            end
            if (Net2.Expression_patterns{exp}{pat}{1}(1)==11) && (size(Net2.Expression_patterns{exp}{pat}{1},2)~=1) % THRESHOLD NOISE
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' has wrong set of threshold noise reference parameters: must be 1'])
                Genotype_legal=0;
            end

            if (Net2.Expression_patterns{exp}{pat}{1}(1)==2) && (size(Net2.Expression_patterns{exp}{pat}{1},2)~=3) % BASIC POWERS
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' has wrong set of basic powers reference parameters: must be 3'])
                Genotype_legal=0;
                if (Net2.Expression_patterns{exp}{pat}{1}(1)==2) && (~ismember(Net2.Expression_patterns{exp}{pat}{1}(2),1:2)) 
                    warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' refers to an unknown basic power feature'])
                    Genotype_legal=0;
                end
            end
        
            if (Net2.Expression_patterns{exp}{pat}{1}(1)==3) && (size(Net2.Expression_patterns{exp}{pat}{1},2)~=3) % DELAYS
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' has wrong set of delays reference parameters: must be 3'])
                Genotype_legal=0;
            else
                if (Net2.Expression_patterns{exp}{pat}{1}(1)==3) && (~ismember(Net2.Expression_patterns{exp}{pat}{1}(2),1:2)) 
                    warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' refers to an unknown delays feature'])
                    Genotype_legal=0;
                end
            end
            if isfield(Net2,"Plasticity")
                if (Net2.Expression_patterns{exp}{pat}{1}(1)==6) && (size(Net2.Expression_patterns{exp}{pat}{1},2)~=4) % PLASTICITY PATTERNS
                    warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' has wrong set of plasticity pattern reference parameters: must be 4'])
                    Genotype_legal=0;
                else
                    if (Net2.Expression_patterns{exp}{pat}{1}(1)==6) && (~ismember(Net2.Expression_patterns{exp}{pat}{1}(2),1:8)) 
                        warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' refers to an unknown plasticity pattern feature'])
                        Genotype_legal=0;
                    end
                    if (Net2.Expression_patterns{exp}{pat}{1}(1)==6) && (~ismember(Net2.Expression_patterns{exp}{pat}{1}(3),1:2)) 
                        warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' refers to an unknown plasticity pattern synaptic addressing'])
                    	Genotype_legal=0;
                    end
                end
            end

            if (Net2.Expression_patterns{exp}{pat}{1}(1)==7) && (size(Net2.Expression_patterns{exp}{pat}{1},2)~=3) % PSP FEATURES
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' has wrong set of PSP reference parameters: must be 3'])
                Genotype_legal=0;
            else
                if (Net2.Expression_patterns{exp}{pat}{1}(1)==7) && (~ismember(Net2.Expression_patterns{exp}{pat}{1}(2),1:3)) 
                    warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' refers to an unknown parameter feature'])
                    Genotype_legal=0;
                end
            end
            
            if (Net2.Expression_patterns{exp}{pat}{1}(1)==8) && (size(Net2.Expression_patterns{exp}{pat}{1},2)~=5) % MODULATORY CONNECTIONS
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' has wrong set of modulatory connections reference parameters: must be 4'])
                Genotype_legal=0;
            else
                if (Net2.Expression_patterns{exp}{pat}{1}(1)==8) && (~ismember(Net2.Expression_patterns{exp}{pat}{1}(2),1:4)) 
                    warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' refers to an unknown parameter feature'])
                    Genotype_legal=0;
                end
            end

            if (Net2.Expression_patterns{exp}{pat}{1}(1)==9) && (size(Net2.Expression_patterns{exp}{pat}{1},2)~=6) % MODULATORY effects
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' has wrong set of modulatory effects reference parameters: must be 5'])
                Genotype_legal=0;
            else
                if (Net2.Expression_patterns{exp}{pat}{1}(1)==9) && (~ismember(Net2.Expression_patterns{exp}{pat}{1}(2),1:6)) 
                    warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' refers to an unknown parameter feature'])
                    Genotype_legal=0;
                end
                if (Net2.Expression_patterns{exp}{pat}{1}(1)==9) && (~ismember(Net2.Expression_patterns{exp}{pat}{1}(6),1:2)) 
                    warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' refers to an unknown parameter feature'])
                    Genotype_legal=0;
                end
            end

        

        %%%%% THE LEGALITY WITH RESPECT TO THE OTHER PARTS OF GENOTYPE

            if (Net2.Expression_patterns{exp}{pat}{1}(1)==1) && (~ismember(Net2.Expression_patterns{exp}{pat}{1}(3),1:n_subclasses)) % CONNECTIONS AND AFFINITIES, IONOTROPIC
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' refers to an undefined subclass'])
                Genotype_legal=0;
            end
            if (Net2.Expression_patterns{exp}{pat}{1}(1)==1) && (ismember(Net2.Expression_patterns{exp}{pat}{1}(2),[1 3])) && (ismember(Net2.Expression_patterns{exp}{pat}{1}(3),n_nonmod+1:n_subclasses))
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' affects postsynaptic connections/affinity pattern for connection received from modulatory subclass'])
                Genotype_legal=0;
            end
            if (Net2.Expression_patterns{exp}{pat}{1}(1)==1) && (ismember(Net2.Expression_patterns{exp}{pat}{1}(2),[2 4])) && (ismember(Net2.Expression_patterns{exp}{1}(1),n_nonmod+1:n_subclasses))
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' affects presynaptic connections/affinity pattern for connections sent by modulatory subclass'])
                Genotype_legal=0;
            end

            if (Net2.Expression_patterns{exp}{pat}{1}(1)==2) && (~ismember(Net2.Expression_patterns{exp}{pat}{1}(3),1:n_subclasses)) % BASIC POWERS
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' refers to an undefined subclass'])
                Genotype_legal=0;
            end
            if (Net2.Expression_patterns{exp}{pat}{1}(1)==2) && (ismember(Net2.Expression_patterns{exp}{pat}{1}(2),[1])) && (ismember(Net2.Expression_patterns{exp}{pat}{1}(3),n_nonmod+1:n_subclasses))
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' affects postsynaptic basic powers for connection received from modulatory subclass'])
                Genotype_legal=0;
            end
            if (Net2.Expression_patterns{exp}{pat}{1}(1)==2) && (ismember(Net2.Expression_patterns{exp}{pat}{1}(2),[2])) && (ismember(Net2.Expression_patterns{exp}{1}(1),n_nonmod+1:n_subclasses))
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' affects presynaptic basic powers for connections sent by modulatory subclass'])
                Genotype_legal=0;
            end

            if (Net2.Expression_patterns{exp}{pat}{1}(1)==3) && (~ismember(Net2.Expression_patterns{exp}{pat}{1}(3),1:n_subclasses)) % DELAYS
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' refers to an undefined subclass'])
                Genotype_legal=0;
            end
            if (Net2.Expression_patterns{exp}{pat}{1}(1)==3) && (ismember(Net2.Expression_patterns{exp}{pat}{1}(2),[1])) && (ismember(Net2.Expression_patterns{exp}{pat}{1}(3),n_nonmod+1:n_subclasses))
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' affects postsynaptic delays for connection received from modulatory subclass'])
                Genotype_legal=0;
            end
            if (Net2.Expression_patterns{exp}{pat}{1}(1)==3) && (ismember(Net2.Expression_patterns{exp}{pat}{1}(2),[2])) && (ismember(Net2.Expression_patterns{exp}{1}(1),n_nonmod+1:n_subclasses))
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' affects presynaptic delays for connections sent by modulatory subclass'])
                Genotype_legal=0;
            end
            
            if isfield(Net2,"Plasticity")
%             [exp pat]
%             Net2.Expression_patterns{exp}{pat}{1}
%             n_plast
            if (Net2.Expression_patterns{exp}{pat}{1}(1)==6) && (~ismember(Net2.Expression_patterns{exp}{pat}{1}(4),1:n_plast)) % PLASTICITY PATTERNS
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' refers to an unknown plasticity pattern'])
                Genotype_legal=0;
            else
                if (Net2.Expression_patterns{exp}{pat}{1}(1)==6) && (ismember(Net2.Expression_patterns{exp}{pat}{1}(2),[7 8])) && ismember(Net2.Plasticity{Net2.Expression_patterns{exp}{pat}{1}(4)}{1},[1 1000])
                    warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' refers to observation frame for plasticity pattern that does not have this feature'])
                    Genotype_legal=0;
                end
                if (Net2.Expression_patterns{exp}{pat}{1}(1)==6) && (Net2.Expression_patterns{exp}{pat}{1}(3)==2) && (ismember(Net2.Expression_patterns{exp}{1}(1),n_nonmod+1:n_subclasses))
                    warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' refers presynaptic feature of a plasticity pattern, but the neuron class is modulatory'])
                    Genotype_legal=0;

                end
                if (Net2.Expression_patterns{exp}{pat}{1}(1)==6) && (Net2.Expression_patterns{exp}{pat}{1}(3)==2) && ismember(Net2.Plasticity{Net2.Expression_patterns{exp}{pat}{1}(4)}{1},[1 2 3]) && (Net2.Plasticity{Net2.Expression_patterns{exp}{pat}{1}(4)}{3}(2) ~= Net2.Expression_patterns{exp}{1})
                    warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' refers to presynaptic feature of a plasticity pattern, but is not expressed by a presynaptic neuron of the plasticity displaying pair'])
                    Genotype_legal=0;
                end
                if (Net2.Expression_patterns{exp}{pat}{1}(1)==6) && (Net2.Expression_patterns{exp}{pat}{1}(3)==1) && ismember(Net2.Plasticity{Net2.Expression_patterns{exp}{pat}{1}(4)}{1},[1 2 3]) && (Net2.Plasticity{Net2.Expression_patterns{exp}{pat}{1}(4)}{3}(1) ~= Net2.Expression_patterns{exp}{1})
                    warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' refers to postsynaptic feature of a plasticity pattern, but is not expressed by a postsynaptic neuron of the plasticity displaying pair'])
                    Genotype_legal=0;
                end
            end
            end
            
            if (Net2.Expression_patterns{exp}{pat}{1}(1)==7) && (~ismember(Net2.Expression_patterns{exp}{pat}{1}(3),1:n_subclasses)) % PSP FEATURES
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' refers to an undefined subclass'])
                Genotype_legal=0;
            end
            if (Net2.Expression_patterns{exp}{pat}{1}(1)==7) && (ismember(Net2.Expression_patterns{exp}{pat}{1}(3),n_nonmod+1:n_subclasses)) 
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' affects PSP features, but source subclass is modulatory'])
                Genotype_legal=0;
            end
        
            
            
            
% [8 1 x y z] - Post_connection_pattern for modulation of the synapse x<-y and the modulatory input z 
% [8 2 x y z] - Pre_connection_pattern for modulation of the synapse x<-y and the modulatory input z 
% [8 3 x y z] - Post_affinity_pattern for modulation of the synapse x<-y and the modulatory input z 
% [8 4 x y z] - Pre_affinity_pattern for modulation of the synapse x<-y and the modulatory input z 
% [9 1 x y z 1/2] - A singular modulaton effect. The fifth place: 1 - postsynapse (any subclass, either x or y, modulatory y is forbidden), 2 presynapse (mod). All for synapse x<-y
% [9 2 x y z 1/2] - Max number of triggers. The fifth place: 1 - postsynapse (any subclass, either x or y, modulatory y is forbidden), 2 presynapse (mod)
% [9 3 x y z 1/2] - Delay. The fifth place: 1 - postsynapse (any subclass, either x or y, modulatory y is forbidden), 2 presynapse (mod)
% [9 4 x y z 1/2] - Onset duration. The fifth place: 1 - postsynapse (any subclass, either x or y, modulatory y is forbidden), 2 presynapse (mod)
% [9 5 x y z 1/2] - Plateau duration. The fifth place: 1 - postsynapse (any subclass, either x or y, modulatory y is forbidden), 2 presynapse (mod)
% [9 6 x y z 1/2] - Offset duration. The fifth place: 1 - postsynapse (any subclass, either x or y, modulatory y is forbidden), 2 presynapse (mod)         
            
            
% CHECK IF IT IS ALL CORRECT            
            
            
            
            if (Net2.Expression_patterns{exp}{pat}{1}(1)==8) && (~ismember(Net2.Expression_patterns{exp}{pat}{1}(3),1:n_subclasses)) % MODULATORY CONNECTIONS
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' refers to an undefined subclass (postsynapse of the modulation target)'])
                Genotype_legal=0;
            end
            if (Net2.Expression_patterns{exp}{pat}{1}(1)==8) && (~ismember(Net2.Expression_patterns{exp}{pat}{1}(4),1:n_subclasses)) 
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' refers to an undefined subclass (presynapse of the modulation target)'])
                Genotype_legal=0;
            end
            if (Net2.Expression_patterns{exp}{pat}{1}(1)==8) && (ismember(Net2.Expression_patterns{exp}{pat}{1}(4),n_nonmod+1:n_subclasses)) 
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' refers to a synapse where presynapse is modulatory itself'])
                Genotype_legal=0;
            end
            if (Net2.Expression_patterns{exp}{pat}{1}(1)==8) && (~ismember(Net2.Expression_patterns{exp}{pat}{1}(5),n_nonmod+1:n_subclasses)) 
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' refers to a nonmodulatory subclass while it has to refer to a modulatory subclass'])
                Genotype_legal=0;
            end
            if (Net2.Expression_patterns{exp}{pat}{1}(1)==8) && (ismember(Net2.Expression_patterns{exp}{pat}{1}(2),[2 4])) && (~ismember(Net2.Expression_patterns{exp}{1}(1),n_nonmod+1:n_subclasses))
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' refers to presynaptic connection feature of the modulatory connection, but the expressing subclass is not modulatory'])
                Genotype_legal=0;
            end
            if (Net2.Expression_patterns{exp}{pat}{1}(1)==8) && (ismember(Net2.Expression_patterns{exp}{pat}{1}(2),[1 3])) && ((Net2.Expression_patterns{exp}{1}(1)~=Net2.Expression_patterns{exp}{pat}{1}(3)) && (Net2.Expression_patterns{exp}{1}(1)~=Net2.Expression_patterns{exp}{pat}{1}(4)))
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' refers to postsynaptic connection feature of the modulatory connection, but the expressing subclass does not participate in the connection'])
                Genotype_legal=0;
            end
            if (Net2.Expression_patterns{exp}{pat}{1}(1)==8) &&  (ismember(Net2.Expression_patterns{exp}{pat}{1}(2),[1 3])) && ((Net2.Expression_patterns{exp}{1}~=Net2.Expression_patterns{exp}{pat}{1}(3)) && (Net2.Expression_patterns{exp}{1}~=Net2.Expression_patterns{exp}{pat}{1}(4)))
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' is expressed in the subclass ' num2str(Net2.Expression_patterns{exp}{1}) ' and is postsynaptic, but the expressing subclass is not a member of the modulation target pair' ])
                Genotype_legal=0;
            end
            if (Net2.Expression_patterns{exp}{pat}{1}(1)==8) &&  (ismember(Net2.Expression_patterns{exp}{pat}{1}(2),[2 4])) && (Net2.Expression_patterns{exp}{1}~=Net2.Expression_patterns{exp}{pat}{1}(5))
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' is expressed in the subclass ' num2str(Net2.Expression_patterns{exp}{1}) ' and is presynaptic, but the expressing subclass is not a source of modulation ' ])
                Genotype_legal=0;
            end    
        
            if (Net2.Expression_patterns{exp}{pat}{1}(1)==9) && (~ismember(Net2.Expression_patterns{exp}{pat}{1}(3),1:n_subclasses)) % MODULATORY EFFECTS
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' refers to an undefined subclass (postsynapse of the modulation target)'])
                Genotype_legal=0;
            end
            if (Net2.Expression_patterns{exp}{pat}{1}(1)==9) && (~ismember(Net2.Expression_patterns{exp}{pat}{1}(4),1:n_subclasses)) 
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' refers to an undefined subclass (presynapse of the modulation target)'])
                Genotype_legal=0;
            end
            if (Net2.Expression_patterns{exp}{pat}{1}(1)==9) && (ismember(Net2.Expression_patterns{exp}{pat}{1}(4),n_nonmod+1:n_subclasses)) 
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' refers to a synapse where presynapse is modulatory itself'])
                Genotype_legal=0;
            end
            if (Net2.Expression_patterns{exp}{pat}{1}(1)==9) && (~ismember(Net2.Expression_patterns{exp}{pat}{1}(5),n_nonmod+1:n_subclasses)) 
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' refers to a nonmodulatory subclass while it has to refer to a modulatory subclass'])
                Genotype_legal=0;
            end
            if (Net2.Expression_patterns{exp}{pat}{1}(1)==9) && (ismember(Net2.Expression_patterns{exp}{pat}{1}(6),[1])) && ((Net2.Expression_patterns{exp}{1}(1)~=Net2.Expression_patterns{exp}{pat}{1}(3)) && (Net2.Expression_patterns{exp}{1}(1)~=Net2.Expression_patterns{exp}{pat}{1}(4)))
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' refers to a postsynaptic effect, but the expressing neuron subclass is not postsynaptic'])
                Genotype_legal=0;
            end
            if (Net2.Expression_patterns{exp}{pat}{1}(1)==9) && (ismember(Net2.Expression_patterns{exp}{pat}{1}(6),[2])) && (~ismember(Net2.Expression_patterns{exp}{1}(1),n_nonmod+1:n_subclasses))
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' refers to a presynaptic effect, but the expresing neuron subclass is not modulatory'])
                Genotype_legal=0;
            end
            if (Net2.Expression_patterns{exp}{pat}{1}(1)==9) &&  (ismember(Net2.Expression_patterns{exp}{pat}{1}(6),[1])) && ((Net2.Expression_patterns{exp}{1}~=Net2.Expression_patterns{exp}{pat}{1}(3)) && (Net2.Expression_patterns{exp}{1}~=Net2.Expression_patterns{exp}{pat}{1}(4)))
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' is expressed in the subclass ' num2str(Net2.Expression_patterns{exp}{1}) ' and is postsynaptic, but the expressing subclass is not a member of the modulation target pair' ])
                Genotype_legal=0;
            end
            if (Net2.Expression_patterns{exp}{pat}{1}(1)==9) &&  (ismember(Net2.Expression_patterns{exp}{pat}{1}(6),[2])) && (Net2.Expression_patterns{exp}{1}~=Net2.Expression_patterns{exp}{pat}{1}(5))
                warning(['Genotype is illegal: the expression pattern ' num2str(exp) ' the exact pattern ' num2str(pat-1) ' is expressed in the subclass ' num2str(Net2.Expression_patterns{exp}{1}) ' and is presynaptic, but the expressing subclass is not a source of modulation ' ])
                Genotype_legal=0;
            end    
            
        end
    end
end
% Net2.Expression_patterns{1}={[3],{[1 1 3],[0 1 3]},{[1 2 4],[1 1 0]},{[2 2 4],[1 1 1],{[2 1 4],[1 2 3]}}};
% Net2.Expression_patterns{2}={[3],{[1 1 3],[0 1 2]},{[3 1 3],[3 -1 -2]},{[3 2 3],[3 -1 -2]},{[2 2 3],[3 -1 -2]}};
% Net2.Expression_patterns{3}={[4],{[2 2 4],[0 1 -3]}};
% Net2.Expression_patterns{4}={[4],{[3 1 3],[0 0 -10]},{[2 1 3],[0 5 0]}};
% Net2.Expression_patterns{5}={[4],{[6 1 1 4],[4 20 30]}}; % The pattern to test the plasticity routine (postsynaptic)
% Net2.Expression_patterns{6}={[4],{[6 1 2 6],[1 10 0]},{[2 2 4],[0 0 -10]}}; % The pattern to test the plasticity routine (presynaptic)
% Net2.Expression_patterns{7}={[5],{[9 1 2 6],[1 10 0]},{[8 2 3 4],[0 0 -10]}}; % To test the modulation

% Neuron-specific targets
% [1 1 y] - Post_connection_pattern from the class y
% [1 2 y] - Pre_connection_pattern from the class y
% [1 3 y] - Post_affinity_pattern  from the class y
% [1 4 y] - Pre_affinity_pattern from the class y
% [4] - Thresholds pattern 
% [5] - AbsRefract pattern 


% Synapse-specific targets
% [2 1 y] - Basic_Powers pattern from the class y
% [2 2 y] - Basic_Powers pattern to the class y
% [3 1 y] - Delay pattern from the class y
% [3 2 y] - Delay pattern to the class y 
% [6 1 1/2 t] - For the plasticity pattern with ID t, a single act maximal effect. Third place: 1 for postsynaptic, 2 for presynaptic
% [6 2 1/2 t] - For the plasticity pattern with ID t, Maximal simultaneous actions. Third place: 1 for postsynaptic, 2 for presynaptic
% [6 3 1/2 t] - For the plasticity pattern with ID t, a delay. Third place: 1 for postsynaptic, 2 for presynaptic
% [6 4 1/2 t] - For the plasticity pattern with ID t, an onset duration. Third place: 1 for postsynaptic, 2 for presynaptic
% [6 5 1/2 t] - For the plasticity pattern with ID t, a plateau duration. Third place: 1 for postsynaptic, 2 for presynaptic
% [6 6 1/2 t] - For the plasticity pattern with ID t, an offset duration. Third place: 1 for postsynaptic, 2 for presynaptic
% [6 7 1/2 t] - For the plasticity pattern with ID t, a triggering condition frame beginning. Third place: 1 for postsynaptic, 2 for presynaptic
% [6 8 1/2 t] - For the plasticity pattern with ID t, a triggering condition frame end. Third place: 1 for postsynaptic, 2 for presynaptic
% [7 1 y] - For the PSP from the class y, an onset duration. NOT DONE YET!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% [7 2 y] - For the PSP from the class y, a plateau duration. 
% [7 3 y] - For the PSP from the class y, an offset duration. 

% Modulation targets
% Modulation can be affected by the genes of all three participants: post, pre_ion and pre_mod
% [8 1 x y] - Post_connection_pattern for modulation of the synapse x<-y  
% [8 2 x y] - Pre_connection_pattern for modulation of the synapse x<-y
% [8 3 x y] - Post_affinity_pattern for modulation of the synapse x<-y
% [8 4 x y] - Pre_affinity_pattern for modulation of the synapse x<-y
% [9 1 x y 1/2] - A singular modulaton effect. The fifth place: 1 - postsynapse (any subclass, either x or y, modulatory y is forbidden), 2 presynapse (mod). All for synapse x<-y
% [9 2 x y 1/2] - Max number of triggers. The fifth place: 1 - postsynapse (any subclass, either x or y, modulatory y is forbidden), 2 presynapse (mod)
% [9 3 x y 1/2] - Delay. The fifth place: 1 - postsynapse (any subclass, either x or y, modulatory y is forbidden), 2 presynapse (mod)
% [9 4 x y 1/2] - Onset duration. The fifth place: 1 - postsynapse (any subclass, either x or y, modulatory y is forbidden), 2 presynapse (mod)
% [9 5 x y 1/2] - Plateau duration. The fifth place: 1 - postsynapse (any subclass, either x or y, modulatory y is forbidden), 2 presynapse (mod)
% [9 6 x y 1/2] - Offset duration. The fifth place: 1 - postsynapse (any subclass, either x or y, modulatory y is forbidden), 2 presynapse (mod)
if Genotype_legal==1
    disp('the genotype is legal')
else error('STOP AND CHECK THE LEGALITY!')
end



end