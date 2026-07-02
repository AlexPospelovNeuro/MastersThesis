
% A separate function that makes the Runtime variable with all the trackers
% based on the Inst and shapes databases. Since it is fully deterministic,
% does not use the seed. Paired with unfold_v2 function. Can be called
% separately, for example for the resimulation. 


function Runtime=make_runtime(Inst,PSP_db,Plast_db)

Runtime.Thrsh=[];
Runtime.ThrshNoise=[];
Runtime.AbsRefrRef=[];
for post_cl=1:size(Inst.Synapses.PSPparam,1) 
    for pre_cl=1:size(Inst.Synapses.PSPparam,2) 
        for syn=1:size(Inst.Synapses.PSPparam{post_cl,pre_cl},1)
            Runtime.PSP{post_cl,pre_cl}{syn}=fliplr([PSP_db.PSP.in{Inst.Synapses.PSPparam{post_cl,pre_cl}(syn,1)} PSP_db.PSP.mid{Inst.Synapses.PSPparam{post_cl,pre_cl}(syn,2)} PSP_db.PSP.out{Inst.Synapses.PSPparam{post_cl,pre_cl}(syn,3)}]);
            Runtime.PSPlength{post_cl,pre_cl}{syn}=size(Runtime.PSP{post_cl,pre_cl}{syn},2);
            
        end
    end
    Runtime.Thrsh=[Runtime.Thrsh Inst.Thresholds{post_cl}]; % The thresholds in a long list
    Runtime.ThrshNoise=[Runtime.ThrshNoise Inst.ThreshNoise{post_cl}]; % The threshold noises in a long list
    Runtime.AbsRefrRef=[Runtime.AbsRefrRef Inst.AbsRefract{post_cl}]; % The absolute refracterity reference
end

for post_cl=1:size(Inst.Synapses.Plast,1) 
    for pre_cl=1:size(Inst.Synapses.Plast,2) 
        % The plasticity shapes and the track tables (Except Relative refracterity, type 1000)
        for Plst_type=1:size(Inst.Synapses.Plast{post_cl,pre_cl},2)
            for Plst_pat=1:size(Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type},2)
                for syn=1:size(Inst.Synapses.Powers{post_cl,pre_cl},2) % For each synapse between between an ordered pair of subclasses
                    Runtime.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{syn}=[zeros(1,Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{4}(1,syn)) Plast_db.plast.in{Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{4}(2,syn)} ...
                        Plast_db.plast.mid{Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{4}(3,syn)}  Plast_db.plast.out{Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{4}(4,syn)} ]; % A shape template (plateau=1)
                    Runtime.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{syn}=Runtime.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{syn}*Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{2}(1,syn); % A final shape (amplitude considered)
                    Runtime.PlastTrack{post_cl,pre_cl}{Plst_type}{Plst_pat}{syn}(1,1:size(Runtime.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{syn},2))=0;
                    Runtime.PlastN{post_cl,pre_cl}{Plst_type}{Plst_pat}{syn}=Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{3}(1,syn);
                   
                end
            end 
        end
        for mod_cl=1:size(Inst.Synapses.Mod,3)
            if ~isempty(Inst.Synapses.Mod{post_cl,pre_cl,mod_cl})
                for syn=1:size(Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{1},2)
                    Runtime.Mod{post_cl,pre_cl,mod_cl}{syn}=[zeros(1,Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{3}(1,syn)) Plast_db.plast.in{Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{3}(2,syn)}...
                        Plast_db.plast.mid{Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{3}(3,syn)} Plast_db.plast.out{Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{3}(4,syn)}];
                    Runtime.Mod{post_cl,pre_cl,mod_cl}{syn}=Runtime.Mod{post_cl,pre_cl,mod_cl}{syn}*Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{1}(1,syn);
                    Runtime.ModTrack{post_cl,pre_cl,mod_cl}{syn}(1, 1:size(Runtime.Mod{post_cl,pre_cl,mod_cl}{syn},2))=0;
                    Runtime.ModN{post_cl,pre_cl,mod_cl}{syn}=Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{2}(1,syn);
                end
            end
        end
    end 
    for Plst_pat=1:size(Inst.RelRefract{post_cl},2)
        for cell=1:Inst.N_cells_vector(post_cl) % For each cell of the postsynaptic subclass
            Runtime.RelRef{post_cl}{Plst_pat}{cell}=[zeros(1,Inst.RelRefract{post_cl}{Plst_pat}{4}(1,cell)) Plast_db.plast.in{Inst.RelRefract{post_cl}{Plst_pat}{4}(2,cell)} ...
                Plast_db.plast.mid{Inst.RelRefract{post_cl}{Plst_pat}{4}(3,cell)}  Plast_db.plast.out{Inst.RelRefract{post_cl}{Plst_pat}{4}(4,cell)} ]; % A shape template (plateau=1)
            Runtime.RelRef{post_cl}{Plst_pat}{cell}=Runtime.RelRef{post_cl}{Plst_pat}{cell}*Inst.RelRefract{post_cl}{Plst_pat}{2}(1,cell); % A final shape (amplitude considered)
            Runtime.RefTrack{post_cl}{Plst_pat}{cell}(1,1:size(Runtime.RelRef{post_cl}{Plst_pat}{cell},2))=0;            
            Runtime.RefN{post_cl}{Plst_pat}{cell}=Inst.RelRefract{post_cl}{Plst_pat}{3}(1,cell);
        end       
    end
end

Runtime.AbsRefCounter=zeros(1,Inst.N_cells_total);           
                
                

end