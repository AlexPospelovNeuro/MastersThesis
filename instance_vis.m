% A function to visualize the circuit (the unfolded instance) and display
% as much information from it as possible.
% INPUT - Instance data structure
% OUTPUT - figure with the graphic representation.

% TO ADD: Handle as an optional output, possibly settings to adjust the
% representation


function instance_vis(ax,Inst)

% Rearrangement procedures: it is handy to keep "reserved" subclasses (inputs and outputs) together for the simulation purposes, but is seems to be more appealing to push the outputs to the end of the list for visualization

Scl_forplot_order=[1:size(Inst.Source.Cells.Input,2) (size(Inst.Source.Cells.Input,2)+size(Inst.Source.Cells.Output,2)+1):Inst.N_subclasses (size(Inst.Source.Cells.Input,2)+1):(size(Inst.Source.Cells.Input,2)+size(Inst.Source.Cells.Output,2))];
Scl_forplot_NONMOD_order=[1:size(Inst.Source.Cells.Input,2) (size(Inst.Source.Cells.Input,2)+size(Inst.Source.Cells.Output,2)+1):Inst.N_nonmod_subclasses (size(Inst.Source.Cells.Input,2)+1):(size(Inst.Source.Cells.Input,2)+size(Inst.Source.Cells.Output,2))];
forplot_N_cells_vector=Inst.N_cells_vector(Scl_forplot_order);
forplot_N_cells_vector_incremental=cumsum(forplot_N_cells_vector);
forplot_N_cells_NONMOD_vector=Inst.N_cells_vector(Scl_forplot_NONMOD_order);
forplot_N_cells_NONMOD_vector_incremental=cumsum(forplot_N_cells_NONMOD_vector);
forplot_class=repmat(1,[1 size(Inst.Source.Cells.Input,2)]); % Inputs
forplot_class=[forplot_class repmat(2,[1 size(Inst.Source.Cells.Ion,2)])]; % Ionotropic
forplot_class_NONMOD=forplot_class;
forplot_class=[forplot_class repmat(3,[1 size(Inst.Source.Cells.Mod,2)])]; % Modulatory
forplot_class=[forplot_class repmat(4,[1 size(Inst.Source.Cells.Output,2)])]; % Outputs
forplot_class_NONMOD=[forplot_class_NONMOD repmat(4,[1 size(Inst.Source.Cells.Output,2)])];

forplot_cell=repmat(1,[1 sum(Inst.Source.Cells.Input(1,:))]); % Inputs
forplot_cell=[forplot_cell repmat(2,[1 sum(Inst.Source.Cells.Ion(1,:))])];  % Ionotropic
forplot_cell_NONMOD=forplot_cell;
forplot_cell=[forplot_cell repmat(3,[1 sum(Inst.Source.Cells.Mod(1,:))])]; % Modulatory
forplot_cell=[forplot_cell repmat(4,[1 sum(Inst.Source.Cells.Output(1,:))])]; % Outputs
forplot_cell_NONMOD=[forplot_cell_NONMOD repmat(4,[1 sum(Inst.Source.Cells.Output(1,:))])];

Orig_index=1:forplot_N_cells_vector_incremental(end);
New_index=[Orig_index(1:sum(Inst.Source.Cells.Input(1,:))) Orig_index(sum(Inst.Source.Cells.Input(1,:))+sum(Inst.Source.Cells.Output(1,:))+1:Orig_index(end)) Orig_index(sum(Inst.Source.Cells.Input(1,:))+1:sum(Inst.Source.Cells.Input(1,:))+sum(Inst.Source.Cells.Output(1,:)))];
Orig_index_NONMOD=1:forplot_N_cells_NONMOD_vector_incremental(end);
New_index_NONMOD=[Orig_index_NONMOD(1:sum(Inst.Source.Cells.Input(1,:))) Orig_index_NONMOD(sum(Inst.Source.Cells.Input(1,:))+sum(Inst.Source.Cells.Output(1,:))+1:sum(Inst.Source.Cells.Input(1,:))+sum(Inst.Source.Cells.Output(1,:))+sum(Inst.Source.Cells.Ion(1,:))) Orig_index_NONMOD(sum(Inst.Source.Cells.Input(1,:))+1:sum(Inst.Source.Cells.Input(1,:))+sum(Inst.Source.Cells.Output(1,:)))];


%f=figure;
%ax=axes(f,'Position',[0 0 1 1]);
Hat=0.1;
ax.XTick=[];
ax.YTick=[];
axis([-Hat 1 -Hat 1])
hold on
line(ax,[0 0],[-Hat 1],'Color','black')
line(ax,[-Hat 1],[1-Hat 1-Hat],'Color','black')
StepPost=1/Inst.N_cells_total;
StepPre=1/Inst.N_cells_nonmod;
Colors={[0 0 0],[0 1 0],[0 0 1],[0.5 0.5 0.5]};

%The subclass grid
for postcl=1:Inst.N_subclasses % Horizontal
    line(ax,[-Hat 1],[1-Hat-forplot_N_cells_vector_incremental(postcl)*StepPost 1-Hat-forplot_N_cells_vector_incremental(postcl)*StepPost],'Color','red')
end
for precl=1:Inst.N_nonmod_subclasses % Vertical
    line(ax,[forplot_N_cells_NONMOD_vector_incremental(precl)*StepPre forplot_N_cells_NONMOD_vector_incremental(precl)*StepPre],[-Hat 1],'Color','red')
end

% Cell icons of a proper color
for n_post=1:size(forplot_cell,2) % Postsynaptic
    Cell_post_handle(n_post) = nsidedpoly(3, 'Center', [-Hat/3 1-Hat-n_post*StepPost+StepPost/2], 'Radius', Hat/10);
    plot(ax,Cell_post_handle(n_post),'FaceColor',Colors{forplot_cell(n_post)});
    hold on
end
for n_pre=1:size(forplot_cell_NONMOD,2) % Presynaptic
    Cell_pre_handle(n_pre) = nsidedpoly(3, 'Center', [(n_pre-1)*StepPre+StepPre/2 1-2*Hat/3], 'Radius', Hat/10);
    plot(ax,Cell_pre_handle(n_pre),'FaceColor',Colors{forplot_cell_NONMOD(n_pre)});
    hold on
end

% The connections 

for post_cell=1:size(forplot_cell,2)
    for pre_cell=1:size(forplot_cell_NONMOD,2)
        if Inst.Connection_matrix(New_index(post_cell),New_index_NONMOD(pre_cell))==1
            Power=Inst.Synapses.Powers{Inst.Cellclass(New_index(post_cell)),Inst.Cellclass(New_index_NONMOD(pre_cell))}(Inst.Synapses.IDs((Inst.Synapses.IDs(:,4)==New_index(post_cell))&(Inst.Synapses.IDs(:,5)==New_index_NONMOD(pre_cell)),6));
            Contact_handle(post_cell,pre_cell)=nsidedpoly(12, 'Center', [(pre_cell-1)*StepPre+StepPre/2 1-Hat-post_cell*StepPost+StepPost/2], 'Radius', log10(2+abs(Power))*Hat/10);
            if Power>=0
                plot(ax,Contact_handle(post_cell,pre_cell),'FaceColor','green','edgecolor','green');
                text((pre_cell-1)*StepPre+StepPre/2,1-Hat-post_cell*StepPost+StepPost/2,num2str(Inst.Synapses.Delays{Inst.Cellclass(New_index(post_cell)),Inst.Cellclass(New_index_NONMOD(pre_cell))}(Inst.Synapses.IDs((Inst.Synapses.IDs(:,4)==New_index(post_cell))&(Inst.Synapses.IDs(:,5)==New_index_NONMOD(pre_cell)),6))))
            else
                plot(ax,Contact_handle(post_cell,pre_cell),'FaceColor','red','edgecolor','red');
                text((pre_cell-1)*StepPre+StepPre/2,1-Hat-post_cell*StepPost+StepPost/2,num2str(Inst.Synapses.Delays{Inst.Cellclass(New_index(post_cell)),Inst.Cellclass(New_index_NONMOD(pre_cell))}(Inst.Synapses.IDs((Inst.Synapses.IDs(:,4)==New_index(post_cell))&(Inst.Synapses.IDs(:,5)==New_index_NONMOD(pre_cell)),6))))
            
            end
        end
    end
end


end



