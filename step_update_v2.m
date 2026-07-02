% The function that updates the raster by one timestep. 







function [Upd,Runtime]=step_update_v2(Slice,Runtime,Inst,t)
    Upd(1:Inst.N_cells_total)=Slice(:,end);
    for cell=Inst.N_input_cells+1:Inst.N_cells_total % for each neuron in the circuit except the inputs; Inputs are exluded by this routine 
        %[Raster, Runtime]=single_cell_update(Raster,Runtime,Inst,t,cell);
        %[Raster, Runtime]=single_cell_update_v3(Raster,Runtime,Inst,t,cell);
        [Upd(cell), Runtime]=single_cell_update_v3_1(Slice,Runtime,Inst,cell);
    end

end


