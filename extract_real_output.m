% A small function to extract a "real_output" output from the raster (for
% the fitness calculation)


function Real_output=extract_real_output(Raster,Inst)
        class_pointer=1;
        for out_cl=1:size(Inst.Source.Cells.Output,2) % for each output subclass
            Real_output{out_cl}=Raster(sum(Inst.Source.Cells.Input(1,:))+class_pointer:sum(Inst.Source.Cells.Input(1,:))+class_pointer-1+Inst.Source.Cells.Output(1,out_cl),:);
            class_pointer=class_pointer+Inst.Source.Cells.Output(1,out_cl);
        end
end