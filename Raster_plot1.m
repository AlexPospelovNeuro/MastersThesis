



function Raster_plot1(ax,Inst,Raster,Desired_output,weights,text)





%figure
    Raster_forplot=Raster;
    N_Input=sum(Inst.Source.Cells.Input(1,:));
    N_Output=sum(Inst.Source.Cells.Output(1,:));
    N_Ion=sum(Inst.Source.Cells.Ion(1,:));
    N_Mod=sum(Inst.Source.Cells.Mod(1,:));
    cell_inclass=[];
    for cl=1:size(Inst.N_cells_vector,2)
        cell_inclass=[cell_inclass 1:Inst.N_cells_vector(cl)];
    end
    for cell=1:Inst.N_cells_total
        Raster_forplot(cell,Raster(cell,:)==1)=Inst.N_cells_total+1-cell;
        Raster_forplot(cell,Raster(cell,:)==0)=NaN;
        if cell<=N_Input
            %plot(Raster_forplot(cell,:),['*' Col{class_col(cell)}]);
            plot(ax,Raster_forplot(cell,:),'*k');
            hold on
        elseif cell<=N_Output+N_Input
            
            for t=1:size(weights{1}{cell-N_Input},2)
                col=1-weights{1}{cell-N_Input}(t)/5;
                rectangle(ax,'Position',[t Inst.N_cells_total-cell+0.5 1 1],'facecolor',[col col col],'edgecolor',[col col col])
            end
            hold on
            plot(ax,Raster_forplot(cell,:),'*b');
            %%%%%%%%% A tentative, only for this specific test
            Desired_output{1}{cell-N_Input}(Desired_output{1}{cell-N_Input}==1)=Inst.N_cells_total+1.25-cell;      
            Desired_output{1}{cell-N_Input}(Desired_output{1}{cell-N_Input}==0)=NaN;
            plot(ax,Desired_output{1}{cell-N_Input},'.r')
            %%%%%%%%%
        elseif cell<=N_Ion+N_Output+N_Input
            plot(ax,Raster_forplot(cell,:),'*g');
            hold on
        else
            plot(ax,Raster_forplot(cell,:),'og');
            hold on
        end
        line(ax,[5 5+1+Inst.AbsRefract{Inst.Cellclass(cell)}(cell_inclass(cell))],[Inst.N_cells_total+1-cell Inst.N_cells_total+1-cell],'color','k','linewidth',3) % The visual representation of the absolute refracterity length
        line(ax,[2.5 3.5],[Inst.N_cells_total+1-cell Inst.N_cells_total+1-cell],'color','k','linewidth',0.5)
        if Inst.Thresholds{Inst.Cellclass(cell)}(cell_inclass(cell))>=0
            line(ax,[3 3],[Inst.N_cells_total+1-cell Inst.N_cells_total+1-cell+Inst.Thresholds{Inst.Cellclass(cell)}(cell_inclass(cell))/10],'color','g','linewidth',3)
        else
            line(ax,[3 3],[Inst.N_cells_total+1-cell Inst.N_cells_total+1-cell+Inst.Thresholds{Inst.Cellclass(cell)}(cell_inclass(cell))/10],'color','r','linewidth',3)
        end
        line(ax,[6 6],[Inst.N_cells_total+1-cell-Inst.ThreshNoise{Inst.Cellclass(cell)}(cell_inclass(cell))/10 Inst.N_cells_total+1-cell+Inst.ThreshNoise{Inst.Cellclass(cell)}(cell_inclass(cell))/10],'color','b','linewidth',3)
        
    end
    for Neur=1:size(Inst.N_cells_vector_incremental_pointer,2)
        line(ax,[0 size(Raster_forplot,2)],[Inst.N_cells_total-Inst.N_cells_vector_incremental_pointer(Neur)+1.5 Inst.N_cells_total-Inst.N_cells_vector_incremental_pointer(Neur)+1.5],'Linestyle','--','Color','Black')
    end
    axis([-1 size(Raster,2)+10 -1 Inst.N_cells_total+2])
    
    title(text);
    
end