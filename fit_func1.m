


function ft=fit_func1(Desired_output1,Real_output1,weights,Inst,Raster)
ft1=[];
for cl=1:size(Real_output1,2) % for each output subclass
    % Desired_output1{cl}
    % size(Desired_output1{cl})
    Desired_output{cl}=[Desired_output1{cl}(:,weights{cl}>0) repmat([0 1],1,size(Desired_output1{cl},1))]; % Add [0 1 to the each vector so that the variance is never 0]
    % Real_output1{cl}
    % size(Real_output1{cl})
    Real_output{cl}=[Real_output1{cl}(:,weights{cl}>0) repmat([0 1],1,size(Real_output1{cl},1))];
    
    for ce=1:size(Desired_output1{cl},1)
        ft1=[ft1 corr(Desired_output{cl}(ce,:)',Real_output{cl}(ce,:)','Type','Spearman')];
    end
end

%     if (sum(Desired_output)==0)&&(sum(Real_output)==0)
%         ft=1;
%     elseif (sum(Desired_output)==size(Desired_output,2))&&(sum(Real_output)==size(Real_output,2))
%         ft=1;
%     elseif (sum(Desired_output)==size(Desired_output,2))&&(sum(Real_output)==0)
%         ft=-1;
%     elseif (sum(Desired_output)==0)&&(sum(Real_output)==size(Real_output,2))
%         ft=-1;
%     elseif var(Real_output)==0
%         ft=-1;
%     else
        %ft=corr(Desired_output',Real_output','Type','Spearman');
%     end
    % A structural penalty: 0.0001 for each non-reserved neuron.
    ft=mean(ft1);
    N_reserved_cells=sum(Inst.Source.Cells.Input(1,:))+sum(Inst.Source.Cells.Output(1,:));
    st_pen=(Inst.N_cells_total-N_reserved_cells)/10000;
    % A metabolic penalty: A total fraction of the activity divided by 1000
    met_pen=(sum(sum(Raster(N_reserved_cells+1:end,:)))/(size(Raster,2)*(Inst.N_cells_total-N_reserved_cells+1)))/1000;
    
    ft=ft-st_pen-met_pen;
end


