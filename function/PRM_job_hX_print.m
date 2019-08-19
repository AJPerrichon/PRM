function [] = PRM_job_hX_print(command,job_hX)

if command.print==1
    % Print output of h(MSD)
    filename=sprintf(['MSD',strrep(command.pathway,'dat\XDATCAR',''),'_',command.print_n,'.txt']);
    fid=fopen(filename,'w');
    fprintf(fid,'%12s %12s\r\n','distance (A)','MSD');
    fprintf(fid,'%12.8e %12.8e\r\n',[job_hX.od_vec;job_hX.od_void]);
    fclose(fid);
    clear filename fid;
    
    % Print outcome of h(POS)
    filename=sprintf(['POS',strrep(command.pathway,'dat\XDATCAR',''),'_rot%g','_',command.print_n,'.txt'],command.angle);
    fid=fopen(filename,'w');
    if command.sym==1 && strcmp(command.sym_op,'m-3m')
        if rem(round(command.supercell(1)*1000)/1000,1)==0 || rem(round(command.supercell(2)*100)/100,1)==0
            dir_n='x=y=z';
            fprintf(fid,'%12s %12s %12s\r\n','distance (A)',dir_n);
            fprintf(fid,'%12.8e %12.8e %12.8e\r\n',[job_hX.id_vec_x;job_hX.id_x]);
        else
            dir_n='x=y=z'; dir_m='m'; 
            fprintf(fid,'%12s %12s %12s\r\n','distance (A)',dir_n,dir_m');
            fprintf(fid,'%12.8e %12.8e %12.8e\r\n',[job_hX.id_vec_x;job_hX.id_x;job_hX.id_z]);
        end
    elseif command.sym==1 && strcmp(command.sym_op,'4/mmm')
        if rem(round(command.supercell(1)*1000)/1000,1)==0 || rem(round(command.supercell(2)*100)/100,1)==0;...
                dir_n='x=y'; else dir_n='m'; end
        fprintf(fid,'%12s %12s %12s\r\n','distance (A)',dir_n,'Z');
        fprintf(fid,'%12.8e %12.8e %12.8e\r\n',[job_hX.id_vec_x;job_hX.id_x;job_hX.id_z]);
    else
        if rem(round(command.supercell(1)*1000)/1000,1)==0 || rem(round(command.supercell(2)*100)/100,1)==0;...
                dir_n='x'; dir_m='y'; else dir_n='m'; dir_m='m"'; end
        fprintf(fid,'%12s %12s %12s %12s\r\n','distance (A)',dir_n,dir_m,'Z');
        fprintf(fid,'%12.8e %12.8e %12.8e %12.8e\r\n',[job_hX.id_vec_x;job_hX.id_x;job_hX.id_y;job_hX.id_z]);
    end
    fclose(fid);
    clear filename fid;
    
else
end


end

