function [] = PRM_job_PSD_print(command,job_PSD)

% Print outcome of job PSD

if command.print==1
    
    if strcmp(command.print_n,'')
        filename=sprintf(['PowerSpectrum',strrep(command.pathway,'dat\XDATCAR',''),'.txt']);
    else
        filename=sprintf(['PowerSpectrum',strrep(command.pathway,'dat\XDATCAR',''),'_',command.print_n,'.txt']);
    end
    fid=fopen(filename,'w');
    
    if strcmp(command.unit,'cm') || strcmp(command.unit,'cm-1')
        [~,encut]=min(abs(job_PSD.ve_en-5000)); job_PSD.ve_en=job_PSD.ve_en(1:encut); job_PSD.ve_ps=job_PSD.ve_ps(1:encut,:);
        en_label='Energy(cm-1)';
    else
        [~,encut]=min(abs(job_PSD.ve_en-600)); job_PSD.ve_en=job_PSD.ve_en(1:encut); job_PSD.ve_ps=job_PSD.ve_ps(1:encut,:);
        en_label='Energy(meV)';
    end
    
    if strcmp(command.PSD_partial,'elt')
        % Get atom labels
        command.elt_tmp=cell(length(command.elt),1);
        for k=1:length(command.elt); command.elt_tmp{k,1}=command.elt{1,k}{1,1}; end; clear k;
        [~,idx]=unique(command.elt_tmp,'first'); command.elt_rn=command.elt_tmp(sort(idx));
        %
        nelt=length(command.elt_rn); len=length(job_PSD.ve_en);
        fprintf(fid,'%12s   %12s',en_label,'PSiso(a.u.)');
        for k=1:nelt
            fprintf(fid,'   %12s',command.elt_rn{k,1}); 
        end; clear k;
        fprintf(fid,'\n');
        
        for h=1:len
            fprintf(fid,' %12.8e',[job_PSD.ve_en(h,1)';job_PSD.ve_ps(h,:)']);
            fprintf(fid,'\n');
        end; clear h;
        
    elseif strcmp(command.PSD_partial,'dir')
        fprintf(fid,'%12s  %12s %12s   %12s   %12s\r\n',en_label,'PSiso(a.u.)','PSx(a.u.)','PSy(a.u.)','PSz(a.u.)');
        fprintf(fid,'%12.8e %12.8e %12.8e %12.8e %12.8e\r\n',[job_PSD.ve_en';job_PSD.ve_ps']);

    else
        fprintf(fid,'%12s  %12s\r\n',en_label,'PSiso(a.u.)');
        fprintf(fid,'%12.8e %12.8e\r\n',[job_PSD.ve_en';job_PSD.ve_ps(:,1)']);
    end
    
    fclose(fid);
    clear filename fid;
    
else
end

end

