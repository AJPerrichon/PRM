function [job_hX] = PRM_job_hX(command,structure)

% Job: Histogram of mean square displacements h(MSD) and positions h(POS)

% set box size for h(MSD)
od_vec_max=sqrt(power(structure.cell_parameters(1),2)+...
    power(structure.cell_parameters(2),2)+...
    power(structure.cell_parameters(3),2));
od_vec=(0:command.PRM_resolution:od_vec_max);
od_void=zeros(1,length(od_vec)); od_msd=zeros(length(command.at),structure.length-1);

% set box size for h(POS)
id_vec_x=(0:command.PRM_resolution:max(structure.cell_parameters));
id_vec_y=(0:command.PRM_resolution:max(structure.cell_parameters));
id_vec_z=(0:command.PRM_resolution:max(structure.cell_parameters));
id_x=zeros(1,length(id_vec_x)); id_y=zeros(1,length(id_vec_y)); id_z=zeros(1,length(id_vec_z));

% job h(MSD)
if command.nrun>1
    for n=1:command.nrun
        for h=1:length(command.at)
            for k=2:structure.length
                od_msd(h,k,n)=sqrt(power(((structure.displacements(k,1,h,n)-0.5)*structure.cell_parameters(1)),2)+...
                    power(((structure.displacements(k,2,h,n)-0.5)*structure.cell_parameters(2)),2)+...
                    power(((structure.displacements(k,3,h,n)-0.5)*structure.cell_parameters(3)),2));
                [~,od_cd]=min(abs(od_msd(h,k,n)-od_vec));
                od_void(od_cd)=od_void(od_cd)+1;
            end
        end
    end; clear n h k od_cd od_msd;
else
    for h=1:length(command.at)
        for k=2:structure.length
            od_msd(h,k)=sqrt(power(((structure.displacements(k,1,h)-0.5)*structure.cell_parameters(1)),2)+...
                power(((structure.displacements(k,2,h)-0.5)*structure.cell_parameters(2)),2)+...
                power(((structure.displacements(k,3,h)-0.5)*structure.cell_parameters(3)),2));
            [~,od_cd]=min(abs(od_msd(h,k)-od_vec));
            od_void(od_cd)=od_void(od_cd)+1;
        end
    end; clear h k od_cd od_msd;
end

% job h(POS)
if command.nrun>1
    for n=1:command.nrun
        for h=1:length(command.at)
            for k=2:structure.length
                [~,id_cd_x]=min(abs(((abs(structure.displacements(k,1,h,n)-0.5))*structure.cell_parameters(1))-id_vec_x));
                [~,id_cd_y]=min(abs(((abs(structure.displacements(k,2,h,n)-0.5))*structure.cell_parameters(2))-id_vec_y));
                [~,id_cd_z]=min(abs(((abs(structure.displacements(k,3,h,n)-0.5))*structure.cell_parameters(3))-id_vec_z));
                id_cd=100+command.cartesian_vectors;
                if id_cd_y<id_cd && id_cd_z<id_cd; id_x(id_cd_x)=id_x(id_cd_x)+1; else end
                if id_cd_x<id_cd && id_cd_z<id_cd; id_y(id_cd_y)=id_y(id_cd_y)+1; else end
                if id_cd_x<id_cd && id_cd_y<id_cd; id_z(id_cd_z)=id_z(id_cd_z)+1; else end
            end
        end
    end; clear n h k id_cd_x id_cd_y id_cd_z;
else
    for h=1:length(command.at)
        for k=2:structure.length
            [~,id_cd_x]=min(abs(((abs(structure.displacements(k,1,h)-0.5))*structure.cell_parameters(1))-id_vec_x));
            [~,id_cd_y]=min(abs(((abs(structure.displacements(k,2,h)-0.5))*structure.cell_parameters(2))-id_vec_y));
            [~,id_cd_z]=min(abs(((abs(structure.displacements(k,3,h)-0.5))*structure.cell_parameters(3))-id_vec_z));
            id_cd=100+command.cartesian_vectors; % 100 % cutoff 10 boxes + 'command.cartesian_vectors' offset (about 0.3+ AA) for perpendicular directions
            if id_cd_y<id_cd && id_cd_z<id_cd; id_x(id_cd_x)=id_x(id_cd_x)+1; else end % non-consistent yet clearer
            if id_cd_x<id_cd && id_cd_z<id_cd; id_y(id_cd_y)=id_y(id_cd_y)+1; else end
            if id_cd_x<id_cd && id_cd_y<id_cd; id_z(id_cd_z)=id_z(id_cd_z)+1; else end
        end; clear k;
    end; clear h id_cd_x id_cd_y id_cd_z;
end
id_x(1)=id_x(1).*2; id_y(1)=id_y(1).*2; id_z(1)=id_z(1).*2; % correct double counting

% Deleting empty end data
od_n=floor(smooth(+(od_void'==0),11))==1; 
od_void(od_n)=[]; od_vec(od_n)=[]; 
clear od_n od_vec_max;
id_n=max(+(mean(cat(2,id_x(:),id_y(:),id_z(:)),2)~=0).*(1:1:length(id_x))');
id_x(id_n+10:end)=[]; 
id_y(id_n+10:end)=[]; 
id_z(id_n+10:end)=[];
id_vec_x(id_n+10:end)=[]; 
id_vec_y(id_n+10:end)=[]; 
id_vec_z(id_n+10:end)=[]; 
clear id_n;

% Normalisation
od_void=od_void./structure.total_count;
id_x=id_x./structure.total_count;
id_y=id_y./structure.total_count;
id_z=id_z./structure.total_count;

% Symmetry expansion h(POS)
if command.sym==1
    if strcmp(command.sym_op,'4/mmm')==1
        id_x=mean(cat(1,id_x,id_y),1); id_y=id_x;
    elseif strcmp(command.sym_op{1},'m-3m')==1
        id_x=mean(cat(1,id_x,id_y,id_z),1); id_y=id_x; id_z=id_x;
    else
    end
else
end

% Rename
job_hX.od_vec=od_vec;
job_hX.od_void=od_void;
job_hX.id_x=id_x;
job_hX.id_y=id_y;
job_hX.id_z=id_z;
job_hX.id_vec_x=id_vec_x;
job_hX.id_vec_y=id_vec_y;
job_hX.id_vec_z=id_vec_z;

end

