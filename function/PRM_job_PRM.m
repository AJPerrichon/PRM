function [command,job_PRM] = PRM_job_PRM(command,structure)

% Perform Positional Recurrence Map (PRM) job, which consists, for each
% atom, for each frame of the trajectory, in locating in which 'box' is the
% atom. The PRM cannot be fully parallelized, because it would require too
% much memory.

% set box size for PRM (PRM resolution)
bs_xvec=(0:command.PRM_resolution/structure.cell_parameters(1):1)';
bs_yvec=(0:command.PRM_resolution/structure.cell_parameters(2):1)';
bs_zvec=(0:command.PRM_resolution/structure.cell_parameters(3):1)';
bs_void=zeros(length(bs_xvec),length(bs_yvec),length(bs_zvec));

% job 1
if command.nrun>1
    [len,~,~,~]=size(structure.displacements);
    for n=1:command.nrun
        for h=1:length(command.at)
            for k=2:len
                [~,id_x]=min(abs(structure.displacements(k,1,h,n)-bs_xvec));
                [~,id_y]=min(abs(structure.displacements(k,2,h,n)-bs_yvec));
                [~,id_z]=min(abs(structure.displacements(k,3,h,n)-bs_zvec));
                bs_void(id_x,id_y,id_z)=bs_void(id_x,id_y,id_z)+1;
            end
        end
    end; clear n k h id_x id_y id_z len;
else
    [len,~,~]=size(structure.displacements);
    for h=1:length(command.at)
        for k=2:len
            [~,id_x]=min(abs(structure.displacements(k,1,h)-bs_xvec));
            [~,id_y]=min(abs(structure.displacements(k,2,h)-bs_yvec));
            [~,id_z]=min(abs(structure.displacements(k,3,h)-bs_zvec));
            bs_void(id_x,id_y,id_z)=bs_void(id_x,id_y,id_z)+1;
        end
    end; clear h k id_x id_y id_z len;
end

% Deleting empty planes to save memory for symmetry expansion
cd_x=smooth(+squeeze(mean(mean(bs_void(:,:,:),2),3)==0),5);
cd_x=floor(mean([cd_x,flip(cd_x,1)],2));
cd_y=smooth(+squeeze(mean(mean(bs_void(:,:,:),3),1)==0)',5);
cd_y=floor(mean([cd_y,flip(cd_y,1)],2));
cd_z=smooth(+squeeze(mean(mean(bs_void(:,:,:),1),2)==0),5);
cd_z=floor(mean([cd_z,flip(cd_z,1)],2));
if strcmp(command.sym_op,'4/mmm') || strcmp(command.sym_op,'4/m')
    cd_x=floor(mean([cd_x,cd_y],2)); cd_y=cd_x;
elseif strcmp(command.sym_op,'m-3m') || strcmp(command.sym_op,'m-3')
    cd_x=floor(mean([cd_x,cd_y,cd_z],2)); cd_y=cd_x; cd_z=cd_x;
end
bs_void(cd_x==1,:,:)=[]; bs_void(:,cd_y==1,:)=[]; bs_void(:,:,cd_z==1)=[];
bs_xvec(cd_x==1)=[]; bs_yvec(cd_y==1)=[]; bs_zvec(cd_z==1)=[];

% Normalisation
bs_void=bs_void./structure.total_count;
command.supercell(command.PRM_direction)=[];
command.supercell=fliplr(command.supercell);

% Rename
job_PRM.bs_void=bs_void;
job_PRM.xvec=bs_xvec;
job_PRM.yvec=bs_yvec;
job_PRM.zvec=bs_zvec;
job_PRM.xvec_length=length(job_PRM.xvec);
job_PRM.yvec_length=length(job_PRM.yvec);
job_PRM.zvec_length=length(job_PRM.zvec);

end

