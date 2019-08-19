function [command,structure] = PRM_bender(command,structure)

% Prepare data for job 1 PRM and job 3 POS
% Bend the space from "atomic positions" to "atomic displacements"
% Set command.cartesian_vectors, an offset for job 3

if strcmp(command.space,'ub')
    map=structure.pos; command.cartesian_vectors=0;

else
    if command.nrun>1
        command.cartesian_vectors=ceil(norm(sum(sum(abs(structure.cartesian),4),3)./length(command.at).*structure.cell_parameters)/command.PRM_resolution);
        shi=zeros(1,3,length(command.at),command.nrun); shi(1,:,:,:)=0.5-structure.pos(1,:,:,:);
%         shi=zeros(1,3,length(command.at),command.nrun); shi(1,:,:,:)=0.5-repmat(structure.pos(1,:,:,1),[1,1,1,command.nrun]);
        if strcmp(command.space,'real')
            map=structure.pos+repmat(shi,[structure.length 1 1 1])+repmat(structure.cartesian,[structure.length 1 1 1]);
        else
            map=structure.pos+repmat(shi,[structure.length 1 1 1])+repmat(structure.cartesian,[structure.length 1 1 command.nrun]);
        end
        
    else
        command.cartesian_vectors=ceil(norm(sum(abs(structure.cartesian),3)./length(command.at).*structure.cell_parameters)/command.PRM_resolution);
        shi=zeros(1,3,length(command.at)); shi(1,:,:)=0.5-structure.pos(1,:,:);
        map=structure.pos+repmat(shi,[structure.length 1 1])+repmat(structure.cartesian,[structure.length 1 1]);
    end; clear shi;
end

% Rename
structure.displacements=map;

end

