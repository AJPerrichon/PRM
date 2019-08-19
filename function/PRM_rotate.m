function [command,structure]=PRM_rotate(command,structure)

% NOTE: End result of the transformation MUST BE orthogonal, as the
% trajectories are stored into orthogonal tensors. Which means basically
% that function is made when the SC is rotated by 45° with a transformation
% [1 1 0; -1 1 0; 0 0 1], and you want to project into the standard [1 0 0,
% 0 1 0; 0 0 1] setup.

% Adapt cell
rmat=[cosd(command.angle) sind(command.angle) 0;-sind(command.angle) cosd(command.angle) 0;0 0 1];
latb=structure.cell_parameters; latb(1)=latb(1)*sqrt(2); latb(2)=latb(2)*sqrt(2);

% Rotate pos vectors
posb=zeros(size(structure.pos));
if command.nrun>1
    for n=1:command.nrun
        for h=1:length(command.at)
            for k=1:structure.length
                rval=((rmat*(structure.pos(k,:,h,n)'-.5).*structure.cell_parameters')./latb')+.5;
                posb(k,:,h,n)=rval';
            end; clear k rval;
        end; clear h; 
    end; clear n;
else
    for h=1:length(command.at)
        for k=1:structure.length
            rval=((rmat*(structure.pos(k,:,h)'-.5).*structure.cell_parameters')./latb')+.5;
            posb(k,:,h)=rval';
        end; clear k rval;
    end; clear h; 
end
structure.pos=posb;

% Rotate rea vectors
if strcmp(command.space,'real')
    reab=zeros(size(structure.cartesian));
    if command.nrun>1
        for n=1:command.nrun
            for h=1:length(command.at)
                rval=(rmat*(structure.cartesian(1,:,h,n)').*structure.cell_parameters')./latb';
                reab(1,:,h,n)=rval';
            end; clear h rval; 
        end; clear n;
    else
        for h=1:length(command.at)
            rval=(rmat*(structure.cartesian(1,:,h)').*structure.cell_parameters')./latb';
            reab(1,:,h)=rval';
        end; clear h rval; 
    end
    structure.cartesian=reab;
else
end

% Rename
structure.cell_parameters=latb; 
command.supercell(1)=command.supercell(1)*sqrt(2); 
command.supercell(2)=command.supercell(2)*sqrt(2);

end

