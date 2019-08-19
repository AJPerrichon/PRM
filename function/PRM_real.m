function [structure] = PRM_real(command,structure)

% command.space_grid for command.space='real': define a 3D high-symmetry grid format
% [Ox,Oy,Oz;V1x,V1y,V1z;V2x,V2y,V2z;V3x,V3y,V3z;Ox',Oy',Oz';V1x',... as a
% Oa(x,y,z) node and 3 vectors V1a(x,y,z), V2a(x,y,z), V3a(x,y,z) this can
% be extended with additional high-sym. grids with any number of 12
% parametes Ob(x,y,z), V1b(x,y,z), etc. If the structural model is best
% described in 2D (freedom of position in one of the cartesian direction),
% one can replace Ox/Oy/Oz by the '555' tag to use the SGS position
% instead. This 'trick' works only for a single-trajectory (command.nrun=1). For
% multi-trajectory or more complex cases, one should use the tag
% command.space='poscar' instead and tweak the POSCAR_* file accordingly.

% Define 3D grid (linear combination)
nsub=length(command.space_grid)/4; grid_lim=10;
grid_true=zeros(power(2*grid_lim+1,3),3,nsub,length(command.at));
command.space_grid_m=repmat(command.space_grid,1,1,length(command.at));

for h=1:length(command.at)
    for v=1:nsub
        
        % check 555 condition, see header
        if command.nrun>1 && sum(sum(command.space_grid==555,2))>0; error('Check /function/PRM_real header regarding the "555" tag'); else end
        if command.space_grid_m(1+4*(v-1),1,h)==555; command.space_grid_m(1+4*(v-1),1,h)=structure.pos(1,1,h); else end
        if command.space_grid_m(1+4*(v-1),2,h)==555; command.space_grid_m(1+4*(v-1),2,h)=structure.pos(1,2,h); else end
        if command.space_grid_m(1+4*(v-1),3,h)==555; command.space_grid_m(1+4*(v-1),3,h)=structure.pos(1,3,h); else end
        
        % build grid
        clear grid_xyz; 
        grid_xyz=zeros(2*grid_lim+1,2*grid_lim+1,2*grid_lim+1,3);
        for x=-grid_lim:1:grid_lim
            for y=-grid_lim:1:grid_lim
                for z=-grid_lim:1:grid_lim
                    grid_xyz(x+grid_lim+1,y+grid_lim+1,z+grid_lim+1,1)=command.space_grid_m(1+4*(v-1),1,h)...
                        +x*command.space_grid_m(2+4*(v-1),1,h)+y*command.space_grid_m(3+4*(v-1),1,h)+z*command.space_grid_m(4+4*(v-1),1,h);
                    grid_xyz(x+grid_lim+1,y+grid_lim+1,z+grid_lim+1,2)=command.space_grid_m(1+4*(v-1),2,h)...
                        +x*command.space_grid_m(2+4*(v-1),2,h)+y*command.space_grid_m(3+4*(v-1),2,h)+z*command.space_grid_m(4+4*(v-1),2,h);
                    grid_xyz(x+grid_lim+1,y+grid_lim+1,z+grid_lim+1,3)=command.space_grid_m(1+4*(v-1),3,h)...
                        +x*command.space_grid_m(2+4*(v-1),3,h)+y*command.space_grid_m(3+4*(v-1),3,h)+z*command.space_grid_m(4+4*(v-1),3,h);
                end
            end
        end; clear x y z;
        
        grid_true(:,1,v,h)=reshape(grid_xyz(:,:,:,1),power(2*grid_lim+1,3),1);
        grid_true(:,2,v,h)=reshape(grid_xyz(:,:,:,2),power(2*grid_lim+1,3),1);
        grid_true(:,3,v,h)=reshape(grid_xyz(:,:,:,3),power(2*grid_lim+1,3),1);
        
    end; clear v;
end; clear h;

grid_true=permute(grid_true,[2 1 3 4]);
grid_true=reshape(grid_true,3,length(grid_true)*nsub,length(command.at),1);
grid_true=permute(grid_true,[2 1 3]);

if command.nrun>1
    rea=zeros(1,3,length(command.at),command.nrun);
    for k=1:command.nrun
        for h=1:length(command.at)
            clear dl; dl=zeros(length(grid_true),1);
            for x=1:length(grid_true)
                dt=structure.pos(1,:,h,k)-grid_true(x,:,h);
                dl(x)=sqrt(power(dt(1),2)+power(dt(2),2)+power(dt(3),2));
            end; clear x;
            [~,tm]=min(dl(:));
            rea(1,:,h,k)=structure.pos(1,:,h,k)-grid_true(tm,:,h);
        end; clear h;
    end; clear k;
else
    rea=zeros(1,3,length(command.at));
    for h=1:length(command.at)
        clear dl; dl=zeros(length(grid_true),1);
        for x=1:length(grid_true)
            dt=structure.pos(1,:,h)-grid_true(x,:,h);
            dl(x)=sqrt(power(dt(1),2)+power(dt(2),2)+power(dt(3),2));
        end; clear x;
        [~,tm]=min(dl(:));
        rea(1,:,h)=structure.pos(1,:,h)-grid_true(tm,:,h);
    end
end

structure.cartesian=rea;

end

