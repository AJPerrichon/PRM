function [job_Kin] = PRM_job_Kin(command,structure)

% Calculate effective temperature of a set of selected atoms through their
% kinetic energy (velocities being calculated from trajectories)

% Constants
kb_eVK=8.6173324*10^-5;     % Boltzmann constant in eV/K
N=6.022140857*10^23;        % Avogadro constant in 1/mol
J_to_eV=6.241509*10^18;     % Joule to eV conversion

% Get kinetic energy in eV
if command.nrun>1
    ki_void=zeros(structure.length-2,length(command.at),command.nrun);
    for n=1:command.nrun
        for h=1:length(command.at)
            if strcmp(command.elt{1,h},'H'); m=1.00794; else end % mass in atomic unit
            if strcmp(command.elt{1,h},'O'); m=15.9994; else end
            if strcmp(command.elt{1,h},'Zr'); m=91.224; else end
            if strcmp(command.elt{1,h},'In'); m=114.818; else end
            if strcmp(command.elt{1,h},'Ba'); m=137.327; else end
            for k=1:structure.length-2
                ve=(structure.pos(k+2,:,h,n)-structure.pos(k+1,:,h,n))./command.time_step; % velocity in r.l.u./fs
                ki_void(k,h,n)=(power(ve(1)*structure.cell_parameters(1)*(10^5),2)...
                    +power(ve(2)*structure.cell_parameters(2)*(10^5),2)...
                    +power(ve(3)*structure.cell_parameters(3)*(10^5),2))...
                    *(1/2)*(m*(10^-3)/N)*J_to_eV;
            end; clear k;
        end; clear h m;
    end; clear n;
    ki_void=mean(ki_void,3);
else
    ki_void=zeros(structure.length-2,length(command.at));
    for h=1:length(command.at)
        if strcmp(command.elt{1,h},'H'); m=1.00794; else end % mass in atomic unit
        if strcmp(command.elt{1,h},'O'); m=15.9994; else end
        if strcmp(command.elt{1,h},'Zr'); m=91.224; else end
        if strcmp(command.elt{1,h},'In'); m=114.818; else end
        if strcmp(command.elt{1,h},'Ba'); m=137.327; else end
        for k=1:structure.length-2
            ve=(structure.pos(k+2,:,h)-structure.pos(k+1,:,h))./command.time_step; % velocity in r.l.u./fs
            ki_void(k,h)=(power(ve(1)*structure.cell_parameters(1)*(10^5),2)...
                +power(ve(2)*structure.cell_parameters(2)*(10^5),2)...
                +power(ve(3)*structure.cell_parameters(3)*(10^5),2))...
                *(1/2)*(m*(10^-3)/N)*J_to_eV;
        end; clear k;
    end; clear h m;
end

% Conversion in K
ki_T=mean(ki_void,2)*(2/3)/kb_eVK;
ki_t=command.time_step:command.time_step:(structure.length-2)*command.time_step;

% Rename
job_Kin.ki_t=ki_t;
job_Kin.ki_T=ki_T;

end

