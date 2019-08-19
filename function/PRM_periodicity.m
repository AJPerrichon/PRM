function [structure] = PRM_periodicity(command,structure)

% Correct data from periodic boundary folding; can't be parallelized (step
% by step domino effect). Older versions of this function were much faster,
% but could not manage multiple diffusion events in the same direction,
% which is an issue e.g. for H at high temperature. The only proper
% solution is through a while loop, which is particularly slow. The
% numerical value of |d| is arbitrarily set to 0.5 (in r.l.u.), which
% basically means an atom needs to move by a distance of 0.5 * (super)cell
% parameter, typically by crossing the cell boundary and appearing on the
% other side of the cell due to periodic conditions. Technically -0.8/0.8
% should be enough (which corresponds to a distance of 0.2 * (super)cell
% parameter.

if command.nrun>1
    for k=2:structure.length
        yL=(structure.pos(k,:,:,:)-structure.pos(k-1,:,:,:))<-0.5;      % cross the 1 rlu boundary
        zL=(structure.pos(k,:,:,:)-structure.pos(k-1,:,:,:))>0.5;       % cross the 0 rlu boundary
        if sum(sum(sum(yL,4),3),2)>0; y=1; else y=0; end                % set flag
        if sum(sum(sum(zL,4),3),2)>0; z=1; else z=0; end                % set flag
        while y==1
            structure.pos(k,:,:,:)=structure.pos(k,:,:,:)+yL;           % logical index used as double
            yL=(structure.pos(k,:,:,:)-structure.pos(k-1,:,:,:))<-0.5;
            if sum(sum(sum(yL,4),3),2)==0; y=0; else end
        end
        while z==1
            structure.pos(k,:,:,:)=structure.pos(k,:,:,:)-zL;           % logical index used as double
            zL=(structure.pos(k,:,:,:)-structure.pos(k-1,:,:,:))>0.5;
            if sum(sum(sum(zL,4),3),2)==0; z=0; else end
        end
    end
else
    for k=2:structure.length
        yL=(structure.pos(k,:,:)-structure.pos(k-1,:,:))<-0.5;          % cross the 1 rlu boundary
        zL=(structure.pos(k,:,:)-structure.pos(k-1,:,:))>0.5;           % cross the 0 rlu boundary
        if sum(sum(yL,3),2)>0; y=1; else y=0; end                       % set flag
        if sum(sum(zL,3),2)>0; z=1; else z=0; end                       % set flag
        while y==1
            structure.pos(k,:,:)=structure.pos(k,:,:)+yL;               % logical index used as double
            yL=(structure.pos(k,:,:)-structure.pos(k-1,:,:))<-0.5;
            if sum(sum(yL,3),2)==0; y=0; else end
        end
        while z==1
            structure.pos(k,:,:)=structure.pos(k,:,:)-zL;               % logical index used as double
            zL=(structure.pos(k,:,:)-structure.pos(k-1,:,:))>0.5;
            if sum(sum(zL,3),2)==0; z=0; else end
        end
    end
end

end

