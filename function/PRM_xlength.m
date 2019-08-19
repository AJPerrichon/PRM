function [c] = PRM_xlength(elt,command)

% Neutron cross-sections
% Get scattering cross sections; http://www.ncnr.nist.gov/resources/n-lengths/

% type should be 'xc' for the cross-section, or 'b' for the scattering length

if strcmp(command.bweight,'n')
    c=1;
    
else
    
    if strcmp(elt,'H')
        if strcmp(command.bweight,'i')
            c=80.26;
        elseif strcmp(command.bweight,'c')
            c=1.7568;
        elseif strcmp(command.bweight,'t')
            c=82.02;
        else
        end

    elseif strcmp(elt,'D')
        if strcmp(command.bweight,'i')
            c=2.05;
        elseif strcmp(command.bweight,'c')
            c=5.592;
        elseif strcmp(command.bweight,'t')
            c=7.64;
        else
        end

    elseif strcmp(elt,'O')
        if strcmp(command.bweight,'i')
            c=0.0008;
        elseif strcmp(command.bweight,'c')
            c=4.232;
        elseif strcmp(command.bweight,'t')
            c=4.232;
        else
        end

    elseif strcmp(elt,'Sc')
        if strcmp(command.bweight,'i')
            c=4.5;
        elseif strcmp(command.bweight,'c')
            c=19;
        elseif strcmp(command.bweight,'t')
            c=23.5;
        else
        end
        
    elseif strcmp(elt,'Ni')
        if strcmp(command.bweight,'i')
            c=5.2;
        elseif strcmp(command.bweight,'c')
            c=13.3;
        elseif strcmp(command.bweight,'t')
            c=18.5;
        else
        end

    elseif strcmp(elt,'Ti')
        if strcmp(command.bweight,'i')
            c=2.87;
        elseif strcmp(command.bweight,'c')
            c=1.485;
        elseif strcmp(command.bweight,'t')
            c=4.35;
        else
        end

    elseif strcmp(elt,'Y')
        if strcmp(command.bweight,'i')
            c=0.15;
        elseif strcmp(command.bweight,'c')
            c=7.55;
        elseif strcmp(command.bweight,'t')
            c=7.7;
        else
        end

    elseif strcmp(elt,'Zr')
        if strcmp(command.bweight,'i')
            c=0.02;
        elseif strcmp(command.bweight,'c')
            c=6.44;
        elseif strcmp(command.bweight,'t')
            c=6.46;
        else
        end

    elseif strcmp(elt,'In')
        if strcmp(command.bweight,'i')
            c=0.54;
        elseif strcmp(command.bweight,'c')
            c=2.08;
        elseif strcmp(command.bweight,'t')
            c=2.62;
        else
        end

    elseif strcmp(elt,'Ba')
        if strcmp(command.bweight,'i')
            c=0.15;
        elseif strcmp(command.bweight,'c')
            c=3.23;
        elseif strcmp(command.bweight,'t')
            c=3.38;
        else
        end
    
    elseif strcmp(elt,'Nd')
        if strcmp(command.bweight,'i')
            c=9.2;
        elseif strcmp(command.bweight,'c')
            c=7.43;
        elseif strcmp(command.bweight,'t')
            c=16.6;
        else
        end

    else
    end
    
%     if strcmp(type,'b'); c=sqrt(c/4/pi*100); end % xc to b

end

end

