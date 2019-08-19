function [m] = PRM_mass(elt)

% Molar mass

if strcmp(elt,'H')
    m=1.00794;
    
elseif strcmp(elt,'D')
    m=2.0141;
    
elseif strcmp(elt,'O')
    m=15.9994;
    
elseif strcmp(elt,'Sc')
    m=44.955908;
    
elseif strcmp(elt,'Ni')
    m=58.6934;
    
elseif strcmp(elt,'Ti')
    m=47.867;
    
elseif strcmp(elt,'Y')
    m=88.90584;
    
elseif strcmp(elt,'Zr')
    m=91.224;
    
elseif strcmp(elt,'In')
    m=114.818;
    
elseif strcmp(elt,'Ba')
    m=137.327;
    
elseif strcmp(elt,'Nd')
    m=144.242;
    
else
end


end

