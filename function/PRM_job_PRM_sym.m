function [job_PRM] = PRM_job_PRM_sym(command,job_PRM)

% Crystal symmetry defined in Hermann-Mauguin notations. Since the aim is
% to relate to diffraction data, symmetry is only defined for Laue classes
% (it might work for non-centrosymmetric symmetry, but then equivalent
% atoms should be selected carefully). And yet only for ortho/tetra/cubic
% (mmm, 4/m, 4/mmm, m-3, m-3m) due to how the code is build (filling matlab
% 3d tensors with orthogonal dimensions). Data should be 'oriented' in
% crystallographic format, as described in the International Tables (i.e.
% e.g. for tetragonal classes, the fourfold axis should be in the c
% direction). Symmetry is implemented as follow: each pixel of the bs_void
% var is considered as an atom on the highest multiplicity (x,y,z) position
% (e.g. m=48 for m-3m), hence we use the set of positions defined in the
% International Tables. Operations are performed on the 3D tensor: e.g.
% (xyz)>(xy-z) is basically flip(...,3). If your structure is not in
% crystallographic base (e.g. 4-axis in T is along b ..?), please use
% permute(bs_void,[3 1 2]) to set in right prior to the symmetry
% application, and permute(bs_sym,[2 3 1]) to rotate it back afterward (or
% you can also set a new 'tag' as (m4/mm) and code positions accordingly).
% There is a risk of RAM memory overload due to the size of bs_sym; if that
% happens, set "bs_sym4(:,:,:)" instead of "bs_sym(:,:,:,4)" to use HD
% instead of RAM.

if strcmp(command.sym_op,'mmm')==1
    
    bs_sym=zeros(length(job_PRM.xvec),length(job_PRM.yvec),length(job_PRM.zvec),8);
    bs_sym(:,:,:,1)=job_PRM.bs_void;                                              % (1) x,y,z
    bs_sym(:,:,:,2)=flip(flip(job_PRM.bs_void,1),2);                              % (2) -x,-y,z
    bs_sym(:,:,:,3)=flip(flip(job_PRM.bs_void,1),3);                              % (3) -x,y,-z
    bs_sym(:,:,:,4)=flip(flip(job_PRM.bs_void,2),3);                              % (4) x,-y,-z
    bs_sym(:,:,:,5)=flip(flip(flip(job_PRM.bs_void,1),2),3);                      % (5) -x,-y,-z
    bs_sym(:,:,:,6)=flip(job_PRM.bs_void,3);                                      % (6) x,y,-z
    bs_sym(:,:,:,7)=flip(job_PRM.bs_void,2);                                      % (7) x,-y,z
    bs_sym(:,:,:,8)=flip(job_PRM.bs_void,1);                                      % (8) -x,y,z
    bs_sym=mean(bs_sym,4);
    
elseif strcmp(command.sym_op,'4/mmm')==1
    
    bs_sym=zeros(length(job_PRM.xvec),length(job_PRM.yvec),length(job_PRM.zvec),16);
    bs_sym(:,:,:,1)=job_PRM.bs_void;                                              % (1) x,y,z
    bs_sym(:,:,:,2)=flip(flip(job_PRM.bs_void,1),2);                              % (2) -x,-y,z
    bs_sym(:,:,:,3)=permute(flip(job_PRM.bs_void,2),[2 1 3]);                     % (3) -y,x,z
    bs_sym(:,:,:,4)=permute(flip(job_PRM.bs_void,1),[2 1 3]);                     % (4) y,-x,z
    bs_sym(:,:,:,5)=flip(flip(job_PRM.bs_void,1),3);                              % (5) -x,y,-z
    bs_sym(:,:,:,6)=flip(flip(job_PRM.bs_void,2),3);                              % (6) x,-y,-z
    bs_sym(:,:,:,7)=permute(flip(job_PRM.bs_void,3),[2 1 3]);                     % (7) y,x,-z
    bs_sym(:,:,:,8)=permute(flip(flip(flip(job_PRM.bs_void,1),2),3),[2 1 3]);     % (8) -y,-x,-z
    bs_sym(:,:,:,9)=flip(flip(flip(job_PRM.bs_void,1),2),3);                      % (9) -x,-y,-z
    bs_sym(:,:,:,10)=flip(job_PRM.bs_void,3);                                     % (10) x,y,-z
    bs_sym(:,:,:,11)=permute(flip(flip(job_PRM.bs_void,1),3),[2 1 3]);            % (11) y,-x,-z
    bs_sym(:,:,:,12)=permute(flip(flip(job_PRM.bs_void,2),3),[2 1 3]);            % (12) -y,x,-z
    bs_sym(:,:,:,13)=flip(job_PRM.bs_void,2);                                     % (13) x,-y,z
    bs_sym(:,:,:,14)=flip(job_PRM.bs_void,1);                                     % (14) -x,y,z
    bs_sym(:,:,:,15)=permute(flip(flip(job_PRM.bs_void,1),2),[2 1 3]);            % (15) -y,-x,z
    bs_sym(:,:,:,16)=permute(job_PRM.bs_void,[2 1 3]);                            % (16) y,x,z
    bs_sym=mean(bs_sym,4);
    
elseif strcmp(command.sym_op{1},'m-3m')==1
    
    bs_sym=zeros(length(job_PRM.xvec),length(job_PRM.yvec),length(job_PRM.zvec),48);
    bs_sym(:,:,:,1)=job_PRM.bs_void;                                              % (1) x,y,z
    bs_sym(:,:,:,2)=flip(flip(job_PRM.bs_void,1),2);                              % (2) -x,-y,z
    bs_sym(:,:,:,3)=flip(flip(job_PRM.bs_void,1),3);                              % (3) -x,y,-z
    bs_sym(:,:,:,4)=flip(flip(job_PRM.bs_void,2),3);                              % (4) x,-y,-z
    bs_sym(:,:,:,5)=permute(job_PRM.bs_void,[3 1 2]);                             % (5) z,x,y
    bs_sym(:,:,:,6)=permute(flip(flip(job_PRM.bs_void,1),2),[3 1 2]);             % (6) z,-x,-y
    bs_sym(:,:,:,7)=permute(flip(flip(job_PRM.bs_void,1),3),[3 1 2]);             % (7) -z,-x,y
    bs_sym(:,:,:,8)=permute(flip(flip(job_PRM.bs_void,2),3),[3 1 2]);             % (8) -z,x,-y
    bs_sym(:,:,:,9)=permute(job_PRM.bs_void,[2 3 1]);                             % (9) y,z,x
    bs_sym(:,:,:,10)=permute(flip(flip(job_PRM.bs_void,1),2),[2 3 1]);            % (10) -y,z,-x
    bs_sym(:,:,:,11)=permute(flip(flip(job_PRM.bs_void,1),3),[2 3 1]);            % (11) y,-z,-x
    bs_sym(:,:,:,12)=permute(flip(flip(job_PRM.bs_void,2),3),[2 3 1]);            % (12) -y,-z,x
    bs_sym(:,:,:,13)=permute(flip(job_PRM.bs_void,3),[2 1 3]);                    % (13) y,x,-z
    bs_sym(:,:,:,14)=permute(flip(flip(flip(job_PRM.bs_void,1),2),3),[2 1 3]);    % (14) -y,-x,-z
    bs_sym(:,:,:,15)=permute(flip(job_PRM.bs_void,1),[2 1 3]);                    % (15) y,-x,z
    bs_sym(:,:,:,16)=permute(flip(job_PRM.bs_void,2),[2 1 3]);                    % (16) -y,x,z
    bs_sym(:,:,:,17)=permute(flip(job_PRM.bs_void,2),[1 3 2]);                    % (17) x,z,-y
    bs_sym(:,:,:,18)=permute(flip(job_PRM.bs_void,1),[1 3 2]);                    % (18) -x,z,y
    bs_sym(:,:,:,19)=permute(flip(flip(flip(job_PRM.bs_void,1),2),3),[1 3 2]);    % (19) -x,-z,-y
    bs_sym(:,:,:,20)=permute(flip(job_PRM.bs_void,3),[1 3 2]);                    % (20) x,-z,y
    bs_sym(:,:,:,21)=permute(flip(job_PRM.bs_void,1),[3 2 1]);                    % (21) z,y,-x
    bs_sym(:,:,:,22)=permute(flip(job_PRM.bs_void,2),[3 2 1]);                    % (22) z,-y,x
    bs_sym(:,:,:,23)=permute(flip(job_PRM.bs_void,3),[3 2 1]);                    % (23) -z,y,x
    bs_sym(:,:,:,24)=permute(flip(flip(flip(job_PRM.bs_void,1),2),3),[3 2 1]);    % (24) -z,-y,-x
    bs_sym(:,:,:,25)=flip(flip(flip(job_PRM.bs_void,1),2),3);                     % (25) -x,-y,-z
    bs_sym(:,:,:,26)=flip(job_PRM.bs_void,3);                                     % (26) x,y,-z
    bs_sym(:,:,:,27)=flip(job_PRM.bs_void,2);                                     % (27) x,-y,z
    bs_sym(:,:,:,28)=flip(job_PRM.bs_void,1);                                     % (28) -x,y,z
    bs_sym(:,:,:,29)=permute(flip(flip(flip(job_PRM.bs_void,1),2),3),[3 1 2]);    % (29) -z,-x,-y
    bs_sym(:,:,:,30)=permute(flip(job_PRM.bs_void,3),[3 1 2]);                    % (30) -z,x,y
    bs_sym(:,:,:,31)=permute(flip(job_PRM.bs_void,2),[3 1 2]);                    % (31) z,x,-y
    bs_sym(:,:,:,32)=permute(flip(job_PRM.bs_void,1),[3 1 2]);                    % (32) z,-x,y
    bs_sym(:,:,:,33)=permute(flip(flip(flip(job_PRM.bs_void,1),2),3),[2 3 1]);    % (33) -y,-z,-x
    bs_sym(:,:,:,34)=permute(flip(job_PRM.bs_void,3),[2 3 1]);                    % (34) y,-z,x
    bs_sym(:,:,:,35)=permute(flip(job_PRM.bs_void,2),[2 3 1]);                    % (35) -y,z,x
    bs_sym(:,:,:,36)=permute(flip(job_PRM.bs_void,1),[2 3 1]);                    % (36) y,z,-x
    bs_sym(:,:,:,37)=permute(flip(flip(job_PRM.bs_void,1),2),[2 1 3]);            % (37) -y,-x,z
    bs_sym(:,:,:,38)=permute(job_PRM.bs_void,[2 1 3]);                            % (38) y,x,z
    bs_sym(:,:,:,39)=permute(flip(flip(job_PRM.bs_void,2),3),[2 1 3]);            % (39) -y,x,-z
    bs_sym(:,:,:,40)=permute(flip(flip(job_PRM.bs_void,1),3),[2 1 3]);            % (40) y,-x,-z
    bs_sym(:,:,:,41)=permute(flip(flip(job_PRM.bs_void,1),3),[1 3 2]);            % (41) -x,-z,y
    bs_sym(:,:,:,42)=permute(flip(flip(job_PRM.bs_void,2),3),[1 3 2]);            % (42) x,-z,-y
    bs_sym(:,:,:,43)=permute(job_PRM.bs_void,[1 3 2]);                            % (43) x,z,y
    bs_sym(:,:,:,44)=permute(flip(flip(job_PRM.bs_void,1),2),[1 3 2]);            % (44) -x,z,-y
    bs_sym(:,:,:,45)=permute(flip(flip(job_PRM.bs_void,2),3),[3 2 1]);            % (45) -z,-y,x
    bs_sym(:,:,:,46)=permute(flip(flip(job_PRM.bs_void,1),3),[3 2 1]);            % (46) -z,y,-x
    bs_sym(:,:,:,47)=permute(flip(flip(job_PRM.bs_void,1),2),[3 2 1]);            % (47) z,-y,-x
    bs_sym(:,:,:,48)=permute(job_PRM.bs_void,[3 2 1]);                            % (48) z,y,x
    bs_sym=mean(bs_sym,4);
    
else
end

% Rename
job_PRM.bs_sym=bs_sym;

end

