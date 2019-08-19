% Post-processing of molecular dynamics trajectories
    % Job PRM,	Positional recurrence maps (2D histogram)
    % Job hX,	Histogram of mean square displacements and atomic positions
    % Job PSD,	Power spectral density
    % Job MSD,  Mean square displacements & Diffusion constants
    % Job Kin,  Effective temperature from kinetic energy
    % Job OHG,	O-H geometry (selected atoms: single or multi H, all O, optional 3rd atom type)
    % Job PSDf, Power spectral density with filtered velocities
    % command.space for job PRM and hX
        % 'real' for real cartesian space (positions), where PRM got static & dynamical components 
        % 'com' for center-of-mass space (displacements), where PRM got dynamical component only
        % 'ub' for unbent space, a 'standard' nuclear density maps without symmetry
        % 'poscar' use structure from POSCAR instead of the first step of the MD (*)
        % (*) usefull in 'com' when MD doesn't start from the optimized structure, e.g. 
        %     (Opt. geom. > NVT thermalization > NVE generation)
        % (*) can also be used as 'real' (or mixture com/real) by putting the 'model' in POSCAR format...
    % see function headers for further descriptions /function/*
    % see readme.txt for version


%%
%%%%%%%%%%%%%%%%%%%%
% Command parameters
%%%%%%%%%%%%%%%%%%%%

% Initialisation
close all;
clear;
clc;
addpath function;

% MD trajectory data
command.importmode='default';        % see /function/PRM_import* % default/merged
command.name='XDATCAR_6b';      % name of XDATCAR file in /dat/
command.nrun=[];                     % [] % only use the first nrun=integer trajectories
command.time_step=0.5;               % time step of the molecular dynamics in fs
command.time_start=500;              % number of steps to remove at the beginning of simulation
command.time_stop=[];                % [] % step number to stop at

% Structural data
command.at=num2cell(9:36);   % selected atoms (XDATCAR lines)
command.supercell=[2,2,2];           % supercell dimension for figure scaling
command.angle=0;                     % z-rotation of the data % see /function/PRM_rotate
command.space='com';                  % real/com/poscar/ub % see header
command.space_grid=[1/2,1/4,1/4;1/4,1/4,0;-1/4,1/4,0;1/4,0,1/4]; % see /function/PRM_real

% Job parameters (overall)
command.PRM_resolution=0.03;         % global resolution % default: 0.03
command.sym=1;                       % 0/1 % without/with symmetry average
command.sym_op={'m-3m'};              % see /function/PRM_sym
command.cut=2048;                    % reduced trajectory length in fs % 1024/2048 (2^n)
command.conv_str=[0.1125,0.025];          % gaussian width for convolution % [flat,scaling] % IN1 = [0.1125,0.025];
command.bweight='i';                 % cross section; 'i' incoherent, 'c' coherent, 't' total, 'n' none
command.unit='mev';                  % Energy unit; mev/cm
% Job PRM
command.PRM_direction=1;             % 1 for x-proj, 2 for y-proj, 3 for z-proj
% Job PSD
command.PSD_nts=[];                  % default: [] == command.time_step
command.PSD_partial='';              % '' no partial / 'elt' per elements / 'dir' per directions / 'gen' per generation;
command.PSD_instrument='IN1';        % to get Q-dependence for DW estimate
% Job MSD
command.MSD_sub=0;                   % Subtract closest 2nd element (H/O only for now)
command.MSD_filter=0;                % Subtract lattice contribution with pass-band filter
% Job OHG
command.OHG_at=num2cell([11,12,15,16,19,20,23,24]);       % selected reference atoms amongst command.at % [] for all
command.OHG_utl=[];                  % upper time limit for time-average % default: []
command.OHG_3rdM=0;                  % 0/1 % without/with a 3rd atom type
command.OHG_demod=0;                 % demodulation routine
command.OHG_load=0;                  % 0 run OHG / 1 load data_OHG
% Job PSDf
command.PSDf_filter=1;               % 0/1 without/with filter (moving average)
command.PSDf_window='none';          % 'none'/'slepian'/'kaiser'/'tukey' % default: 'none'

% DO Job & Plot Figures
command.plot=zeros(7,2);             % (job_num,options) % see header
command.plot(1,:) = [1,0];           % Job PRM, [job&plot,plot3d]
command.plot(2,1) = 1;               % Job hX, h(MSD) and h(POS)
command.plot(3,1) = 1;               % Job PSD, power spectral density and Sqw estimate
command.plot(4,1) = 0;               % Job MSD, mean square displacement and diffusion constant
command.plot(5,1) = 1;               % Job Kin, effective temperature from kinetic energy
command.plot(6,1) = 0;               % Job OHG, geometry of the hydroxide group
command.plot(7,1) = 0;               % Job PSDf, PSD on filtered trajectories (requires job OHG)

% Printing
command.print=0;                     % 0/1 % print outcome in .txt files
command.print_n='';                  % label of selected atoms % Oap, H1, res5meV, etc...

%%
%%%%%%%%%%%%%
% Import data
%%%%%%%%%%%%%

disp('Import data, started...'); tic;
[command,structure]=PRM_import(command);
% call your own function to get {structure.pos,structure.lattice,command.elt} if data are not in XDATCAR format
disp('... done.'); toc; disp(' ');

%%
%%%%%%%%%%%%%%
% Prepare data 
%%%%%%%%%%%%%%

% Correct data from periodic boundaries 
disp('Correct data from periodic boundaries, started...'); tic;
structure=PRM_periodicity(command,structure); disp('... done.'); toc; disp(' ');

% Convert lattice from matrice to norm+angle
if exist('structure.lat_val','var')==0 && exist('structure.lat_ang','var')==0; [structure]=PRM_lattice(structure); end
if sum(structure.cell_angles~=90)~=0; command.ortho=0; else; command.ortho=1; end

% Convert non-orthogonal crystal to 'virtually' orthogonal crystal
% not implemented yet
if command.ortho==0; end

% Convert from center-of-mass space to cartesian
if strcmp(command.space,'real')
    disp('Convert data from center-of-mass to cartesian space, started...'); tic
    structure=PRM_real(command,structure); disp('... done.'); toc; disp(' ');
else
    structure.cartesian=zeros(1,3,length(command.at));
end

% Rotate data
if command.angle~=0
    disp('Rotate data, started...'); tic;
    [command,structure]=PRM_rotate(command,structure);
    disp('... done.'); toc; disp(' ');
else
end

% Bend space, for job 1:3
if command.plot(1,1)==1 || command.plot(2,1)==1 || command.plot(3,1)==1
	[command,structure]=PRM_bender(command,structure);
else
end

% Normalize total count; Set figure flag
structure.total_count=(structure.length-1)*length(command.at)*command.nrun;
fig_flag=1;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Processing, plot, print jobs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% job PRM, positional recurrence maps
if command.plot(1,1)==1
    disp('Job, Positional Recurrence Maps, started...'); tic
    % Run job
    [command,job_PRM]=PRM_job_PRM(command,structure);
    % Symmetry expansion
    if command.sym==1; job_PRM=PRM_job_PRM_sym(command,job_PRM); end
    % Projection and log scale
    if command.sym==1; job_PRM.bs_plan=sum(job_PRM.bs_sym(:,:,:),command.PRM_direction); else; job_PRM.bs_plan=sum(job_PRM.bs_void(:,:,:),command.PRM_direction); end
    if command.PRM_direction==1; job_PRM.bs_plan=reshape(job_PRM.bs_plan,[job_PRM.yvec_length,job_PRM.zvec_length]); end
    if command.PRM_direction==2; job_PRM.bs_plan=reshape(job_PRM.bs_plan,[job_PRM.xvec_length,job_PRM.zvec_length]); end
    if command.PRM_direction==3; job_PRM.bs_plan=reshape(job_PRM.bs_plan,[job_PRM.xvec_length,job_PRM.yvec_length]); end
    job_PRM.bs_plan=permute(job_PRM.bs_plan,[2 1]); job_PRM.bs_log=log(job_PRM.bs_plan);
    % Plot
    fig_flag=PRM_job_PRM_plot  (command,structure,job_PRM,fig_flag);
    fig_flag=PRM_job_PRM_plot3D(command,structure,job_PRM,fig_flag); 
    disp('... done.'); toc; disp(' ');
end

% job hX, histogram of mean square displacements and atomic positions
if command.plot(2,1)==1
    disp('Job, Histogram of MSD and Atomic Positions, started...'); tic;
    % Run job
    job_hX=PRM_job_hX(command,structure);
    % Plot
    fig_flag=PRM_job_hX_plot(command,job_hX,fig_flag);
    % Print
    PRM_job_hX_print(command,job_hX);
    disp('... done.'); toc; disp(' ');
end

% job PSD, power spectral density
if command.plot(3,1)==1
    disp('Job, Power Spectral Density, started...'); tic;
    % Run job
    if isempty(command.PSD_nts); command.PSD_nts=command.time_step; end
    [command,job_PSD]=PRM_job_PSDe(command,structure);
    % Plot
    fig_flag=PRM_job_PSD_plot(command,job_PSD,fig_flag);
    % Print
    PRM_job_PSD_print(command,job_PSD);
    disp('... done.'); toc; disp(' ');
end

% job MSD, mean square displacements and diffusion constant
if command.plot(4,1)==1
    disp('Job, MSD, started...'); tic;
    % Run job & Plot
    job_MSD=PRM_job_MSD(command,structure);
    disp('... done.'); toc; disp(' ');
end

% job Kin, effective temperature from kinetic energy
if command.plot(5,1)==1
    disp('Job, Effective temperature, started...'); tic;
    % Run job
    job_Kin=PRM_job_Kin(command,structure);
    % Plot
    fig_flag=PRM_job_Kin_plot(job_Kin,fig_flag);
    disp('... done.'); toc; disp(' ');
end

% job OHG, geometry of the O-H--O bonds
if command.plot(6,1)==1
    disp('Job, O-H geometry, started...'); tic;
    % Check
    if command.angle~=0; error('Job OHG: command.angle must be 0 due to SC expansion'); end
    % Run job
    if command.OHG_load==0
        job_OHG=PRM_job_OHGe(command,structure);
        save data_OHG job_OHG;
    else
        load data_OHG;
    end
    % Plot
    fig_flag=PRM_job_OHGe_plot(command,job_OHG,fig_flag);
    disp('... done.'); toc; disp(' ');
end

% job PSDf, PSD on filtered velocities (requires job OHG)
if command.plot(7,1)==1
    disp('Job, PSDf, started...'); tic;
    % Check
    if command.plot(6,1)==1; else; error('Job OHG must be performed for Job PSDf to work'); end
	% Run Job
    job_PSDf=PRM_job_PSDf(command,structure,job_OHG);
	% Plot
	fig_flag=PRM_job_PSDf_plot(command,job_PSDf,fig_flag);
    disp('... done.'); toc; disp(' ');
end
