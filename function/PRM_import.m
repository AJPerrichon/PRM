function [command,structure] = PRM_import(command)

% Trajectories should be stored in a \dat subfolder. Import function to get
% lattice parameters and trajectories from Vasp' XDATCAR file, and set up
% properly the "t=0" model depending on the command.space tag. Works for
% XDATCAR 'old' format (VASP 5.2), and 'new' format (VASP 5.3+); for the
% new format, string lines 'Direct configuration= n' are skipped using
% 'CommentStyle','Direct' in textscan function. If several XDATCAR are
% provided for a single command.name, in the form [command.name,'_x'] with
% x an integer - in that case command.name 'by itself' should not exist
% (check lines 18-27) - then all these equivalent yet strictly
% non-concatenable trajectories will be processed independently. This means
% that a fourth dimension is added to the 'structure.pos' variable (note
% that all generation MD must have the same length, as they are stored into
% tensors and not cells for time-efficiency reasons - they will be cutoff
% to the shortest multi-trajectory). For test purposes, a 'merged' option
% to command.importmode is available, to concatenate multiple trajectories;
% note that it's okay-ish when the atomic positions are used, but
% extremelly wrong for function that use the derivative (velocities) or
% time-lag (correlations).

if strcmp(command.importmode,'default') || strcmp(command.importmode,'merged')

    % Check if single or multiple XDATCAR
    command.pathway=['dat\',command.name];
    if strcmp(command.nrun,'auto')==1 || isempty(command.nrun)==1
        a=dir('dat'); b=struct2cell(a);
        if any(ismember(b(1,:),command.name))
            command.nrun=1;
        else
            c=zeros(1,100);
            for k=1:100
                c(k)=any(ismember(b(1,:),[command.name,sprintf('_%.0f',k)]));
            end
            command.nrun=sum(c);
        end; clear a b c k;
    else
    end

    % Build multiple XDATCAR paths
    if command.nrun>1
        for k=1:command.nrun; command.truepath{1,k}=sprintf('%s_%1.0f',command.pathway,k); end; clear k;
    else
    end

    % Read number of atoms & Number of inequivalent atoms
    if command.nrun>1; fid=fopen(strrep(command.truepath{1,1},'\',filesep)); else fid=fopen(strrep(command.pathway,'\',filesep)); end
    linenum=7; line_tmp=textscan(fid,'%s',1,'delimiter','/n','headerlines',linenum-1);
    n_in=sscanf(line_tmp{1,1}{1,1},'%f'); n_to=sum(n_in);
    fclose(fid); clear fid linenum line_tmp;

    % Read atom type
    if command.nrun>1; fid=fopen(strrep(command.truepath{1,1},'\',filesep)); else fid=fopen(strrep(command.pathway,'\',filesep)); end
    label_tmp=textscan(fid,'%s%s%s%s%s%[^\n\r]',1,'Delimiter',' ','MultipleDelimsAsOne',true,'HeaderLines',5,'ReturnOnError',false);
    command.label=[label_tmp{1:end-1}];
    fclose(fid); clear fid ans label_tmp;

    % Identify atoms to command.label
    plabel=zeros(max(n_in),length(n_in));
    for k=1:length(n_in); plabel(1:n_in(k),k)=k; end; clear k;
    plabel=plabel(:); plabel(plabel==0)=[];

    % Read cell parameters
    if command.nrun>1; fid=fopen(strrep(command.truepath{1,1},'\',filesep)); else fid=fopen(strrep(command.pathway,'\',filesep)); end
    linenum=3; line_tmp=textscan(fid,'%s',3,'delimiter','/n','headerlines',linenum-1);
    lat=[sscanf(line_tmp{1,1}{1,1},'%f'),sscanf(line_tmp{1,1}{2,1},'%f'),sscanf(line_tmp{1,1}{3,1},'%f')];
    fclose(fid); clear fid linenum line_tmp;

    % Import trajectories
    if command.nrun>1
        for k=1:command.nrun
            fid=fopen(strrep(command.truepath{1,k},'\',filesep));
            S_tmp=textscan(fid,'%f %f %f','headerlines',7,'CommentStyle','Direct');
            S{1,k}=[S_tmp{1,1}(:,1),S_tmp{1,2}(:,1),S_tmp{1,3}(:,1)];
            fclose(fid); clear fid S_tmp ans;
            len_tmp(k)=length(S{1,k});
        end; clear k;
        % cutoff to homogeneous len for cell2mat conversion
        len=min(len_tmp);
        for k=1:command.nrun
            S{1,k}=S{1,k}(1:len,:);
        end; clear k len_tmp;
        S=reshape(cell2mat(S),[len,3,command.nrun]);
    else
        fid=fopen(strrep(command.pathway,'\',filesep));
        S_tmp=textscan(fid,'%f %f %f','headerlines',7,'CommentStyle','Direct');
        S=[S_tmp{1,1}(:,1),S_tmp{1,2}(:,1),S_tmp{1,3}(:,1)];
        fclose(fid); clear fid S_tmp ans;
        len=length(S);
    end

    % Extract data for selected atoms
    if command.nrun>1
        for k=1:command.nrun
            if strcmp(command.space,'com') || strcmp(command.space,'real') || strcmp(command.space,'ub')
                SGS{1,k}=S(1:n_to,:,k);
            elseif strcmp(command.space,'poscar')
                command.pathtmp=sprintf('%s%s%s_%1.0f',command.pathway(1:4),'POSCAR',command.pathway(12:end),1);
                fid=fopen(strrep(command.pathtmp,'\',filesep));
                SGS_tmp=textscan(fid,'%f %f %f','headerlines',8);
                SGS{1,k}=[SGS_tmp{1,1}(:,1),SGS_tmp{1,2}(:,1),SGS_tmp{1,3}(:,1)];
                fclose(fid); clear fid SGS_tmp;
            else
            end
        end
        if strcmp(command.space,'poscar'); command.space='com'; else end
        SGS=reshape(cell2mat(SGS),[n_to,3,command.nrun]);
    else
        if strcmp(command.space,'com') || strcmp(command.space,'real') || strcmp(command.space,'ub')
            SGS=S(1:n_to,:);
        elseif strcmp(command.space,'poscar')
            command.pathtmp=sprintf('%s%s%s',command.pathway(1:4),'POSCAR',command.pathway(12:end));
            fid=fopen(strrep(command.pathtmp,'\',filesep));
            SGS_tmp=textscan(fid,'%f %f %f','headerlines',8);
            SGS=[SGS_tmp{1,1}(:,1),SGS_tmp{1,2}(:,1),SGS_tmp{1,3}(:,1)];
            fclose(fid); clear fid SGS_tmp; command.space='com';
        else
        end
    end

    if isempty(command.time_stop); command.time_stop=len/n_to; end
    if command.time_stop>(len/n_to); command.time_stop=len/n_to; end

    if command.nrun>1
        S=[SGS;S(n_to*(command.time_start+1)+1:n_to*(command.time_stop),:,:)];
        len=length(S);
        S=reshape(permute(S,[2 1 3]),[3,n_to,len/n_to,command.nrun]); S=permute(S,[3 1 2 4]);
        pos=S(:,:,cell2mat(command.at)-8,:);
    else
        S=[SGS;S(n_to*(command.time_start+1)+1:n_to*(command.time_stop),:)];
        len=length(S);
        S=reshape(permute(S,[2 1]),[3,n_to,len/n_to]); S=permute(S,[3 1 2]);
        pos=S(:,:,cell2mat(command.at)-8);
    end

    % Get atom labels
    for k=1:n_to; command.elt_full{1,k}=command.label(plabel(k,1)); end % full structure
    plabel=plabel(cell2mat(command.at)-8);
    for k=1:length(command.at); command.elt{1,k}=command.label(plabel(k,1)); end % selected atoms
    clear SGS len n_to;
    
    % Concatenate all trajectories into one
    if command.nrun>1 && strcmp(command.importmode,'merged')
        ren=pos(:,:,:,1);
        for k=2:command.nrun
            ren=cat(1,ren,pos(:,:,:,k));
        end; clear k;
        pos=ren;
        command.nrun=1;
    else
    end
    
    % Get trajectory length
    if command.nrun>1
        [len,~,~,~]=size(pos);
    else
        [len,~,~]=size(pos);
    end
    
    % Rename
    structure.pos=pos;
    structure.pos_full=S;
    structure.lattice=lat;
    structure.length=len;
    
else
end

end

