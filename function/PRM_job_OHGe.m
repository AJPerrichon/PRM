function [job_OHG] = PRM_job_OHGe(command,structure)

% Updated OHG function. Now atoms of interest AND reference are to specify
% in command.at, and reference atoms from which to compute the
% distances/etc are ALSO to specify in command.OHG_at. If none is specified
% in command.OHG_at, then command.at atoms will be used as reference (ie
% self term) BEWARE: Cell expansion and reduction steps are extremely slow.
% These steps are performed to find closest OHG_at neighbor from at atoms
% (which might be in a different cell from periodic boundaries).


%% Correct trajectories from periodic boundaries

% Identify the S(ample) and the R(eference)
n_T=length(command.elt_full);
if isempty(command.OHG_at)
    n_S=cell2mat(command.at)-8;
    n_R=n_S;
%     n_A=n_S;
else
    n_S=cell2mat(command.at)-8;
    n_R=cell2mat(command.OHG_at)-8;
%     n_A=union(n_S,n_R);
end

% Cell to supercell expansion
structure.pos_ohg=structure.pos_full(:,:,:,:)-floor(structure.pos_full(:,:,:,:)); % fold trajectory
% structure.pos_ohg=structure.pos_full(:,:,n_A,:)-floor(structure.pos_full(:,:,n_A,:)); % fold trajectory
nrep=1; ncell=((nrep*2)+1)^3; % BEWARE: nrep scales hard
% pos_sc=zeros(structure.length,3,length(n_A),command.nrun,ncell);
pos_sc=zeros(structure.length,3,n_T,command.nrun,ncell);
t=1;
for h=[0,-nrep:1:-1,1:1:nrep] % 0 should be first as (t=1 == initial structure) is required later
    for k=[0,-nrep:1:-1,1:1:nrep]
        for l=[0,-nrep:1:-1,1:1:nrep]
            for m=1:n_T
                for n=1:command.nrun
                    pos_sc(:,:,m,n,t)=cat(2,structure.pos_ohg(:,1,m,n)+h,structure.pos_ohg(:,2,m,n)+k,structure.pos_ohg(:,3,m,n)+l);
                end
            end
            t=t+1;
        end
    end
end
clear t h k l m n;

% Supercell to cell reduction 
% BEWARE: Slow AF % Note: This is not elegant, but it works
pos_red=zeros(structure.length,3,command.nrun,length(n_S),length(n_R)+1);
for n=1:command.nrun
    for h=1:length(n_S) % for all selected atoms
        pos_red(:,:,n,h,1)=pos_sc(:,:,n_S(h),n,1); % store its self-structure in 5-dim index 1
        for k=1:length(n_R) % then for each reference atoms
            for x=1:structure.length  % for each step, search through all cells the shortest distance between R(k) and S(h)
                [~,tt]=min(sqrt(sum((permute(pos_sc(x,:,n_R(k),n,:)-repmat(pos_sc(x,:,n_S(h),n,1),...
                    [1 1 1 1 ncell]),[5 2 1 3 4]).*repmat(structure.cell_parameters,[ncell 1])).^2,2)));
                pos_red(x,:,n,h,k+1)=pos_sc(x,:,n_R(k),n,tt);  % and extract it
            end
        end
    end
end
clear n h k x tt pos_sc ;

%%
% Calculate distances and angle
doh_tmp=zeros(structure.length,command.nrun,length(n_S),length(n_R));
doh_abc=zeros(structure.length-1,command.nrun,length(n_S),8);
    % 1 % Shortest S-R distance (for S=H and R=O, it's covalent bond length)
    % 2 % Second shortest S-R distance (for S=H and R=O, it's HB length 1)
    % 3 % Third shortest S-R distance (for S=H and R=O, it's HB length 2)
    % 4 % Angle between the 1st and 2nd nearest neighbor (for S=H and R=O, it's HB angle 1)
    % 5 % Angle between the 1st and 3rd nearest neighbor (for S=H and R=O, it's HB angle 2)
    % 6 % Sorted 4/5, smallest angle
    % 7 % Sorted 4/5, biggest angle
    % 8 % Distance between the two R closest to S

for n=1:command.nrun
    for h=1:length(n_S)

        % Calculate all distances (will be used to sort n_R atoms)
        for k=2:length(n_R)+1
            for x=1:structure.length
                doh_tmp(x,n,h,k-1)=norm(cat(2,...
                    (pos_red(x,1,n,h,k)-pos_red(x,1,n,h,1)).*structure.cell_parameters(1),...
                    (pos_red(x,2,n,h,k)-pos_red(x,2,n,h,1)).*structure.cell_parameters(2),...
                    (pos_red(x,3,n,h,k)-pos_red(x,3,n,h,1)).*structure.cell_parameters(3)));
            end
        end
        clear k x;

        % At each time step, sort closest neighbors and compute angles
        for x=1:structure.length-1
            
            [M,I]=sort(permute(doh_tmp(x+1,n,h,:),[4 3 2 1]),1);
            vOa=(pos_red(x+1,:,n,h,I(1)+1)-pos_red(x+1,:,n,h,1))';
            vOb=(pos_red(x+1,:,n,h,I(2)+1)-pos_red(x+1,:,n,h,1))';
            vOc=(pos_red(x+1,:,n,h,I(3)+1)-pos_red(x+1,:,n,h,1))';
            doh_abc(x,n,h,1)=M(1); % d(Oa-H)
            doh_abc(x,n,h,2)=M(2); % d(Ob-H)
            doh_abc(x,n,h,3)=M(3); % d(Oc-H)
            doh_abc(x,n,h,4)=rad2deg(atan2(norm(cross(vOa,vOb)),dot(vOa,vOb))); % Oa-H-Ob angle
            doh_abc(x,n,h,5)=rad2deg(atan2(norm(cross(vOa,vOc)),dot(vOa,vOc))); % Oa-H-Oc angle
            
            if x==1 % fix the angle !!!!!!!!!!!!!!!!!!!!!
                if doh_abc(1,n,h,4)<90; doh_abc(1,n,h,4)=180-doh_abc(1,n,h,4); end
                if doh_abc(1,n,h,5)<90; doh_abc(1,n,h,5)=180-doh_abc(1,n,h,5); end
            else%if x>1
%                 if (doh_abc(x,n,h,4)-doh_abc(x-1,n,h,4))>30; doh_abc(x,n,h,4)=doh_abc(x,n,h,4)-45; else end;
%                 if (doh_abc(x,n,h,4)-doh_abc(x-1,n,h,4))<30; doh_abc(x,n,h,4)=doh_abc(x,n,h,4)+45; else end;
%                 if (doh_abc(x,n,h,5)-doh_abc(x-1,n,h,5))>30; doh_abc(x,n,h,5)=doh_abc(x,n,h,5)-45; else end;
%                 if (doh_abc(x,n,h,5)-doh_abc(x-1,n,h,5))>30; doh_abc(x,n,h,5)=doh_abc(x,n,h,5)+45; else end;
            end
            
            if doh_abc(x,n,h,4)<doh_abc(x,n,h,5) % sort by angle instead of distance
                doh_abc(x,n,h,6)=doh_abc(x,n,h,4); % smallest Oa-H_(Ob/Oc)
                doh_abc(x,n,h,7)=doh_abc(x,n,h,5);
            else
                doh_abc(x,n,h,6)=doh_abc(x,n,h,5);
                doh_abc(x,n,h,7)=doh_abc(x,n,h,4);
            end
            
            doh_abc(x,n,h,8)=norm((vOb-vOa).*structure.cell_parameters');
            
        end; clear x M MM I vOa vOb vOc;
    end; clear h;
end; clear n;


% Time-average with or without upper bound
if isempty(command.OHG_utl)
    doh_av=cat(2,permute(mean(mean(doh_abc,2),1),[4 3 2 1]),permute(std(mean(doh_abc,2),1),[4 3 2 1]));
    doh_avr=cat(2,permute(mean(doh_abc,1),[4 3 2 1]),permute(std(doh_abc,1),[4 3 2 1]));
else
    command.OHG_utl=command.OHG_utl/command.time_step;
    n=length(1:command.OHG_utl:structure.length-1)-1;
    doh_tmp=zeros(command.OHG_utl,command.nrun,n_H,n,7);
    for k=1:command.nrun
        for h=1:n_H
            for x=1:n
                doh_tmp(:,k,h,x,:)=doh_abc((x-1)*command.OHG_utl+1:x*command.OHG_utl,k,h,:);
            end
        end
    end; clear l h n k x;
    doh_av=cat(2,permute(mean(mean(doh_tmp,2),1),[5 3 2 4 1]),...
        permute(std(mean(doh_tmp,2),1),[5 3 2 4 1]));
end

%% Histogram of distances
% H--O distances // doh_abc(:,:,:,2)
doh_x=0:0.005:6;
doh_void_1st=zeros(length(doh_x),1);
doh_void_2nd=zeros(length(doh_x),1);
for n=1:command.nrun
    for h=1:length(n_S)
        for x=1:structure.length-1
            [~,doh_L1]=min(abs(doh_x-doh_abc(x,n,h,1)));
            [~,doh_L2]=min(abs(doh_x-doh_abc(x,n,h,2)));
            doh_void_1st(doh_L1,1)=doh_void_1st(doh_L1,1)+1;
            doh_void_2nd(doh_L2,1)=doh_void_2nd(doh_L2,1)+1;
        end
    end
end
doh_void_1st=doh_void_1st./trapz(doh_void_1st);
doh_void_2nd=doh_void_2nd./trapz(doh_void_2nd);
g=@(x,t) x(1).*exp(-log(2).*((t-x(2))./(x(3))).^2); % fitting a gaussian
[doh_void_1st_P,~,~,~,~]=lsqcurvefit(g,[0.02 1.00 0.5],doh_x',doh_void_1st);
[doh_void_2nd_P,~,~,~,~]=lsqcurvefit(g,[0.02 1.95 0.5],doh_x',doh_void_2nd);
clear n h x doh_L g;

%% FFT of distances fluctuations to get characteristic frequencies
% same code as used in the PSD function, see PSD header

% Input parameters
[len,~,~,ncol]=size(doh_abc);
thz_to_mev=4.135665538536;
mev_to_cm=8.06554;
command.cut=command.cut/command.time_step;
command.jres=3;
slep=dpss(command.cut,command.jres);
dslep=ceil(0.7*command.cut/sqrt(command.jres+0.3));
vslep=1:dslep:len; 
nslep=sum((vslep+command.cut-1)<=len);

% Reshape into segments of length command.cut
doh_true=zeros([command.cut,command.nrun,length(n_S),ncol,nslep]);
for x=1:nslep
    doh_true(:,:,:,:,x)=doh_abc(vslep(x):vslep(x)+command.cut-1,:,:,:);
end; clear x;

% FFT
doh_fft=zeros([command.cut,command.nrun,length(n_S),8,nslep]);
for n=1:command.nrun
    for h=1:length(n_S)
        for y=1:ncol
            for x=1:nslep
                yg=(doh_true(:,n,h,y,x)-mean(doh_true(:,n,h,y,x))).*slep(:,1);
                doh_fft(:,n,h,y,x)=fft(yg,command.cut);
            end
        end
    end
end; clear x y h n ncol yg;

% Power sum
doh_ps=permute(mean(mean(doh_fft,5),2),[1 4 3 2]);
doh_ps=abs(doh_ps(1:command.cut/2+1,:,:));
doh_en=((((1/command.time_step*10^3)/2*linspace(0,1,command.cut/2+1)).*thz_to_mev)');
if strcmp(command.unit,'cm') || strcmp(command.unit,'cm-1'); doh_en=doh_en.*mev_to_cm; end


%% Demodulation
if command.OHG_demod==1
    
    % Demodulation of O-H covalent bond distances
    doh_demod=zeros(len,command.nrun,n_H,5);
    for n=1:command.nrun
        for h=1:n_H
            [Mtmp,Itmp]=findpeaks(doh_abc(:,n,h,1));
            doh_demod(:,n,h,1)=spline(Itmp,Mtmp,(1:1:len)'); % upper envelope
            [Mtmp,Itmp]=findpeaks(doh_abc(:,n,h,1).*(-1));
            doh_demod(:,n,h,2)=spline(Itmp,Mtmp.*(-1),(1:1:len)'); % lower envelope
            doh_demod(:,n,h,3)=(doh_demod(:,n,h,1)-doh_demod(:,n,h,2))./2; % envelope amplitude
            doh_demod(:,n,h,4)=(doh_demod(:,n,h,1)+doh_demod(:,n,h,2))./2; % envelope center
            doh_demod(:,n,h,5)=(doh_abc(:,n,h,1)-doh_demod(:,n,h,4))./doh_demod(:,n,h,3); % pure carrier
        end
    end; clear n h x Mtmp Itmp;
    
	% Zero-crossing of carrier and associated time-lag
    doh_zc_tmp=cell(1,command.nrun,n_H);
    for n=1:command.nrun
        for h=1:n_H
            [~,Izc]=findpeaks(abs(doh_demod(:,n,h,5)).*(-1)); % abs to fold from -1:1 to 0:1; .*(-1) to turn minima as peaks;
            % Dzc average dO-H over delta Izc
            Dzc=zeros(length(Izc)-1,1);
            for k=2:length(Izc) % Izc is index of zero crossing
                Dzc(k-1,1)=mean(doh_demod(Izc(k-1,1):Izc(k,1),n,h,4),1);
            end
            % Tzc half period; time between each zero crossing
            Tzc=(Izc(2:end,1)-Izc(1:length(Izc)-1)).*command.time_step;
            doh_zc_tmp{1,n,h}(:,1)=Dzc;
            doh_zc_tmp{1,n,h}(:,2)=Tzc;
        end
    end; clear n h k Izc Dzc Tzc;
    T=(command.time_step:command.time_step:100)'; % grid of time
    D=(min(min(min(doh_abc(:,:,:,1)))):0.001:max(max(max(doh_abc(:,:,:,1)))))'; % grid of dO-H
    doh_zc=zeros(command.nrun,n_H,length(T),length(D));
    for n=1:command.nrun
        for h=1:n_H
            for k=1:length(doh_zc_tmp{1,n,h}(:,1))
                [~,Dx]=min(abs(doh_zc_tmp{1,n,h}(k,1)-D));
                [~,Tx]=min(abs(doh_zc_tmp{1,n,h}(k,2)-T));
                doh_zc(n,h,Tx,Dx)=doh_zc(n,h,Tx,Dx)+1; % histogram
            end
        end
    end; clear n h k Dx Tx doh_zc_tmp;
    doh_zc=permute(sum(sum(doh_zc,1),2),[3 4 2 1]);
    
    % PSD of demodulated signal
    [~,~,~,ncol]=size(doh_demod);
    doh_demod_tmp=zeros([command.cut,command.nrun,n_H,ncol,nslep]);
    for x=1:nslep
        doh_demod_tmp(:,:,:,:,x)=doh_demod(vslep(x):vslep(x)+command.cut-1,:,:,:);
    end; clear x;
    doh_demod_fft=zeros([command.cut,command.nrun,n_H,3,nslep]);
    X=(1:1:command.cut)';
    for n=1:command.nrun
        for h=1:n_H
            for x=1:nslep
                [a1,~]=polyfit(X,doh_demod_tmp(:,n,h,3,x),1);
                yg=(doh_demod_tmp(:,n,h,3,x)-(X.*a1(1)+a1(2))).*slep(:,1); % amplitude
                [a2,~]=polyfit(X,doh_demod_tmp(:,n,h,4,x),1);
                yh=(doh_demod_tmp(:,n,h,4,x)-(X.*a2(1)+a2(2))).*slep(:,1); % center
                [a3,~]=polyfit(X,doh_demod_tmp(:,n,h,5,x),1);
                yk=(doh_demod_tmp(:,n,h,5,x)-(X.*a3(1)+a3(2))).*slep(:,1); % carrier
                doh_demod_fft(:,n,h,1,x)=fft(yg,command.cut);
                doh_demod_fft(:,n,h,2,x)=fft(yh,command.cut);
                doh_demod_fft(:,n,h,3,x)=fft(yk,command.cut);
            end
        end
    end; clear n h x a1 a2 a3 X yg yh yk;
    doh_demod_ps=permute(mean(mean(doh_demod_fft,5),2),[1 4 3 2]);
    doh_demod_ps=abs(doh_demod_ps(1:command.cut/2+1,:,:));
    
else
    doh_demod=[]; doh_demod_ps=[]; doh_zc=[];
end

% Rename
job_OHG.doh_abc=doh_abc;
job_OHG.doh_av=doh_av;
job_OHG.doh_avr=doh_avr;
job_OHG.doh_en=doh_en;
job_OHG.doh_ps=doh_ps;
%
job_OHG.doh_x=doh_x;
job_OHG.doh_void_1st=doh_void_1st;
job_OHG.doh_void_2nd=doh_void_2nd;
job_OHG.doh_void_1st_P=doh_void_1st_P;
job_OHG.doh_void_2nd_P=doh_void_2nd_P;
%
job_OHG.doh_demod=doh_demod;
job_OHG.doh_demod_ps=doh_demod_ps;
job_OHG.doh_zc=doh_zc;

end
