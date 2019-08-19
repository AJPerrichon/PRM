function [command,job_PSD] = PRM_job_PSD(command,structure)

% Compute the power spectral density (PSD) of selected atoms. Option to
% convolute PS with a gaussian with constant and/or scaling width with
% command.conv_str to reproduce experimental resolution: basically the
% gaussian width is defined by b+ax, with b==command.conv_str(1) in meV
% (0-5 meV) and a==command.conv_str(2) a rate (0-0.1).
% Velocities are weighted with scattering lengths (thus several atom types
% can be selected at the same time), and each partial power spectrum is
% normalized to 1. The choice of weighting the velocities by the scattering
% lengths, instead of the power spectra by the cross section, is motivated
% by the mixing of the modes (distinct atoms contribute to a single mode,
% which cannot be modeled easily).
% In order to avoid spurious correlation in the FFT, both the 'augmented'
% sampling frequency (command.nts, in fs) and the 'reduced' length of data
% (command.cut, in fs), must be powers of 2: 2^k, typically with k=-1:2 for
% command.nts and k=9:12 for command.cut. Note that using an 'augmented'
% time step has absolutely no impact on the variance, and since the
% trajectories are short, mostly no impact on computational time (which is
% proportional to sqrt(N); the gain is only a few milliseconds); however,
% it can be used as a 'trick' to get a time step which is a power of two:
% e.g. 'native' time step is 0.4 fs, merging 5 points leads to a 2^1
% 'augmented' time step. Concering the length of the simulation, one should
% play arround with command.cut to get a spacing in energy short enough to
% have well defined peaks. Using this technique, the variance is reduced by
% a factor K, where K is the number of segments.
% In order to fix the issues of discrete FFT of discrete data, we use a
% slepian window function; to insure no frequency leakage during the FFT
% process, we rather use the strategy of taking only the smallest
% eigenvector (k=0) of the Slepian and overlapping data segments, rather
% than taking the first 2-3 eigenvectors of the Slepian without
% overlapping. We strongly suggest to use command.jres=3 for the halfwidth
% of the Slepian, which give the best compromise between loss in resolution
% of the main peak (3 times the resolution - which should still be smaller
% than any kind of experimental resolution), and no leakage at all. If the
% resolution is too large, command.jres can be decreased to 2 (while the
% leakage won't be strictly 0, the highest sidelobe level will only be
% 10^-5 of the main). Note that there is actually no reason to take
% command.jres>3, if you want to increase the resolution to match the
% experimental value, it is 'cleaner' to do it though the convolution.
% Partial PSD are available through the tag 'command.part'; the only
% difference from the standard case is that the related dimensions of the
% ve_fft 6-dim tensor are not averaged... 'dir' is divided by 3 for an
% easier visualization, while 'elt' is arbitrarily not furtherly weighted:
% it can easily be changed on line 127.
% Please also note that Matlab deletes singleton dimension if dim>2. In
% order to avoid to code 'exceptions' for singleton dimensions, we chose to
% permute the data in an counter-intuitive way to ensure those potential
% singleton dimensions are preserved.

% Input parameters
    % command.at;        selected atoms
    % command.nts;       augmented time step in fs, default command.nts=command.time_step
    % command.cut;       length of data segments in fs
    % command.time_step; native time step in fs
    % command.conv_str;  gaussian width, see header
    % structure.pos;     reshaped MD trajectories
    % command.nrun;      number of independent trajectories
    % command.elt;       selected atoms label
    % command.bweight;   type of b length: incoherent/total/coherent/none
    % command.unit;      energy unit in mev/cm-1
    % command.part;      enable partial PSD, for elements or directions


% Setting things up
thz_to_mev=4.135665538536;
mev_to_cm=8.06554;
q2toE=2.072124365744764;

%% Prepare
command.jres=3; % see header
[len,~,~,~]=size(structure.pos);

command.PSD_nts=command.PSD_nts/command.time_step; % conversion of command.nts from fs to simulation step
if command.PSD_nts~=2^nextpow2(command.PSD_nts); error('command.nts*command.time_step must be a power of 2'); else end
command.cut=command.cut/command.time_step; % conversion of command.cut from fs to simulation step
if command.cut~=2^nextpow2(command.cut); error('command.cut must be a power of 2'); else end

% Get velocity matrix 
    % here len-2 instead of len-1 because first step of pos should be the
    % structural model needed for job 1-3, and thus velocity from
    % step2-step1 cannot be calculated
structure.pos=permute(structure.pos,[4,3,2,1]);
ve_void(:,:,:,1:len-2)=(structure.pos(:,:,:,3:len)-structure.pos(:,:,:,2:len-1))./command.time_step;
[~,~,~,len]=size(ve_void);

% Average trajectory to match augmented time step
ve_true=zeros([command.PSD_nts,command.nrun,length(command.at),3,floor(len/command.PSD_nts)]);
for x=1:command.PSD_nts
    for y=1:floor(len/command.PSD_nts)
        ve_true(x,:,:,:,y)=ve_void(:,:,:,y*command.PSD_nts-x+1);
    end
end; clear y x;
ve_true=mean(ve_true,1);
[~,~,~,~,len]=size(ve_true);

% Window function (slepian) and overlapping data segments
command.cut=command.cut/command.PSD_nts;
slep=dpss(command.cut,command.jres);
dslep=ceil(0.7*command.cut/sqrt(command.jres+0.3));
vslep=1:dslep:len; nslep=sum((vslep+command.cut-1)<=len);

ve_truer=zeros([nslep,1,command.nrun,length(command.at),3,command.cut]);
for x=1:nslep
    ve_truer(x,:,:,:,:,:)=ve_true(:,:,:,:,vslep(x):vslep(x)+command.cut-1);
end; clear x;

%% from R(t) to G(w) one-phonon incoherent
% Fast Fourier Transform (FFT)
ve_fft=zeros([command.cut,nslep,1,command.nrun,length(command.at),3]);
for x=1:nslep
    for y=1:command.nrun
        for h=1:length(command.at)
            elt=command.elt{1,h}; c=PRM_xlength(elt,command);
            for k=1:3
                yg=permute(ve_truer(x,1,y,h,k,:),[6 1 2 3 4 5]).*slep(:,1).*c;
                ve_fft(:,x,1,y,h,k)=fft(yg,command.cut);
            end
        end
    end
end; clear k h c elt y x;

% Averaging & Normalization; checking if partial PSD are used
if strcmp(command.PSD_partial,'elt') % per element type
    ve_psum=permute(sum(sum(sum(sum(abs(ve_fft).^2,2),3),4),6),[1 5 6 4 3 2]);
    ve_ps=ve_psum(1:command.cut/2+1,:);
    ve_ps_iso=mean(ve_ps,2)./sum(sum(ve_ps,2),1);
    command.elt_tmp=cell(length(command.elt),1);
    for k=1:length(command.elt); command.elt_tmp{k,1}=command.elt{1,k}{1,1}; end
    [~,idx]=unique(command.elt_tmp,'first'); command.elt_rn=command.elt_tmp(sort(idx));
    command.elt_mat=zeros(length(command.elt),length(command.elt_rn));
    for h=1:length(command.elt_rn)
        for k=1:length(command.elt)
            command.elt_mat(k,h)=strcmp(command.elt{1,k},command.elt_rn{h,1});
        end
    end
    ve_ps_proj=zeros(length(ve_ps_iso),length(command.elt_rn));
    for k=1:length(command.elt_rn)
        ve_ps_proj(:,k)=sum(ve_ps(:,command.elt_mat(:,k)==1),2);
    end
    ve_ps=cat(2,ve_ps_iso,ve_ps_proj);
    
elseif strcmp(command.PSD_partial,'dir') % per cartesian direction
    ve_psum=permute(sum(sum(sum(sum(abs(ve_fft).^2,2),3),4),5),[1 6 5 4 3 2]);
    ve_ps=ve_psum(1:command.cut/2+1,:);
    ve_ps_iso=mean(ve_ps,2)./sum(sum(ve_ps,2),1);
    ve_ps_proj=ve_ps./sum(sum(ve_ps,2),1)./3;
    ve_ps=cat(2,ve_ps_iso,ve_ps_proj);
    
elseif strcmp(command.PSD_partial,'gen') % per trajectory
    ve_psum=permute(sum(sum(sum(sum(abs(ve_fft).^2,2),3),5),6),[1 4 6 5 3 2]);
    ve_ps=ve_psum(1:command.cut/2+1,:);
    ve_ps_iso=mean(ve_ps,2)./sum(sum(ve_ps,2),1);
    ve_ps_proj=ve_ps./sum(sum(ve_ps,2),1);
    ve_ps=cat(2,ve_ps_iso,ve_ps_proj);
    
else % normal case, gDOS G(w)
    ve_psum=sum(sum(sum(sum(sum(abs(ve_fft).^2,2),3),4),5),6);
    ve_ps=ve_psum(1:command.cut/2+1,1);
    ve_ps=ve_ps./max(ve_ps);
end

% Energy vector
ve_en=(((1/(command.time_step*command.PSD_nts)*10^3)/2*linspace(0,1,command.cut/2+1)).*thz_to_mev)';

%% Overtones
[~,dim]=size(ve_ps);
ve_pse=zeros([size(ve_ps),7]);
for d=1:dim
    ve_pse(:,d,1)=ve_ps(:,d);
    ve_pse(:,d,2)=interp1(ve_en.*2,ve_pse(:,d,1),ve_en)./(2)./2;
    ve_pse(:,d,3)=interp1(ve_en.*3,ve_pse(:,d,1),ve_en)./(3*2)./3;
    ve_pse(:,d,4)=interp1(ve_en.*4,ve_pse(:,d,1),ve_en)./(4*3*2)./4;
    ve_pse(:,d,5)=interp1(ve_en.*5,ve_pse(:,d,1),ve_en)./(5*4*3*2)./5;
    ve_pse(:,d,6)=interp1(ve_en.*6,ve_pse(:,d,1),ve_en)./(6*5*4*3*2)./6;
    ve_pse(:,d,7)=interp1(ve_en.*7,ve_pse(:,d,1),ve_en)./(7*6*5*4*3*2)./7;
end
ve_pse=sum(ve_pse,3);

%%
% Real space discrete convolution
if command.conv_str(1)~=0 || command.conv_str(2)~=0
    
    % Generate extended energy domain
    command.gspace=199; % gaussian space for convolution, need to be odd (default 199)
    t=(((1/(command.time_step*command.PSD_nts)*10^3)/2*linspace(0,2,command.cut+1)).*thz_to_mev)';
    sp=length(ve_en)+floor(command.gspace/2);
    ve_en=t(1:sp,1); dt=ve_en(2,1)-ve_en(1,1);
    
    % Check if partial PSD are used
    if strcmp(command.PSD_partial,'elt') || strcmp(command.PSD_partial,'dir') || strcmp(command.PSD_partial,'gen'); [~,dim]=size(ve_ps); else dim=1; end 
    ve_psx=zeros(sp,dim); ve_psex=zeros(sp,dim);
    
    % Pre-generate gaussian laser~
    h=zeros(sp); %g=(rand(sp,sp)-0.5)./100;
    for i=1:sp
        [~,u]=min(abs(t-t(i)));
        h(1:sp,i)=(2*sqrt(2*log(2))/((command.conv_str(1)./dt+u*command.conv_str(2))*...
                (sqrt(2*pi))))*exp(-4*log(2)*((1:sp)-i).^2/((command.conv_str(1)./dt+u*command.conv_str(2)).^2));
    end; clear i; % surf(log(h),'lines','none'); view(0,90); axis([1 sp 1 sp]);
   % h=h+g;
    % Convolution; loop through dim, data space, and gaussian space
    for x=1:dim
        ve_psc=zeros(1,sp); ve_psec=zeros(1,sp);
        X=[ve_ps(:,x)',zeros(1,ceil(command.gspace/2))];
        Xe=[ve_pse(:,x)',zeros(1,ceil(command.gspace/2))];
        for i=1:sp
            for j=1:sp
                ve_psc(i)=ve_psc(i)+X(j)*h(j,i);
                ve_psec(i)=ve_psec(i)+Xe(j)*h(j,i);
            end
        end; clear j i;
        ve_psx(:,x)=ve_psc'; ve_psex(:,x)=ve_psec'; clear ve_psc ve_psec;
    end; clear x;
    ve_ps=ve_psx; ve_pse=ve_psex;
end

%%
% Energy units
if strcmp(command.unit,'cm') || strcmp(command.unit,'cm-1')
    ve_en=ve_en.*mev_to_cm;
else
end

%%
% Rename
job_PSD.ve_ps=ve_ps;
job_PSD.ve_pse=ve_pse;
%
job_PSD.ve_en=ve_en;
end