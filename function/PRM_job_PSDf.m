function [job_PSDf] = PRM_job_PSDf(command,structure,job_OHG)

% PSD on filtered velocities
    % Extract trajectories associated with specific HB distances as defined
    % by doh_voidP, using doh_abc (job OHG). Fourier transform velocities
    % using the FFT algorithm. Numerical noise is inherent to data series
    % with gap and cannot be fixed from the algorithm (nor any). As we
    % consider scenario with up to 90% of gap, FFT is prefered over
    % DFT/Z/LS as its results are more consistent at high gap level. Gaps
    % are filed with zeros, and each fragments of the data series is
    % renormalized arround zero. Each fragment can also be weighted by a
    % window function, but from my test is seems only to make the results
    % worse. Theoretically, use the slepian window for optimized peak
    % density, the kaiser window for a good overall attenuation, or the
    % tukey window to maximize data usage. Similarly, the renormalization
    % of each fragment can be done by 'mean' or by 'polyfit' (default); it
    % has only little effect on the results. Compared to 'job PSD' where a
    % proper exact set of parameters can be defined, here the gaps dominate
    % the results, and little can be done about it. End results are quite
    % noisy and can be smoothed (command.PSDf_filter) or convoluted
    % (command.conv_str).

% Build limits from doh_voidP
% X=job_OHG.doh_voidP(2); S=job_OHG.doh_voidP(3);
% PSDf_lim=[X-1.5*S,X-0.5*S,X+0.5*S,X+1.5*S]; clear X S;
% PSDf_lim=[1.725,1.875,2.025,2.175];
% PSDf_lim=[1.4,1.6,1.8,2.0];
PSDf_lim=[2.6,2.7,2.8,2.9];

% PSDf_lim as logical table
[len,command.nrun,n_H,~]=size(job_OHG.doh_abc); frag=length(PSDf_lim)+1;
PSDf_geom=permute(job_OHG.doh_abc(2:len,:,:,:),[4 5 3 2 1]); % dOH/1/n_H/command.nrun/len
PSDf_logical=zeros([length(command.at) 3 n_H command.nrun len-1 frag]);
% PSDf_logical(:,:,:,:,:,1)=repmat(PSDf_geom(2,:,:,:,:)<=PSDf_lim(1),[length(command.at) 3 1 1 1]); %HB filter
PSDf_logical(:,:,:,:,:,1)=repmat(PSDf_geom(8,:,:,:,:)<=PSDf_lim(1),[length(command.at) 3 1 1 1]); %OO filter
for f=2:frag-1
%     PSDf_logical(:,:,:,:,:,f)=repmat(logical((PSDf_geom(2,:,:,:,:)>PSDf_lim(f-1))...
%         .*(PSDf_geom(2,:,:,:,:)<=PSDf_lim(f))),[length(command.at) 3 1 1 1]);
    PSDf_logical(:,:,:,:,:,f)=repmat(logical((PSDf_geom(8,:,:,:,:)>PSDf_lim(f-1))...
        .*(PSDf_geom(8,:,:,:,:)<=PSDf_lim(f))),[length(command.at) 3 1 1 1]);
end; clear f;
PSDf_logical(:,:,:,:,:,frag)=repmat(PSDf_geom(2,:,:,:,:)>PSDf_lim(frag-1),[length(command.at) 3 1 1 1]);
PSDf_logical=logical(PSDf_logical);

%     permute((sum(PSDf_logical(x,1,x,:,:,:),5)),[4 6 5 1 2 3])

% 	permute((PSDf_logical(1,1,1,:,:,3)),[5 4 6 1 2 3])+...
%         permute((PSDf_logical(2,1,2,:,:,3)),[5 4 6 1 2 3])+...
%         permute((PSDf_logical(3,1,3,:,:,3)),[5 4 6 1 2 3])+...
%         permute((PSDf_logical(4,1,4,:,:,3)),[5 4 6 1 2 3])

%     permute(sum((permute((PSDf_logical(1,1,1,:,:,:)),[5 4 6 1 2 3])+...
%         permute((PSDf_logical(2,1,2,:,:,:)),[5 4 6 1 2 3])+...
%         permute((PSDf_logical(3,1,3,:,:,:)),[5 4 6 1 2 3])+...
%         permute((PSDf_logical(4,1,4,:,:,:)),[5 4 6 1 2 3]))==1,1),[2 3 1])

% Compute filtered velocities
PSDf_ve_tmp=repmat(permute((structure.pos(3:len+1,:,:,:)-structure.pos(2:len,:,:,:))./command.time_step,[3 2 5 4 1]),[1 1 n_H 1 1]);
PSDf_ve=zeros(size(PSDf_logical));
for f=1:frag
    for n=1:command.nrun
        for h=1:n_H
            for y=1:3
                for x=1:length(command.at)
                    PSDf_ve(x,y,h,n,:,f)=PSDf_ve_tmp(x,y,1,n,:).*PSDf_logical(1,1,h,n,:,f);
                end
            end
        end
    end
end; clear x y h n f;


% Reshape velocities into nseg segments of length command.cut
command.cut=command.cut/command.time_step;
vseg=1:command.cut:len; nseg=sum((vseg+command.cut-1)<=len);
PSDf_ve_red=zeros([nseg,3,n_H,command.nrun,command.cut,frag]);
for f=1:frag
    for x=1:nseg
        for h=1:n_H
            PSDf_ve_red(x,:,h,:,:,f)=PSDf_ve(h,:,h,:,vseg(x):vseg(x)+command.cut-1,f);
        end
    end
end; clear h x f;
PSDf_ve_red=permute(PSDf_ve_red,[5 2 3 4 1 6]); % command.cut,3,n_H,nrun,nseg,frag


% Data correction: correct velocities from local 'drift'
% for f=1:frag
%     for x=1:nseg
%         for n=1:command.nrun
%             for h=1:n_H
%                 for y=1:3
%                     vt=PSDf_ve_red(:,y,h,n,x,f);
%                     A=logical(vt(2:length(vt)))-logical(vt(1:length(vt)-1));
%                     B1=find(A==1)+1; B2=find(A==-1);
%                     if length(B2)>length(B1); B1=cat(1,1,B1); end
%                     if length(B2)<length(B1); B2=cat(1,B2,length(vt)); end
%                     B=cat(2,B1,B2); [Blen,~]=size(B); % index and number of fragments
%                     for k=1:Blen
%                         % Per fragment, shifting velocities to <v>=0
%                         C=vt(B(k,1):B(k,2),1); Dx=(1:1:length(C))';
%                         if length(Dx)>21 % fit ax+b % min 3 for consistency, can be increase up to 21 without major changes
%                             Cn=polyfit(Dx,C,1); D=C-(Dx.*Cn(1)+Cn(2));
%                         else
% %                             D=C-mean(C,1); % cst
%                             D=C.*0;
%                         end
%                         PSDf_ve_red(B(k,1):B(k,2),y,h,n,x,f)=D;
%                     end
%                     clear vt A B1 B2 B Blen C Cn D Dx;
%                 end
%             end
%         end
%     end
% end; clear k y h n x f;


% FFT
PSDf_fft=zeros([command.cut,3,n_H,command.nrun,nseg,frag]);
for f=1:frag
    for n=1:command.nrun
        for y=1:3
            for h=1:n_H
                for x=1:nseg
                    PSDf_fft(:,y,h,n,x,f)=fft(PSDf_ve_red(:,y,h,n,x,f),command.cut);
                end
            end
        end
    end
end; clear x h y n f A;

PSDf_psum=permute(sum(sum(sum(sum(power(abs(PSDf_fft),2),5),4),3),2),[1 6 2 3 4 5]);
PSDf_ps=PSDf_psum(1:command.cut/2+1,:);

thz_to_mev=4.135665538536;
PSDf_en=((((1/command.time_step*10^3)/2*linspace(0,1,command.cut/2+1)).*thz_to_mev)');
t=(((1/command.time_step*10^3)/2*linspace(0,2,command.cut+1)).*thz_to_mev)'; % pre-built for convolution


% Energy unit
mev_to_cm=8.06554;
if strcmp(command.unit,'cm') || strcmp(command.unit,'cm-1')
    PSDf_en=PSDf_en.*mev_to_cm;
else
end


% Optional smoothing
if command.PSDf_filter==1
    for f=1:frag
        PSDf_ps(:,f)=smooth(PSDf_ps(:,f),9);%9
    end; clear f;
else
end


% Optional convolution
if (((command.conv_str(1)~=0)+(command.conv_str(2)~=0))*(command.PSDf_filter==0))==1

    % Generate extended energy domain
    command.gspace=199; sp=length(PSDf_en)+floor(command.gspace/2);
    PSDf_en=t(1:sp,1); dt=PSDf_en(2,1)-PSDf_en(1,1);
    PSDf_psx=zeros(sp,frag);

    % Pre-generate gaussian laser~
    h=zeros(sp);
    for i=1:sp
        [~,u]=min(abs(t-t(i)));
        h(1:sp,i)=(2*sqrt(2*log(2))/((command.conv_str(1)./dt+u*command.conv_str(2))*...
                (sqrt(2*pi))))*exp(-4*log(2)*((1:sp)-i).^2/((command.conv_str(1)./dt+u*command.conv_str(2)).^2));
    end; clear i u;

    % Convolution; loop through frag, data space, and gaussian space
    for x=1:frag
        PSDf_psc=zeros(1,sp);
        X=[PSDf_ps(:,x)',zeros(1,ceil(command.gspace/2))];
        for i=1:sp
            for j=1:sp
                PSDf_psc(i)=PSDf_psc(i)+X(j)*h(j,i);
            end
        end
        PSDf_psx(:,x)=PSDf_psc'; clear PSDf_psc;
    end; clear j i x;
    PSDf_ps=PSDf_psx;
else
end

% Rename
job_PSDf.PSDf_en=PSDf_en;
job_PSDf.PSDf_ps=PSDf_ps;
job_PSDf.PSDf_lim=PSDf_lim;

end

