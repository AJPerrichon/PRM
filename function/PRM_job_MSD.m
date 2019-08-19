function [job_MSD] = PRM_job_MSD(command,structure)

[len,~,~,~]=size(structure.pos);

% Isolate relevant displacements % Meaning finding at each step the closest
% O from considered H and substracting its position, assuming a general
% case where H is diffusing and moving from O to O (with O that can be
% images from periodic conditions).
if command.MSD_sub==1
    
    % Cell to supercell expansion
    pos_sc=zeros(len,3,length(command.at),command.nrun,27);
    t=1;
    for h=[0,-1,1] % one extra cell on each side
        for k=[0,-1,1] % for all 3 dimensions
            for l=[0,-1,1] % equals 3^3=27 cells
                for m=1:length(command.at)
                    for n=1:command.nrun
                        pos_sc(:,:,m,n,t)=cat(2,structure.pos(:,1,m,n)+h,structure.pos(:,2,m,n)+k,structure.pos(:,3,m,n)+l);
                    end; clear n;
                end; clear m;
                t=t+1;
            end; clear l;
        end; clear k;
    end; clear h t;
    
    % Count atoms
    n_H=zeros(length(command.at),1);
    for k=1:length(command.at)
        n_H(k,1)=strcmp(command.elt{1,k},'H');
    end; clear k;
    n_H=sum(n_H);

    n_O=zeros(length(command.at),1);
    for k=1:length(command.at)
        n_O(k,1)=strcmp(command.elt{1,k},'O');
    end; clear k;
    n_O=sum(n_O);
    
    % Supercell to cell reduction
    pos_red=zeros(len,3,command.nrun,n_H,n_O+1);
    for n=1:command.nrun
        for h=1:n_H
            pos_red(:,:,n,h,1)=pos_sc(:,:,h,n,1); % k=1 is H,
            for k=1:n_O % for each O and M k=2:end
                dsc_tag=zeros(27,1);
                for y=1:27 % search through all cells
                    dsc_tag(y,1)=norm(cat(2,...
                        structure.cell_parameters(1)*(pos_sc(1,1,k+n_H,n,y)-pos_sc(1,1,h,n,1)),...
                        structure.cell_parameters(2)*(pos_sc(1,2,k+n_H,n,y)-pos_sc(1,2,h,n,1)),...
                        structure.cell_parameters(3)*(pos_sc(1,3,k+n_H,n,y)-pos_sc(1,3,h,n,1))));
                end; clear y;
                [~,tt]=min(dsc_tag); % which distance is shortest to H
                pos_red(:,:,n,h,k+1)=pos_sc(:,:,k+n_H,n,tt); % and extract it
            end; clear k dsc_tag tt;
        end; clear h;
    end; clear n pos_sc;
    
    pos_true=zeros(len-1,3,n_H,command.nrun);
    
    doh_tmp=zeros(len,command.nrun,n_H,n_O);
    for n=1:command.nrun
        for h=1:n_H
        % Calculate all H-O distances (will be used to sort O atoms)
            for k=2:n_O+1
                for x=1:len
                    doh_tmp(x,n,h,k-1)=norm(cat(2,...
                        (pos_red(x,1,n,h,k)-pos_red(x,1,n,h,1)).*structure.cell_parameters(1),...
                        (pos_red(x,2,n,h,k)-pos_red(x,2,n,h,1)).*structure.cell_parameters(2),...
                        (pos_red(x,3,n,h,k)-pos_red(x,3,n,h,1)).*structure.cell_parameters(3)));
                end; clear x;
            end; clear k;
            for x=1:len-1
                [~,I]=sort(permute(doh_tmp(x+1,n,h,:),[4 3 2 1]),1);
                pos_true(x,:,h,n)=pos_red(x+1,:,n,h,I(1)+1)-pos_red(x+1,:,n,h,1);
            end; clear x;
        end; clear h;
    end; clear n;
    
    structure.pos=pos_true;
    [len,~,~,~]=size(structure.pos);
    len_at=n_H;
    clear pos_true n_H n_O I doh_tmp;
else
    len_at=length(command.at);
end

% Calculation of the MSD

msd_void=zeros(len_at,command.nrun,len-1);

for n=1:len-1 % time lag
    
    % Compute displacements with lag
    msd_tmp=zeros(len-n,len_at,command.nrun,3);
    for h=1:len_at
        for k=1:command.nrun
            msd_tmp(1:len-n,h,k,:)=(structure.pos(1+n:len,:,h,k)-structure.pos(1:len-n,:,h,k)).^2;
        end; clear k;
    end; clear h;
    
    % Convert r.l.u. to angstroms
    msd_tmp(:,:,:,1)=msd_tmp(:,:,:,1).*(structure.cell_parameters(1)^2);
    msd_tmp(:,:,:,2)=msd_tmp(:,:,:,2).*(structure.cell_parameters(2)^2);
    msd_tmp(:,:,:,3)=msd_tmp(:,:,:,3).*(structure.cell_parameters(3)^2);
    msd_tmp=sum(msd_tmp,4);
    
    % MSD
    for h=1:len_at
        for k=1:command.nrun
            msd_void(h,k,n)=1/(len-n)*sum(msd_tmp(:,h,k));
        end; clear k;
    end; clear h;
    
    clear msd_tmp;
end; clear n;

msd_true=permute(msd_void,[3 1 2]);

% Time vector
msd_time=(command.time_step:command.time_step:(len-1)*command.time_step)';

% Temporary plot
figure;
% hold on
%     plot(msd_time(:,1)./1000,msd_true(:,1,1),'-b','linewidth',2); % per gen
%     plot(msd_time(:,1)./1000,msd_true(:,1,2),'-k','linewidth',2);
%     plot(msd_time(:,1)./1000,msd_true(:,1,3),'-r','linewidth',2);
%     plot(msd_time(:,1)./1000,msd_true(:,1,4),'-g','linewidth',2);
% hold off
% hold on
plot(msd_time(:,1)./1000,mean(msd_true(:,1,:),3),'-k','linewidth',2);
% plot(msd_time_50(:,1)./1000,mean(msd_true_50(:,1,:),3),'-r','linewidth',2);
% hold off
box on; set(gcf,'color','w');
xlabel('Time (ps)'); ylabel('MSD (arb. units)'); set(gca,'fontsize',20)
% xlim([0 10]); ylim([0 80]);
% legend('100% hydrated','50% hydrated');

% Clean up of the lattice part of the MSD using a pass-band filter

if command.MSD_filter==1
    msd_tmp=zeros(size(msd_true));
    for h=1:len_at
        for k=1:command.nrun
            msd_tmp(:,h,k)=smooth(msd_true(:,h,k),74);
        end; clear k;
    end; clear h;
    
    thz_to_mev=4.135665538536;
    msd_fft_highpass=fft(msd_true-msd_tmp); msd_psd_highpass=power(abs(msd_fft_highpass(1:(len-1)/2,1,:)),2);
    msd_fft_lowpass=fft(msd_tmp); msd_psd_lowpass=power(abs(msd_fft_lowpass(1:(len-1)/2,1,:)),2);
    msd_fft_full=fft(msd_true); msd_psd_full=power(abs(msd_fft_full(1:(len-1)/2,1,:)),2);
    msd_en=((1/command.time_step*10^3)/2*linspace(0,1,(len-1)/2).*thz_to_mev)';
else
end

if command.MSD_filter==1
    figure;
    hold on
        plot(msd_en,sum(msd_psd_highpass(:,1,:),3),'-r','linewidth',2);
        plot(msd_en,sum(msd_psd_lowpass(:,1,:),3),'-b','linewidth',2);
        plot(msd_en,sum(msd_psd_full(:,1,:),3),'-k','linewidth',2);
    hold off
    legend('high-pass MSD','low-pass MSD','full MSD');
    xlim([0 150]); ylim([0 2500]);

    figure; %msd_true is full MSD; msd_tmp is filtered MSD
    hold on
        plot(msd_time(:,1),msd_true(:,1,1),'-b'); % per gen
        plot(msd_time(:,1),msd_true(:,1,2),'-k');
        plot(msd_time(:,1),msd_true(:,1,3),'-r');
        plot(msd_time(:,1),msd_true(:,1,4),'-g');
        plot(msd_time(:,1),msd_tmp(:,1,1),'--b');
        plot(msd_time(:,1),msd_tmp(:,1,2),'--k');
        plot(msd_time(:,1),msd_tmp(:,1,3),'--r');
        plot(msd_time(:,1),msd_tmp(:,1,4),'--g');
    hold off

    figure; %difference between low pass filter and full MSD == high frequ.
    hold on
        plot(msd_time(:,1),msd_true(:,1,1)-msd_tmp(:,1,1),'-b');
        plot(msd_time(:,1),msd_true(:,1,2)-msd_tmp(:,1,2),'-k');
        plot(msd_time(:,1),msd_true(:,1,3)-msd_tmp(:,1,3),'-r');
        plot(msd_time(:,1),msd_true(:,1,4)-msd_tmp(:,1,4),'-g');
    hold off
else
end

% Rename
job_MSD.msd_true=msd_true;
job_MSD.msd_time=msd_time;

end

