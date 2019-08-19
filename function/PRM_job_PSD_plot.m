function [fig_flag] = PRM_job_PSD_plot(command,job_PSD,fig_flag)

%% Plot Q=f(E) dependence
f=figure(fig_flag);
set(f,'position',[100 100 500 400]);
hold on
plot(job_PSD.ve_en,job_PSD.ve_qx,'-k','linewidth',2);
plot(job_PSD.ve_en,job_PSD.ve_qxmin,'--k','linewidth',1);
plot(job_PSD.ve_en,job_PSD.ve_qxmax,'--k','linewidth',1);
hold off
axis square; 
xlim([0 800]);
ylim([0 20]);
xlabel('$\hbar\omega$~[meV]','Interpreter','latex'); 
ylabel('$Q$~[1/\AA]','Interpreter','latex');
box on;
set(gcf,'color','w');
set(gca,'fontsize',18);
fig_flag=fig_flag+1;

%% Plot DW=f(E)
f=figure(fig_flag);
set(f,'position',[100 100 500 400]);
A0=mean(mean(mean(mean(mean(job_PSD.ve_dw(:,1,1,1,:,1),6),5),4),3),2);
plot(job_PSD.ve_en,A0,'-k','linewidth',2);
axis square;
xlim([0 800]);
ylim([0 3]);
xlabel('$\hbar\omega$~[meV]','Interpreter','latex'); 
ylabel('$Q^{2}<u^{2}>$','Interpreter','latex');
box on;
set(gcf,'color','w');
set(gca,'fontsize',18);
l=legend({sprintf('$<u^{2}> = %.4f(%2.0f)$~[%s$^{2}$]',mean(job_PSD.ve_u2at),ceil(std(job_PSD.ve_u2at)*10000),'\AA')});
set(l,'location','northwest','box','off','interpreter','latex');
fig_flag=fig_flag+1;

%% Plot ~Tpn=f(E)
f=figure(fig_flag);
set(f,'position',[100 100 900 400]);
at=1; % sum over nslep, nrun, dir for atom at
A1=sum(sum(sum(job_PSD.ve_tp(:,:,1,:,at,:,1),6),4),2);
A2=sum(sum(sum(job_PSD.ve_tp(:,:,1,:,at,:,2),6),4),2);
A3=sum(sum(sum(job_PSD.ve_tp(:,:,1,:,at,:,3),6),4),2);
A4=sum(sum(sum(job_PSD.ve_tp(:,:,1,:,at,:,4),6),4),2);
A5=sum(sum(sum(job_PSD.ve_tp(:,:,1,:,at,:,5),6),4),2);
A6=sum(sum(sum(job_PSD.ve_tp(:,:,1,:,at,:,6),6),4),2);
A7=sum(sum(sum(job_PSD.ve_tp(:,:,1,:,at,:,7),6),4),2);
hold on;  
plot(job_PSD.ve_en,A1,'-k','linewidth',2);
plot(job_PSD.ve_en,A2,'-r','linewidth',2);
plot(job_PSD.ve_en,A3,'-b','linewidth',2);
plot(job_PSD.ve_en,A4,'-g','linewidth',2);
plot(job_PSD.ve_en,A5,'--r','linewidth',2);
plot(job_PSD.ve_en,A6,'--b','linewidth',2);
plot(job_PSD.ve_en,A7,'--g','linewidth',2);
hold off
xlim([0 800]);
box on;
set(gcf,'color','w');
set(gca,'fontsize',18);
set(gca,'YTickLabel',[]);
xlabel('$\hbar\omega$~[meV]','Interpreter','latex');
ylabel('$\tilde{T_{p}}$','Interpreter','latex');
fig_flag=fig_flag+1;

%% Plot Job4 PSD results
f=figure(fig_flag);
set(f,'Position',[100,100,900,400]);

if strcmp(command.PSD_partial,'elt')
    command.elt_tmp=cell(length(command.elt),1);
    for k=1:length(command.elt); command.elt_tmp{k,1}=command.elt{1,k}{1,1}; end; clear k;
    [~,idx]=unique(command.elt_tmp,'first'); command.elt_rn=command.elt_tmp(sort(idx));
    [~,dim]=size(job_PSD.ve_ps); legendinfo=cell(1,dim-1);
    hold on
    plot(job_PSD.ve_en,job_PSD.ve_ps(:,1),'-k','linewidth',2);
    for k=1:dim-1
        plot(job_PSD.ve_en,job_PSD.ve_ps(:,k+1));
        legendinfo{k}=command.elt_rn{k,1};
    end; clear k;
    hold off
    legend(['iso',legendinfo]);
    
elseif strcmp(command.PSD_partial,'dir')
    hold on
    plot(job_PSD.ve_en,job_PSD.ve_ps(:,2),'-b','linewidth',2);
    plot(job_PSD.ve_en,job_PSD.ve_ps(:,3),'-g','linewidth',2);
    plot(job_PSD.ve_en,job_PSD.ve_ps(:,4),'-r','linewidth',2);
    plot(job_PSD.ve_en,job_PSD.ve_ps(:,1),'-ok','linewidth',2);
    hold off
    legend('x','y','z','iso');
    
elseif strcmp(command.PSD_partial,'gen')
    [~,dim]=size(job_PSD.ve_ps);
    hold on
    plot(job_PSD.ve_en,job_PSD.ve_ps(:,1),'-k','linewidth',2);
    for x=2:dim
        plot(job_PSD.ve_en,job_PSD.ve_ps(:,x),'-');
    end; clear x;
    hold off
    
else
    hold on
    Ac=job_PSD.ve_pc(:,1)./trapz(job_PSD.ve_pc(:,1));
    As=job_PSD.ve_ps(:,1)./trapz(job_PSD.ve_ps(:,1));
    plot(job_PSD.ve_en,As,'-k','linewidth',2);
    plot(job_PSD.ve_en,Ac,'-r','linewidth',2);
    hold off
    AM=max([max(As),max(Ac)]);
    ylim([0 1.05*AM]);
    l=legend({'$G(\omega)$','$S(Q,\omega)$'});
    set(l,'box','off','Interpreter','latex');
end

if strcmp(command.unit,'cm') || strcmp(command.unit,'cm-1')
    xlim([0 5000]);
    xlabel('$\nu$~[cm$^{-1}$]','Interpreter','latex');
else
    xlim([0 600]);
    xlabel('$\hbar\omega$~[meV]','Interpreter','latex');
end

ylabel('Intensity~[arb. units]','Interpreter','latex');
set(gca,'YTickLabel',[]);
set(gca,'fontsize',18);

fig_flag=fig_flag+1;
set(gcf,'color','w');
box on;
clear f_k;
    
end

