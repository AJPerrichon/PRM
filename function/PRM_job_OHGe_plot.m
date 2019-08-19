function [fig_flag] = PRM_job_OHGe_plot(command,job_OHG,fig_flag)

% Plot job 7 OHG 
% only plot command.nrun=1 for clarity (yet average values are using all selected generations

% % Plot distance-related figures
% [len,~,n_H,~]=size(job_OHG.doh_abc); h=0;
% t=1:1:len; t=(t.*command.time_step)';
% while h<n_H
% f_m=figure(fig_flag);
% set(f_m,'Position',[80,15,1601,801]);
% %
% subplot(2,2,1);
% plot(t,job_OHG.doh_abc(:,1,h+1,1),'-k');
% xlim([1000 1500]); %xlim([0 max(t)]);
% xlabel('Time (fs)');
% ylabel('Distance (Å)');
% title(sprintf('O-H distance; <d_{O-H}>=%1.4f Å',job_OHG.doh_av(1,h+1)));
% set(gca,'fontsize',20);
% box on;
% %
% subplot(2,2,2);
% hold on
% plot(t,job_OHG.doh_abc(:,1,h+1,2),'-r');
% plot(t,job_OHG.doh_abc(:,1,h+1,3),'-b');
% hold off
% xlim([1000 1500]); %xlim([0 max(t)]);
% xlabel('Time (fs)');
% ylabel('Distance (Å)');
% legend('closest second neighbor','next closest second neighbor','location','north');
% title(sprintf('Sorted H- -O distances; <d_{H- -O}>=%1.4f and %1.4f Å',job_OHG.doh_av(2,h+1),job_OHG.doh_av(3,h+1)));
% set(gca,'fontsize',20);
% box on;
% %
% subplot(2,2,3);
% plot(job_OHG.doh_en,job_OHG.doh_ps(:,1,h+1),'-k','linewidth',2);
% if strcmp(command.unit,'cm') || strcmp(command.unit,'cm-1')
%     xlim([1 4000]); xlabel('Frequency (cm^{-1})');
% else
%     xlim([1 600]); xlabel('Frequency (meV)');
% end
% ylim([0 1.1*max(max(job_OHG.doh_ps(1:end,1,h+1)))]); ylabel('Density (a.u.)');
% title(sprintf('FFT of O-H'));
% set(gca,'fontsize',20);
% box on;
% %
% subplot(2,2,4);
% hold on
% plot(job_OHG.doh_en,job_OHG.doh_ps(:,2,h+1),'-r','linewidth',2);
% plot(job_OHG.doh_en,job_OHG.doh_ps(:,3,h+1),'-b','linewidth',2);
% hold off
% if strcmp(command.unit,'cm') || strcmp(command.unit,'cm-1')
%     xlim([1 4000]); xlabel('Frequency (cm^{-1})');
% else
%     xlim([1 600]); xlabel('Frequency (meV)');
% end
% ylim([0 1.1*max(max(job_OHG.doh_ps(1:end,2:3,h+1)))]); ylabel('Density (a.u.)');
% legend('closest second neighbor','next closest second neighbor','location','north');
% title(sprintf('FFT of H- -O'));
% set(gca,'fontsize',20);
% box on;
% 
% fig_flag=fig_flag+1;
% set(gcf,'color','w');
% clear f_m;
% h=h+1;
% end

%%
X=job_OHG.doh_x'; Xstep=mean(X(2:end)-X(1:end-1));
Y1=job_OHG.doh_void_1st;
Y2=job_OHG.doh_void_2nd;
g=@(x,t) exp(-0.5.*(((t-x(1))./x(2)).^2))./sqrt(2.*pi)./x(2).*Xstep; % fitting a gaussian
[P1,~,~,~,~]=lsqcurvefit(g,[1.0 0.1],X,Y1);
[P2,~,~,~,~]=lsqcurvefit(g,[1.9 0.2],X,Y2);

% Distribution of distances plot
f_n=figure(fig_flag);
set(f_n,'Position',[80,15,501,501]);
hold on
plot(X,Y1,'-ok','linewidth',2);
plot(X,g(P1,X),'-r','linewidth',2);
hold off
xlim([0.5 4.0]); ylim([0 1.15*max(Y1)]);
set(gca,'fontsize',14,'YTick',[],'XTick',1:1:4,'XMinorTick','on');
xlabel('Distance~[\AA]','interpreter','latex');
ylabel('Counts','interpreter','latex');
title('Histogram of distances (1st neighbor)');
text(0.05*3.5+0.5,1.08*max(Y1),sprintf('FWHM=%3.3f~%s',P1(2).*2*sqrt(2.*log(2)),'\AA'),'interpreter','latex','fontsize',14);
text(0.5*3.5+0.5,1.08*max(Y1),sprintf('x$_0$=%3.3f~%s',P1(1),'\AA'),'interpreter','latex','fontsize',14);
box on;
fig_flag=fig_flag+1;
set(gcf,'color','w');
clear f_n;

% Distribution of distances plot
f_n=figure(fig_flag);
set(f_n,'Position',[80,15,501,501]);
hold on
plot(X,Y2,'-ok','linewidth',2);
plot(X,g(P2,X),'-r','linewidth',2);
hold off
xlim([0.5 4]); ylim([0 1.15*max(Y2)]);
set(gca,'fontsize',14,'YTick',[],'XTick',1:1:4,'XMinorTick','on');
xlabel('Distance~[\AA]','interpreter','latex');
ylabel('Counts','interpreter','latex');
title('Histogram of distances (2nd neighbor)');
text(0.05*3.5+0.5,1.08*max(Y2),sprintf('FWHM=%3.3f~%s',P2(2).*2*sqrt(2*log(2)),'\AA'),'interpreter','latex','fontsize',14);
text(0.5*3.5+0.5,1.08*max(Y2),sprintf('x$_0$=%3.3f~%s',P2(1),'\AA'),'interpreter','latex','fontsize',14);
box on;
fig_flag=fig_flag+1;
set(gcf,'color','w');
clear f_n;

% Correlation distance (1,2) and (8)
X_HB=(0:0.01:3)';
X_OO=(2:0.01:4)';
[len,nrun,nat,~]=size(job_OHG.doh_abc);
corr_void_CB=zeros([length(X_HB),length(X_OO)]);
corr_void_HB=zeros([length(X_HB),length(X_OO)]);
for n=1:nrun
    for x=1:len
        for h=1:nat
            [~,iCB]=min(abs(X_HB-job_OHG.doh_abc(x,n,h,1))); % covalent bond
            [~,iHB]=min(abs(X_HB-job_OHG.doh_abc(x,n,h,2))); % hydrogen bond
            [~,iOO]=min(abs(X_OO-job_OHG.doh_abc(x,n,h,8))); % O-O distance
            corr_void_CB(iCB,iOO)=corr_void_CB(iCB,iOO)+1;
            corr_void_HB(iHB,iOO)=corr_void_HB(iHB,iOO)+1;
        end
    end
end; clear x h n;
g=figure;
set(g,'Position',[100 100 500 500])
hold on
surf(X_OO,X_HB,log10(corr_void_CB),'lines','none'); view(0,90); colormap jet;
surf(X_OO,X_HB,log10(corr_void_HB),'lines','none'); view(0,90); colormap jet;
l1=line([0 8],[0 4]); set(l1,'linewidth',1,'linestyle','--','color','k');
% l2=line([2.40 2.40],[0 4]); set(l2,'linewidth',1,'linestyle','--','color','k');
l3=line([2.6 2.6],[0 4]); set(l3,'linewidth',1,'linestyle','-','color','k');
hold off
axis square;
xlim([2.3 3.5]); ylim([0.8 2.8]);
set(gca,'fontsize',14,...
    'XTick',2.4:.2:3.4,'XTickLabel',{'2.4','2.6','2.8','3.0','3.2','3.4'},'XMinorTick','on',...
    'YTick',1.0:.4:2.8,'YTickLabel',{'1.0','1.4','1.8','2.2','2.6'},'YMinorTick','on');
xlabel('O-O distance~[\AA]','Interpreter','Latex');
ylabel('H covalent/hydrogen bond length~[\AA]','Interpreter','Latex');
% text(3.25,1.75,'(2$c$)','Interpreter','Latex','Fontsize',16);
text(2.44,2.65,'(A)','Interpreter','Latex','Fontsize',16);
text(2.66,2.65,'(B)','Interpreter','Latex','Fontsize',16);
box on;
set(gcf,'color','w');
% set(g,'RendererMode','man');
% set(g,'Renderer','painters');
% print(g,'S2','-depsc');
clear g l1 l3;

%%
if command.OHG_demod==1
    % Demodulation and PSD plots
    f_n=figure(fig_flag); set(f_n,'Position',[80,15,1601,801]);
    k=1;
    subplot(4,2,1);
    hold on
    plot((1:1:len)'.*command.time_step,job_OHG.doh_abc(:,1,k,1),'-b','linewidth',2);
    plot((1:1:len)'.*command.time_step,job_OHG.doh_demod(:,1,k,1),'-r','linewidth',2);
    plot((1:1:len)'.*command.time_step,job_OHG.doh_demod(:,1,k,2),'-r','linewidth',2);
    hold off
    title('Raw signal'); set(gca,'fontsize',18); set(gca,'XTickLabel',[]); 
    box on; set(gcf,'color','w'); xlim([1000 2000]); ylim([min(job_OHG.doh_abc(:,1,k,1))-0.002 max(job_OHG.doh_abc(:,1,k,1))+0.002]);
    %
    subplot(4,2,3);
    plot((1:1:len)'.*command.time_step,job_OHG.doh_demod(:,1,k,4),'-g','linewidth',2);
    title('Demodulation: Envelope center'); set(gca,'fontsize',18); set(gca,'XTickLabel',[]); 
    box on; set(gcf,'color','w'); xlim([1000 2000]); ylim([min(job_OHG.doh_demod(:,1,k,4))-0.002 max(job_OHG.doh_demod(:,1,k,4))+0.002]);
    %
    subplot(4,2,5);
    plot((1:1:len)'.*command.time_step,job_OHG.doh_demod(:,1,k,3),'-k','linewidth',2);
    title('Demodulation: Envelope amplitude'); set(gca,'fontsize',18); set(gca,'XTickLabel',[]); 
    box on; set(gcf,'color','w'); xlim([1000 2000]); ylim([0 max(job_OHG.doh_demod(:,1,k,3)).*1.05]);
    %
    subplot(4,2,7);
    plot((1:1:len)'.*command.time_step,job_OHG.doh_demod(:,1,k,5),'-r','linewidth',2);
    title('Demodulation: Pure Carrier'); xlabel('time (fs)'); set(gca,'fontsize',18);
    box on; set(gcf,'color','w'); xlim([1000 2000]); ylim([-1.05 1.05]);
    %
    subplot(4,2,2);
    plot(job_OHG.doh_en,job_OHG.doh_ps(:,1,k),'-b','linewidth',2);
    title('PSD: Raw signal'); set(gca,'YTickLabel',[]); 
    set(gca,'XTickLabel',[]); ylabel('P(\omega)'); set(gca,'fontsize',18);
    box on; set(gcf,'color','w'); xlim([0 500]); ylim([0 max(job_OHG.doh_ps(:,1,k)).*1.05]);
    %
    subplot(4,2,6);
    plot(job_OHG.doh_en,job_OHG.doh_demod_ps(:,1,k),'-k','linewidth',2);
    title('PSD: Envelope amplitude'); set(gca,'YTickLabel',[]); 
    set(gca,'XTickLabel',[]); ylabel('P(\omega)'); set(gca,'fontsize',18);
    box on; set(gcf,'color','w'); xlim([0 500]); ylim([0 max(job_OHG.doh_demod_ps(:,1,k)).*1.05]);
    %
    subplot(4,2,4);
    plot(job_OHG.doh_en,job_OHG.doh_demod_ps(:,2,k),'-g','linewidth',2);
    title('PSD: Envelope center'); set(gca,'YTickLabel',[]); 
    set(gca,'XTickLabel',[]); ylabel('P(\omega)'); set(gca,'fontsize',18);
    box on; set(gcf,'color','w'); xlim([0 500]); ylim([0 max(job_OHG.doh_demod_ps(:,2,k)).*1.05]);
    %
    subplot(4,2,8);
    plot(job_OHG.doh_en,job_OHG.doh_demod_ps(:,3,k),'-r','linewidth',2);
    title('PSD: Pure Carrier'); set(gca,'YTickLabel',[]); 
    xlabel('\omega (meV)'); ylabel('P(\omega)'); set(gca,'fontsize',18);
    box on; set(gcf,'color','w'); xlim([0 500]); ylim([0 max(job_OHG.doh_demod_ps(:,3,k)).*1.05]);
    %
    clear f_n; fig_flag=fig_flag+1;
    
    % Time Lag versus distance plot
    f_n=figure(fig_flag); set(f_n,'Position',[80,15,880,815]);
    T=((command.time_step:command.time_step:100).*2)';
    D=(min(min(min(job_OHG.doh_abc(:,:,:,1)))):0.001:max(max(max(job_OHG.doh_abc(:,:,:,1)))))';
    Te=((command.time_step:command.time_step/10:100).*2)';
    De=(min(min(min(job_OHG.doh_abc(:,:,:,1)))):0.001/10:max(max(max(job_OHG.doh_abc(:,:,:,1)))))';
    [DD,TT]=meshgrid(D,T);
    [DDe,TTe]=meshgrid(De,Te);
    M=log10(job_OHG.doh_zc);
    Me=interpn(TT,DD,M,TTe,DDe);
    surf(De,Te,Me,'lines','none'); view(0,90);
    set(gcf,'color','w'); box on; 
    xlim([0.95 1.1]); ylim([0 30]);
    xlabel('d_{O-H} (Å)'); ylabel('Stretch period (fs)'); 
    set(gca,'fontname','times new roman','fontsize',20);
    clear f_n; fig_flag=fig_flag+1;
    
else
end

end

