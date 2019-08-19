function [fig_flag] = PRM_job_PSDf_plot(command,job_PSDf,fig_flag)

%% Equal weight plot
f_k=figure(fig_flag);
set(f_k,'Position',[80,15,480,480]);
[~,frag]=size(job_PSDf.PSDf_ps);
cc=flip(jet(frag),1);

hold on
for f=1:frag; h(f)=plot(job_PSDf.PSDf_en,job_PSDf.PSDf_ps(:,f)./sum(job_PSDf.PSDf_ps(:,f),1),'-','linewidth',2,'color',cc(f,:)); end; clear f;
plot(job_PSDf.PSDf_en,sum(job_PSDf.PSDf_ps,2)./sum(sum(job_PSDf.PSDf_ps,2),1),'-ok','linewidth',2);
hold off
%
if strcmp(command.unit,'cm') || strcmp(command.unit,'cm-1')
    xlim([0 5000]); 
    xlabel('Wavenumber (cm^{-1})');
else
    xlim([0 600]); 
    xlabel('\omega (meV)');
end
ylim auto; % ylim([0 0.01]); % ylim([0 0.03]);
ylabel('S(\omega) (arb. units)');
set(gca,'fontsize',20); 
set(gca,'fontname','Times New Roman');
set(gcf,'color','w'); box on;
%
fig_flag=fig_flag+1;
clear f_k;

%% True weight plot
% f_k=figure(fig_flag);
% set(f_k,'Position',[80,15,480,480]);
% [~,frag]=size(job_PSDf.PSDf_ps);
% cc=flip(jet(frag),1);
% EN=job_PSDf.PSDf_en;
% %
% hold on
% for f=1:frag; plot(EN,job_PSDf.PSDf_ps(:,f)./sum(sum(job_PSDf.PSDf_ps,2),1),'linestyle','-','linewidth',2,'color',cc(f,:)); end
% clear f;
% plot(EN,sum(job_PSDf.PSDf_ps,2)./sum(sum(job_PSDf.PSDf_ps,2),1),'-ok','linewidth',2);
% hold off
% %
% if strcmp(command.unit,'cm') || strcmp(command.unit,'cm-1')
%     xlim([0 5000]); 
%     xlabel('Wavenumber (cm^{-1})');
% else
%     xlim([0 600]); 
%     xlabel('\omega (meV)');
% end
% ylim auto; % ylim([0 0.01]); % ylim([0 0.03]);
% ylabel('S(\omega) (arb. units)'); 
% set(gca,'fontname','Times New Roman');
% set(gcf,'color','w'); box on;
% %
% fig_flag=fig_flag+1;
% clear f_k;

%% Area plot from corrected true weight

[~,frag]=size(job_PSDf.PSDf_ps);
EN=job_PSDf.PSDf_en;
cd_reverse=[1,0]; % 0 from short to long HB % 1 from long to short HB

% Split data into wag / stretch intervals
[~,iWmin]=min(abs(EN-30));
[~,iWmax]=min(abs(EN-220));
[~,iSmin]=min(abs(EN-230));
[~,iSmax]=min(abs(EN-550));
%
ENW=EN(iWmin:iWmax,1); lENW=length(ENW);
ENS=EN(iSmin:iSmax,1); lENS=length(ENS);
%
YYW=zeros(lENW,frag);
YYS=zeros(lENS,frag);
for f=1:frag
    YYW(:,f)=job_PSDf.PSDf_ps(iWmin:iWmax,f);
    YYS(:,f)=job_PSDf.PSDf_ps(iSmin:iSmax,f);
end
clear f iWmin iWmax iSmin iSmax;

% Evaluate background from wag/stretch interval boundaries
for f=1:frag
    PW=polyfit(ENW([1:3,lENW-2:lENW],1),YYW([1:3,lENW-2:lENW],f),1);
    YYW(:,f)=YYW(:,f)-polyval(PW,ENW);
    %
    PS=polyfit(ENS([1:3,lENS-2:lENS],1),YYS([1:3,lENS-2:lENS],f),1);
    YYS(:,f)=YYS(:,f)-polyval(PS,ENS);
end
clear f;

% Compute cumulative sum
YYWC=zeros(lENW,frag); YYSC=zeros(lENS,frag);
YYWCw=zeros(1,frag); YYSCw=zeros(1,frag);
for f=1:frag
    if cd_reverse(1)==0
        YYWC(:,f)=sum(YYW(:,1:f),2); YYWCw(1,f)=sum(YYW(:,f),1)./sum(sum(YYW,2),1);
    else
        YYWC(:,f)=sum(YYW(:,frag:-1:f),2); YYWCw(1,f)=sum(YYW(:,f),1)./sum(sum(YYW,2),1);
    end
    if cd_reverse(2)==0
        YYSC(:,f)=sum(YYS(:,1:f),2); YYSCw(1,f)=sum(YYS(:,f),1)./sum(sum(YYS,2),1);
    else
        YYSC(:,f)=sum(YYS(:,frag:-1:f),2); YYSCw(1,f)=sum(YYS(:,f),1)./sum(sum(YYS,2),1);
    end
end
clear f;

% INS
% load data_INS;
% [~,iWmin]=min(abs(INS(:,1)-73));%75
% [~,iWmax]=min(abs(INS(:,1)-220));%170
% [~,iSmin]=min(abs(INS(:,1)-230));%340
% [~,iSmax]=min(abs(INS(:,1)-550));%470
% ENW_INS=INS(iWmin:iWmax,1); 
% ENS_INS=INS(iSmin:iSmax,1);
% %
% YYW_INS=INS(iWmin:iWmax,2);
% % YYW_INS=YYW_INS-min(YYW_INS);
% YYW_INS=YYW_INS./max(YYW_INS);
% %
% YYS_INS=INS(iSmin:iSmax,2); 
% % YYS_INS=YYS_INS-min(YYS_INS); 
% YYS_INS=YYS_INS./max(YYS_INS);
% clear iWmin iWmax iSmin iSmax;

% Area plot
cc=flip(jet(frag),1);
%
f_k=figure(fig_flag); % wag interval
hold on
ENWX=[ENW',fliplr(ENW')];
if cd_reverse(1)==0
    fw(1)=area(ENW,YYWC(:,1)./sum(sum(YYW,2),1));
    set(fw(1),'linestyle','-','linewidth',1,'edgecolor','k','facecolor',cc(1,:));
    for f=2:frag
        YYWCX=[(YYWC(:,f-1)'),fliplr((YYWC(:,f)'))]./sum(sum(YYW,2),1);
        fw(f)=fill(ENWX,YYWCX,'b');
        set(fw(f),'linestyle','-','linewidth',1,'edgecolor','k','facecolor',cc(f,:));
    end
else
%     plot(ENW_INS(:,1),YYW_INS.*max(YYWC(:,1)./sum(sum(YYW,2),1)),'ok','linewidth',2);
    for f=2:frag
        YYWCX=[(YYWC(:,f-1)'),fliplr((YYWC(:,f)'))]./sum(sum(YYW,2),1);
        fw(f)=fill(ENWX,YYWCX,'b');
        set(fw(f),'linestyle','-','linewidth',1,'edgecolor','k','facecolor',cc(f-1,:));
    end
    fw(1)=fill([ENW',fliplr(ENW')],[(YYWC(:,frag))',zeros(1,length(YYWC(:,frag)))]./sum(sum(YYW,2),1),'b');
    set(fw(1),'linestyle','-','linewidth',1,'edgecolor','k','facecolor',cc(frag,:));
    plot(ENW,YYWC(:,1)./sum(sum(YYW,2),1),'-k','linewidth',2);
%     plot(ENW_INS(:,1),YYW_INS.*max(YYWC(:,1)./sum(sum(YYW,2),1)),'ok','linewidth',2);
end

hold off
xlim([40 210]); ylim([0 max(max(YYWC))./sum(sum(YYW,2),1).*1.1]);
xlabel('\omega (meV)'); ylabel('P(\omega) (arb. units)'); 
set(gca,'fontname','times new roman','fontsize',18);
set(gca,'XTick',0:50:200);
set(gca,'YTick',[]);
set(gcf,'color','w'); box on;
%
lgd=cell(frag,1);
% lgd{1,:}=sprintf(' INS');
lgd{1,:}=sprintf(' VS');
lgd{2,:}=sprintf(' S');
lgd{3,:}=sprintf(' M');
lgd{4,:}=sprintf(' W');
lgd{5,:}=sprintf(' VW');
lk=legend(lgd); 
set(lk,'box','off','position',[0.1186 0.4374 0.2232 0.4500]);
fig_flag=fig_flag+1; clear f_k;

%
f_k=figure(fig_flag); % stretch interval
hold on
ENSX=[ENS',fliplr(ENS')];
if cd_reverse(2)==0
%     plot(ENS_INS,YYS_INS.*max(YYSC(:,frag)./sum(sum(YYS,2),1)),'ok','linewidth',2);
    fs(1)=area(ENS,YYSC(:,1)./sum(sum(YYS,2),1));
    set(fs(1),'linestyle','-','linewidth',1,'edgecolor','k','facecolor',cc(1,:));
    for f=2:frag
        YYSCX=[(YYSC(:,f-1)'),fliplr((YYSC(:,f)'))]./sum(sum(YYS,2),1);
        fs(f)=fill(ENSX,YYSCX,'b');
        set(fs(f),'linestyle','-','linewidth',1,'edgecolor','k','facecolor',cc(f,:));
    end
%     plot(ENS_INS,YYS_INS.*max(YYSC(:,frag)./sum(sum(YYS,2),1)),'ok','linewidth',2);
    plot(ENS,YYSC(:,frag)./sum(sum(YYS,2),1),'-k','linewidth',2);
else
    for f=2:frag
        YYSCX=[(YYSC(:,f-1)'),fliplr((YYSC(:,f)'))]./sum(sum(YYS,2),1);
        fs(f)=fill(ENSX,YYSCX,'b');
        set(fs(f),'linestyle','-','linewidth',1,'edgecolor','k','facecolor',cc(f-1,:));
    end
    fs(1)=fill([ENS',fliplr(ENS')],[(YYSC(:,frag))',zeros(1,length(YYSC(:,frag)))]./sum(sum(YYS,2),1),'b');
    set(fs(1),'linestyle','-','linewidth',1,'edgecolor','k','facecolor',cc(frag,:));
end
hold off
xlim([270 520]); ylim([0 max(max(YYSC))./sum(sum(YYS,2),1).*1.1]);
xlabel('\omega (meV)'); ylabel('P(\omega) (arb. units)'); 
set(gca,'fontname','times new roman','fontsize',18);
% set(gca,'XTick',250:50:550);
set(gca,'XTick',300:50:500);
set(gca,'YTick',[]);
set(gcf,'color','w'); box on;
%
lgd=cell(frag,1);
% lgd{1,:}=sprintf(' INS');
lgd{1,:}=sprintf(' VS');
lgd{2,:}=sprintf(' S');
lgd{3,:}=sprintf(' M');
lgd{4,:}=sprintf(' W');
lgd{5,:}=sprintf(' VW');
lj=legend(lgd); 
set(lj,'box','off','position',[0.1186 0.47 0.2232 0.4500]);
%
fig_flag=fig_flag+1; clear f_k;

%% Relation between mode frequency and HB length

[~,frag]=size(job_PSDf.PSDf_ps);
EN=job_PSDf.PSDf_en;
ENRRstep=job_PSDf.PSDf_lim(2)-job_PSDf.PSDf_lim(1);
ENRR=[job_PSDf.PSDf_lim-ENRRstep./2,job_PSDf.PSDf_lim(end)+ENRRstep./2]';
ENRRe=(min(ENRR):0.01:max(ENRR))';

% Find maxima
g=fittype('a1./pi./2.*c1./((x-b1).^2+(c1/2).^2)+a2./pi./2.*c2./((x-b2).^2+(c2/2).^2)',...
    'independent','x','dependent','y');
opts=fitoptions('Method','NonlinearLeastSquares');
opts.Display='Off';
opts.StartPoint=[600,600,100,400,50,50];
YY=zeros(length(job_PSDf.PSDf_ps(:,1)),frag); RR=zeros(frag,4); 
for f=1:frag
    YY(:,f)=job_PSDf.PSDf_ps(:,f);
    [fitR,~]=fit(EN,YY(:,f),g,opts);
    ci=confint(fitR); cw=(ci(2,3)-ci(1,3))./2; cs=(ci(2,4)-ci(1,4))./2;
%     figure; plot(fitR,EN,YY(:,f)); xlabel EN; ylabel YY; xlim([0 600]);
    RR(f,:)=[fitR.b1,fitR.b2,cw,cs];
end
clear f;
Pwag=polyfit(ENRR,RR(:,1),2);
Pstr=polyfit(ENRR,RR(:,2),2);

% Plot relation
f_k=figure(fig_flag);
% set(f_k,'Position',[80,15,580,480]);
%
hold on
yyaxis left
set(gca,'ycolor','k');
errorbar(ENRR,RR(:,1),RR(:,3),'ok','linewidth',2);
plot(ENRRe,polyval(Pwag,ENRRe),'-k');
ylim([100 150]); set(gca,'YTick',100:10:160);
ylabel('\delta (O-H) (meV)');
%
yyaxis right
set(gca,'ycolor','r');
errorbar(ENRR,RR(:,2),RR(:,4),'or','linewidth',2);
plot(ENRRe,polyval(Pstr,ENRRe),'-r');
ylim([325 475]); set(gca,'YTick',250:25:500);
ylabel('\nu (O-H) (meV)');
hold off
%
xlabel('Hydrogen bond length (Å)');
xlim([1.6 2.3]); set(gca,'XTick',1.7:0.1:2.2); set(gca,'XTicklabel',{'1.7','1.8','1.9','2.0','2.1','2.2'});
set(gca,'fontname','times new roman','fontsize',18);
set(gcf,'color','w'); box on;
fig_flag=fig_flag+1; clear f_k;

%% Model for legends

% lgd=cell(frag,1);
% lgd{1,:}=sprintf('H··O < x_0-1.5%s Å','\sigma');
% lgd{2,:}=sprintf('H··O [x_0-1.5%s:x_0-0.5%s] Å','\sigma','\sigma');
% lgd{3,:}=sprintf('H··O [x_0-0.5%s:x_0+0.5%s] Å','\sigma','\sigma');
% lgd{4,:}=sprintf('H··O [x_0+0.5%s:x_0+1.5%s] Å','\sigma','\sigma');
% lgd{5,:}=sprintf('H··O > x_0+1.5%s Å','\sigma');
% t=legend(lgd); set(t,'box','off','location','northwest');

% lgd=cell(frag,1);
% lgd{1,:}=sprintf('VS');
% lgd{2,:}=sprintf('S');
% lgd{3,:}=sprintf('M');
% lgd{4,:}=sprintf('W');
% lgd{5,:}=sprintf('VW');
% t=legend(lgd); set(t,'box','off','location','northeast');

% lgd=cell(frag,1);
% lgd{1,:}=sprintf('VS (1%s) : d_{H··O} < 1.602 Å','%');
% lgd{2,:}=sprintf('S (23%s) : d_{H··O} = [1.602; 1.827] Å','%');
% lgd{3,:}=sprintf('M (52%s) : d_{H··O} = [1.827; 2.052] Å','%');
% lgd{4,:}=sprintf('W (23%s) : d_{H··O} = [2.052; 2.277] Å','%');
% lgd{5,:}=sprintf('VW (1%s) : d_{H··O} > 2.277 Å','%');
% t=legend(lgd); set(t,'box','off','location','northwest');
% set(gca,'YTickLabel',[]);
% set(h(1),'color',[204/255 0 0]');
% set(h(2),'color',[1 128/255 0]');
% set(h(3),'color',[204/255 204/255 0]');
% set(h(4),'color',[0 204/255 0]');
% set(h(5),'color',[0 0 1]');
% set(gca,'XMinorTick','on');

% norm=sum(sum(job_PSDf.PSDf_ps,2))/100;
% lgd=cell(frag,1);
% % lgd{1,:}=sprintf('Hydrogen bond < %.2f %s',PSDf_lim(1),char(197));
% lgd{1,:}=sprintf('~%.0f%% of hydrogen bond < %.2f %s',sum(job_PSDf.PSDf_ps(:,1),1)/norm,PSDf_lim(1),char(197));
% for f=2:frag-1
%     lgd{f,:}=sprintf('~%.0f%% of hydrogen bond > %.2f and < %.2f %s',sum(job_PSDf.PSDf_ps(:,f),1)/norm,PSDf_lim(f-1),PSDf_lim(f),char(197));
% end; clear f;
% lgd{frag,:}=sprintf('~%.0f%% of hydrogen bond > %.2f %s',sum(job_PSDf.PSDf_ps(:,frag),1)/norm,PSDf_lim(frag-1),char(197));
% lgd{1,:}=sprintf('H··O < x_0-1.5%s Å','\sigma');
% lgd{2,:}=sprintf('H··O [x_0-1.5%s:x_0-0.5%s] Å','\sigma','\sigma');
% lgd{3,:}=sprintf('H··O [x_0-0.5%s:x_0+0.5%s] Å','\sigma','\sigma');
% lgd{4,:}=sprintf('H··O [x_0+0.5%s:x_0+1.5%s] Å','\sigma','\sigma');
% lgd{5,:}=sprintf('H··O > x_0+1.5%s Å','\sigma');
% legend(lgd);

end

