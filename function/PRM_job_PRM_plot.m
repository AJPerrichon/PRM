function [fig_flag] = PRM_job_PRM_plot(command,structure,job_PRM,fig_flag)

% Plot job 1 PRM

meshnorm=length(0:command.PRM_resolution:max(structure.cell_parameters));

AAAA=job_PRM.bs_log;
% AAAA=job_PRM.bs_plan;

f_e=figure(fig_flag);
if command.PRM_direction==1; surf(job_PRM.yvec,job_PRM.zvec,AAAA,'lines','none'); end
if command.PRM_direction==2; surf(job_PRM.xvec,job_PRM.zvec,AAAA,'lines','none'); end
if command.PRM_direction==3; surf(job_PRM.xvec,job_PRM.yvec,AAAA,'lines','none'); end
colormap jet; set(gcf,'color','w'); axis image; box on; grid off;
xmin=0.5*(1-1/command.supercell(2)); xmax=0.5*(1+1/command.supercell(2));
ymin=0.5*(1-1/command.supercell(1)); ymax=0.5*(1+1/command.supercell(1));
xlim([xmin xmax]); ylim([ymin ymax]);
set(gca,'Xtick',[xmin ((xmax-xmin)*0.25+xmin)...
    0.5 ((xmax-xmin)*0.75+xmin) xmax]);
set(gca,'Ytick',[ymin ((ymax-ymin)*0.25+ymin)...
    0.5 ((ymax-ymin)*0.75+ymin) ymax]);
set(gca,'XTickLabel',{sprintf('%u',round((0.5-(abs(0.5-xmin)*command.supercell(2)))*1000)/1000),...
    sprintf('%0.2f',0.5-(abs(0.5-((xmax-xmin)*0.25+xmin))*command.supercell(2))),...
    '0.5',sprintf('%0.2f',0.5+(abs(0.5-((xmax-xmin)*0.75+xmin))*command.supercell(2))),...
    sprintf('%1.0f',0.5+(abs(0.5-xmax)*command.supercell(2)))});
set(gca,'YTickLabel',{sprintf('%u',round((0.5-(abs(0.5-ymin)*command.supercell(1)))*1000)/1000),...
    sprintf('%0.2f',0.5-(abs(0.5-((ymax-ymin)*0.25+ymin))*command.supercell(1))),...
    '0.5',sprintf('%0.2f',0.5+(abs(0.5-((ymax-ymin)*0.75+ymin))*command.supercell(1))),...
    sprintf('%1.0f',0.5+(abs(0.5-ymax)*command.supercell(1)))});
% title(sprintf('Positional Recurrence Maps in %s space', command.space));
if command.PRM_direction==1; xdir='y'; ydir='z'; end
if command.PRM_direction==2; xdir='x'; ydir='z'; end
if command.PRM_direction==3; xdir='x'; ydir='y'; end
if command.PRM_direction==1 && command.supercell(1)~=command.supercell(2); daspect([1,2,1]); end % to fix
if command.PRM_direction==2 && command.supercell(1)~=command.supercell(2); daspect([1,2,1]); end
% xlabel(sprintf('%s with respect to 1x1x1 supercell (r.l.u.)',xdir));
xlabel(sprintf('%s (r.l.u.)',xdir));
% ylabel(sprintf('%s with respect to 1x1x1 supercell (r.l.u.)',ydir));
ylabel(sprintf('%s (r.l.u.)',ydir));
% caxis([log(0.005*meshnorm/structure.total_count) log(1*meshnorm/structure.total_count)]); % arbitrary, [-13.5 -5.8]
% cb_scale=[log(0.1/structure.total_count),log(1/structure.total_count),log(10/structure.total_count),log(100/structure.total_count),log(1000/structure.total_count),log(10000/structure.total_count)];
% c=colorbar('YTick',cb_scale);
% ax=gca; axpos=ax.Position; cpos =c.Position; cpos(3)=0.5*cpos(3);
% c.Position=cpos; ax.Position=axpos; clear axpos cpos ax;
% set(c,'YTickLabel',{sprintf('%1.0f',0.1),sprintf('%1.0f',1),sprintf('%1.0f',10),...
%     sprintf('%1.0f',100),sprintf('%1.0f',1000),sprintf('%1.0f',10000)});
% ylabel(c,sprintf('Count per d^3r = %1.2e angstrom^3',(power(command.PRM_resolution,3))));
set(gca,'fontname','times new roman','fontsize',22);
set(f_e,'Position',[80,15,880,480]);
if rem(round(command.supercell(1)*1000)/1000,1)==0 || rem(round(command.supercell(2)*100)/100,1)==0;...
        view(0,90); else view(45,90); end
fig_flag=fig_flag+1;
clear xmin xmax ymin ymax cb_scale xdir ydir cb_pos c f_e;

end

