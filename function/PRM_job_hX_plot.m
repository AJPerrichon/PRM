function [fig_flag] = PRM_job_hX_plot(command,job_hX,fig_flag)

% Plot job h(MSD)
f_g=figure(fig_flag);
set(f_g,'Position',[80,15,480,480]);
plot(job_hX.od_vec,job_hX.od_void,'-ob');
title(sprintf('Histogram of MSD in %s space', command.space));
xlabel('Distance (Å)');
ylabel(sprintf('Count per dr=%0.3f Å',command.PRM_resolution));
set(gca,'fontname','times new roman','fontsize',22);
set(gcf,'color','w'); box on;
fig_flag=fig_flag+1;
clear f_g;

% Plot job h(POS)
f_h=figure(fig_flag);
set(f_h,'Position',[80,15,480,480]);
hold on
if command.sym==1
    if strcmp(command.sym_op,'4/mmm') || strcmp(command.sym_op,'4/m')
        plot(job_hX.id_vec_x(:),job_hX.id_x(:),'-ob','linewidth',2);
        plot(job_hX.id_vec_z(:),job_hX.id_z(:),'-or','linewidth',2);
        if rem(round(command.supercell(1)*1000)/1000,1)==0 || rem(round(command.supercell(2)*100)/100,1)==0;...
                legend('x=y','z'); else legend('m','z'); end
    elseif strcmp(command.sym_op,'m-3m') || strcmp(command.sym_op,'m-3')
        plot(job_hX.id_vec_x(:),job_hX.id_x(:),'-ob','linewidth',2);  
        legend('x=y=z');
    elseif  strcmp(command.sym_op,'mmm')
        plot(job_hX.id_vec_x(:),job_hX.id_x(:),'-ob','linewidth',2);
        plot(job_hX.id_vec_y(:),job_hX.id_y(:),'-og','linewidth',2);
        plot(job_hX.id_vec_z(:),job_hX.id_z(:),'-or','linewidth',2);
        legend('x','y','z');
    else
    end
else
    plot(job_hX.id_vec_x,job_hX.id_x(:),'-ob','linewidth',2);
    plot(job_hX.id_vec_y,job_hX.id_y(:),'-og','linewidth',2);
    plot(job_hX.id_vec_z,job_hX.id_z(:),'-or','linewidth',2);
    if rem(round(command.supercell(1)*1000)/1000,1)==0 || rem(round(command.supercell(2)*100)/100,1)==0;...
            legend('x','y','z'); else legend('m_1','m_2','z'); end
end
hold off
title(sprintf('Histogram of positions in %s space',command.space));
xlabel('Distance (Å)');
ylabel(sprintf('Count per dr=%0.3f Å',command.PRM_resolution));
set(gca,'fontname','times new roman','fontsize',22);
xlim([0 1]); ylim([0 0.35]);
set(gcf,'color','w'); box on;
fig_flag=fig_flag+1;
clear f_h;

end

