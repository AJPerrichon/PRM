function [fig_flag] = PRM_job_Kin_plot(job_Kin,fig_flag)

% Plot job Kin

f_l=figure(fig_flag);
set(f_l,'Position',[80,15,480,480]);
plot(job_Kin.ki_t,smooth(job_Kin.ki_T,9),'-b');
ki_Teff=mean(job_Kin.ki_T);
axis([0 max(job_Kin.ki_t) 0 max(job_Kin.ki_T)]);
title(sprintf('Effective temperature; <T>=%3.3f K',ki_Teff));
xlabel('Time (fs)');
ylabel('Temperature (K)');
fig_flag=fig_flag+1;
set(gcf,'color','w');
clear f_l;

end

