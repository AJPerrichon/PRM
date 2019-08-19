function [fig_flag] = PRM_job_PRM_plot3D(command,structure,job_PRM,fig_flag)

% Plot job 1 PRM in 3D

if command.sym==1 && command.plot(1,2)==1
    
    % Test 1
    bs_iso=job_PRM.bs_sym./max(max(max(job_PRM.bs_sym)));
    f_i=figure(fig_flag);
    set(f_i,'Position',[80,15,801,801]);
    [ps_xvec,ps_yvec,ps_zvec]=meshgrid(job_PRM.xvec,job_PRM.yvec,job_PRM.zvec);
    is=patch(isosurface(ps_xvec,ps_yvec,ps_zvec,bs_iso,0.05));
    axis vis3d;
    set(is,'FaceColor','k','EdgeColor','r');
    camlight left;
    camlight;
    lighting gouraud;
    fig_flag=fig_flag+1;
    set(gcf,'color','w');
    clear f_i is;

    % Test 2
    data=log(job_PRM.bs_sym.*structure.total_count.*16); data(data<0)=0; isoval=.2;
    % data=job_PRM.bs_sym.*structure.total_count.*16; isoval=.01;
    data=data./max(max(max(data)));
    xvec=job_PRM.xvec; yvec=job_PRM.yvec; zvec=job_PRM.zvec;
    % data(151:301,:,:)=[]; xvec(151:301)=[];
    % data(:,151:301,:)=[]; yvec(151:301)=[];
    % data(:,:,151:301)=[]; zvec(151:301)=[];
    f_j=figure(fig_flag); set(f_j,'Position',[80,15,801,801]);
    [ps_xvec,ps_yvec,ps_zvec]=meshgrid(yvec,xvec,zvec);
    h=patch(isosurface(ps_xvec,ps_yvec,ps_zvec,data,isoval),...
        'FaceColor',[0 0 0.8],...
        'EdgeColor','none');
    isonormals(data,h);
    patch(isocaps(ps_xvec,ps_yvec,ps_zvec,data,isoval),...
        'FaceColor','interp',...
        'EdgeColor','none');
    colormap(jet(100))
    daspect([1,1,1])
    axis tight
    view(0,90)
    camlight
    fig_flag=fig_flag+1;
    set(gcf,'color','w');
    
    % Test 3
    C=job_PRM.bs_sym.*100./max(max(max(job_PRM.bs_sym)));
    x=job_PRM.xvec; y=job_PRM.yvec; z=job_PRM.zvec;
    [X,Y,Z]=meshgrid(y,x,z);
    Cx=zeros(size(C));
    for n=1:length(z); Cx(:,:,n)=medfilt2(C(:,:,n),[10 10]); end
    Cc=zeros(size(permute(C,[3 2 1])));
    for n=1:length(x); Cc(:,:,n)=medfilt2(permute(Cx(n,:,:),[3 2 1]),[10 1]); end
    Cc=permute(Cc,[3 2 1]);
    S1=isosurface(X,Y,Z,Cc,0.05);
    P1=patch(S1);
    % isonormals(X,Y,Z,C,P1);
    set(P1,'FaceColor','yellow','EdgeColor','none','FaceAlpha',0.1);
    daspect([1,1,1])
    view(3); axis tight;
    camlight; lighting gouraud;
    % zlim([0.375 0.625])
    zlim([0.65 0.85])

else
end

end

