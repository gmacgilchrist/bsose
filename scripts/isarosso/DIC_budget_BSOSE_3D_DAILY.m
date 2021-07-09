
% this is where rdmds and other functions are
addpath /home/irosso/scripts

clear all;  
%run('DIC_budget_BSOSE_3D')
% ----------------------------------------------------------------------- %
% BUDGET FOR THIS AREA 
xel      = 1001:1080; % last is 1080  
yel      = 1:300;     
zel      = 1:52; 

% time interval
tsnap     = 1;   
tmax      = 73*5; % only 2008 
dt        = 86400;
nt        = tmax-tsnap+1; 
TimeStep  = 24:24:43848;

% specify this:
linear_free_surface_fix = 1;

% folder to save the files in
dir_path = '/data/irosso/data/BSOSE/DIC/3D/DAILY/200pts_6/'

% ----------------------------------------------------------------------- % 
% files to read 

path_1day    = '/data/soccom/SO3/optim2012/ITERATION105/RunWithNewCode4BDGTS/ITERATION105_DICbdgt1dy/OUTPUT_DICbdgt1dy/';


filename1 = 'diag_bgc';        % 'TRAC01  ' 'TRAC02  ' 'TRAC03  ' 'TRAC04  ' 'TRAC05  ' 'TRAC06  ' 'BLGPH3D ' 'BLGOMAR '

filename2 = 'diag_state';      % 'THETA   ' 'SALT    ' 'UVEL    ' 'VVEL    ' 'WVEL    ' 'PHIHYD  ' 'DRHODR  '

filename3 = 'diag_surf';       % 'ETAN    ' 'BLGPCO2 ' 'SIarea  ' 'SIheff  ' 'PHIBOT  '

filename4 = 'diag_dic_budget'; % 'ADVxTr01' 'ADVyTr01' 'ADVrTr01' 'DFxETr01' 'DFyETr01' 'DFrITr01' 'UTRAC01 ' 'VTRAC01 ' 'WTRAC01 ' 'ForcTr01'

filename5 = 'diag_dic_snaps';  % 'TRAC01  ' 

filename6 = 'diag_airsea';     % 'TFLUX   ' 'SFLUX   ' 'BLGCFLX ' 'BLGOFLX ' 'Add2EmP'

filename7 = 'diag_dic_BIOC';   % 'BLGBIOC' % before it was 'BLGBIOA'

filenameU = 'UVEL';
filenameV = 'VVEL';
filenamwW = 'WVEL';


% ----------------------------------------------------------------------- % 
% MITgcm grid
load /data/soccom/GRID_3/grid.mat
nx        = length(xel);
ny        = length(yel);
nz        = length(zel);

len_x     = length(XC(:,1));
len_y     = length(YC(1,:));
len_z     = length(RC);

% cell volume, face areas (for flux calculations)
volume    = zeros(nx,ny,nz);
AREAWEST  = zeros(nx+1,ny,nz);
AREASOUTH = zeros(nx,ny+1,nz);
AREACELL  = zeros(nx,ny,nz+1);
for k=1:nz
    volume(:,:,k)    = hFacC(xel,yel,k).*RAC(xel,yel).*DRF(k);
    AREACELL(:,:,k)  = RAC(xel,yel);
    %AREAWEST(:,:,k)  = DYG([xel xel(end)+1],yel).*DRF(k).*hFacW([xel xel(end)+1],yel,k);
    AREASOUTH(:,:,k) = DXG(xel,[yel yel(end)+1]).*DRF(k).*hFacS(xel,[yel yel(end)+1],k);
end
%AREACELL(:,:,nz+1)   = RAC(xel,yel);

DRF       = DRF(1:nz);
DRF1      = permute(repmat(DRF,[1,nx,ny]),[2,3,1]);

if 0
% ----------------------------------------------------------------------- %
% COMPUTE THE TERMS IN DIC EQUATION AND INTEGRATE AND AVERAGE OVER DEPTH AFTER EACH CALCULATION
% 1.Surface fluxes: gas exchange and precipitation

SURF     = zeros(nx,ny,nz,nt,'single');
PRECIP   = zeros(nx,ny,nz,nt,'single');
SURF_CORR= zeros(nx,ny,nz,nt,'single');

tt       = 1;
disp(' ')
disp('..air-sea DIC fluxes..')
for t = tsnap:tmax
    display(sprintf('t = %d',t)) 

    
    % Air-sea flux
    % BLGCFLX is in mol/m2/s
    % DICTFLX is SURC in the model (increases DIC), in mol/m3/s
    field             = rdmds([path_1day filename6],TimeStep(t),'rec',3);
    %field(field==0)=NaN;
    SURF(:,:,1,tt)    = field(xel,yel)./(DRF(1).*hFacC(xel,yel,1));

    % Precipitation!
    field             = rdmds([path_1day filename4],TimeStep(t),'rec',10);   % EmP = ForcTr01 in mol/m3/s
    %field(field==0)=NaN;
    PRECIP(:,:,1,tt)  = field(xel,yel);
    

    field             = rdmds([path_1day filename4],TimeStep(t),'rec',9);
    %field(field==0)=NaN;
    WTRAC             = field(xel,yel,1);
    SURF_CORR(:,:,1,tt)= -WTRAC./ (DRF(1) * hFacC(xel,yel,1));
    
    tt               = tt + 1;
end


% save the term
file_name =  [dir_path 'DIC_flx_BL2_3D.mat']
  save([file_name],'SURF','PRECIP','-v7.3');

file_name =  [dir_path 'DIC_SURF_CORR_BL2_3D.mat']
 save([file_name],'SURF_CORR','-v7.3'); 

clear PRECIP SURF
end

if 0
% ----------------------------------------------------------------------- %
% 2. calculate DIC concentration from 1day snapshots

DIC_snap  = zeros(nx,ny,nz,nt,'single');

tt        = 1;
disp(' ')
disp('..reading DIC...')
for t = tsnap:tmax+1
    display(sprintf('t = %d',t))
    
    % read DIC 1day snapshots
    if t==1
        field = rdmds([path_1day filename5],0,'rec',1);
    else
        field = rdmds([path_1day filename5],TimeStep(t-1),'rec',1);
    end
    %field(field==0)=NaN;
    DIC_snap(:,:,:,tt) = field(xel,yel,zel);
    
    tt = tt + 1;
end

% tendency of DIC 
% (change in DIC content / storage)
dCdt        = (DIC_snap(:,:,:,2:end) - DIC_snap(:,:,:,1:end-1))/dt;
dCdt_corr   = dCdt - SURF_CORR(:,:,:,1:nt); %-1);
clear DIC_snap field

% save the file 
file_name = [dir_path 'dCdt_BL2_3D.mat']
save([file_name],'dCdt', 'dCdt_corr','-v7.3');

clear dCdt dCdt_corr
end

if 0
% ----------------------------------------------------------------------- %
% 3. Advection (div(bar(u)*dic))

ADV_x     = zeros(nx,ny,nz,nt,'single');
ADV_y     = zeros(nx,ny,nz,nt,'single');
ADV_z     = zeros(nx,ny,nz,nt,'single');

tt        = 1;
disp(' ')
disp('..advection..')
for t = tsnap:tmax 
    display(sprintf('t = %d',t)) 
    % read advective flux
    % for calculation of derivatives: pad with 0
    field = zeros(len_x+1, len_y+1, len_z+1, 3); 
    field(1:end-1,1:end-1,1:end-1,:) = rdmds([path_1day filename4],TimeStep(t),'rec',1:3);
    %field(field==0)=NaN;
    FLUXx = field([xel xel(end)+1],yel,zel,1);
    FLUXy = field(xel,[yel yel(end)+1],zel,2);
    FLUXz = field(xel,yel,[zel zel(end)+1],3);
            
    % derivatives 
    ADV_x(:,:,:,tt) = (FLUXx(2:end,:,:)-FLUXx(1:end-1,:,:))./volume;
    ADV_y(:,:,:,tt) = (FLUXy(:,2:end,:)-FLUXy(:,1:end-1,:))./volume;
    ADV_z(:,:,:,tt) = -(FLUXz(:,:,2:end)-FLUXz(:,:,1:end-1))./volume;

    tt = tt + 1;

end
clear field FLUXz FLUXy FLUXx

% terms for the equation
ADVn        = -(ADV_x+ADV_y+ADV_z);
ADV_corr    = ADVn + SURF_CORR;
clear SURF_CORR

ADVn(isnan(ADVn))=0.;
ADV_corr(isnan(ADV_corr))=0.;
ADV_x(isnan(ADV_x))=0.;
ADV_y(isnan(ADV_y))=0.;
ADV_z(isnan(ADV_z))=0.;

% save the terms
file_name =  [dir_path 'ADV_DIC_BL2_3D.mat']
save([file_name],'ADVn','ADV_corr','ADV_x','ADV_y','ADV_z','-v7.3');

clear ADV_corr ADVn ADV_x ADV_y ADV_z
end



if 0
% ----------------------------------------------------------------------- %
% 4. Divergence

% CORRECTED tendency of DIC in z_star
% (change in DIC content / storage)
DIV_t     = zeros(nx,ny,nz,nt,'single');
DIV_x     = zeros(nx,ny,nz,nt,'single');
DIV_y     = zeros(nx,ny,nz,nt,'single');
DIV_z     = zeros(nx,ny,nz,nt,'single');

tt        = 1;
disp(' ')
disp('..divergence..') 
for t = tsnap:tmax
    display(sprintf('t = %d',t))
    % calculate divergence of the flow
    % read state variables
    field = zeros(len_x+1, len_y+1, len_z+1,3);
    field(1:end-1,1:end-1,1:end-1,:) = rdmds([path_1day filename2],TimeStep(t),'rec',3:5);
    %field(field==0)=NaN; 

    U     = field([xel xel(end)+1],yel,zel,1).*AREAWEST;
    V     = field(xel,[yel yel(end)+1],zel,2).*AREASOUTH;
    W     = field(xel,yel,[zel zel(end)+1],3).*AREACELL;
    
    % linear free surface: tracer advection does not "feel" the free
    % surface height changes, hence horizontal divergence leads to an
    % apparent change in volume, which results in a change in tracer
    % concentration (though this is not a real effect, and we should not 
    % use the linear free surface for budget calculations).
    if linear_free_surface_fix == 1
        % Setting w=0 at the surface gives the change in volume
        W(:,:,1) = 0;
    end

    % read DIC 1day averages
    field  = rdmds([path_1day filename1],TimeStep(t),'rec',1);
    DIC    = field(xel,yel,zel);

    % gradients:
    DIV_x(:,:,:,tt)  = ((U(2:end,:,:)-U(1:end-1,:,:))./volume).*DIC;
    DIV_y(:,:,:,tt)  = ((V(:,2:end,:)-V(:,1:end-1,:))./volume).*DIC;
    DIV_z(:,:,:,tt)  = ((W(:,:,1:end-1)-W(:,:,2:end))./volume).*DIC;
      
    % compute DIC*DIV_vec(U)
    DIV_t(:,:,:,tt)  = (DIV_x(:,:,:,tt)+DIV_y(:,:,:,tt)+DIV_z(:,:,:,tt));

    tt = tt + 1;
end

DIV_t(isnan(DIV_t))=0.;
DIV_x(isnan(DIV_x))=0.;
DIV_y(isnan(DIV_y))=0.;
DIV_z(isnan(DIV_z))=0.;

% save the terms
file_name =  [dir_path 'DIV_DIC_BL2_3D.mat']
save([file_name],'DIV_t','DIV_x','DIV_y','DIV_z','-v7.3');
clear DIV_t DIV_x DIV_y DIV_z
end

if 0
% ----------------------------------------------------------------------- %
% 5. Diffusion

DIFF_x     = zeros(nx,ny,nz,nt,'single');
DIFF_y     = zeros(nx,ny,nz,nt,'single');
DIFF_z     = zeros(nx,ny,nz,nt,'single');

tt         = 1;
disp(' ')
disp('..diffusion..')
for t = tsnap:tmax     
    display(sprintf('t = %d',t)) 
    % read diffusive flux
    % for calculation of derivatives: pad with 0
    field   = zeros(len_x+1,len_y+1,len_z+1,3);
    field(1:end-1,1:end-1,1:end-1,:) = rdmds([path_1day filename4],TimeStep(t),'rec',4:6);
    %field(field==0)=NaN;
    FLUXx   = field([xel xel(end)+1],yel,zel,1);
    FLUXy   = field(xel,[yel yel(end)+1],zel,2);
    % There is no explicit vertical diffusive flux, only implicit
    FLUXz   = field(xel,yel,[zel zel(end)+1],3);
    % derivatives
    DIFF_x(:,:,:,tt)  = (FLUXx(2:end,:,:)-FLUXx(1:end-1,:,:))./volume;
    DIFF_y(:,:,:,tt)  = (FLUXy(:,2:end,:)-FLUXy(:,1:end-1,:))./volume;
    DIFF_z(:,:,:,tt)  = -(FLUXz(:,:,2:end)-FLUXz(:,:,1:end-1))./volume;

    tt = tt + 1;
end
clear field FLUXx FLUXy FLUXz

% terms for the equation
DIFF_horiz  = DIFF_x + DIFF_y;
DIFF        = DIFF_horiz + DIFF_z; 
DIFFn       = -DIFF;
clear DIFF_x DIFF_y


% save the terms
file_name  =  [dir_path 'DIFF_DIC_BL2_3D.mat']
  save([file_name],'DIFFn', 'DIFF_horiz','DIFF_z','-v7.3')

%clear DIFFn DIFF_horiz DIFF_z

end

if 0
% ----------------------------------------------------------------------- %
% 6. Biological uptake and remineralization

BIO     = zeros(nx,ny,nz,nt,'single');

tt      = 1;
disp(' ')
disp('..biology..')
for t = tsnap:tmax  
    display(sprintf('t = %d',t)) 
    % Biology
    field = rdmds([path_1day filename7],TimeStep(t),'rec',1);
    %field(field==0)=NaN;
    % BLGBIOA is the change in DIC due to all biological activity (production, remin and calcium carb cycling)
    BIO(:,:,:,tt) = field(xel,yel,zel);
    
    tt = tt + 1;
end
clear field

% save the terms
file_name  =  [dir_path 'BIO_BL2_3D.mat']
save([file_name],'BIO','-v7.3')
%clear BIO

end
% ----------------------------------------------------------------------- %



if 1
% ----------------------------------------------------------------------- %
% PLOT

% plot to check the budget is closed..

disp('plotting')

load( [dir_path 'DIC_flx_BL2_3D.mat'])
load( [dir_path 'dCdt_BL2_3D.mat'])
load( [dir_path 'ADV_DIC_BL2_3D.mat'])
load( [dir_path 'DIFF_DIC_BL2_3D.mat'])
load( [dir_path 'BIO_BL2_3D.mat'])


t = 1;
xx = 5;
yy = 120;

RES = dCdt(:,:,:,t) - ADV_corr(:,:,:,t) - DIFFn(:,:,:,t) - BIO(:,:,:,t) - SURF(:,:,:,t) - PRECIP(:,:,:,t);

figure(1);
plot(squeeze(dCdt(xx,yy,1:33,t)),RC(1:33),'k-')
hold on
plot(squeeze(-ADV_corr(xx,yy,1:33,t)),RC(1:33),'b-')
plot(squeeze(-DIFFn(xx,yy,1:33,t)),RC(1:33),'g-')
plot(squeeze(-SURF(xx,yy,1,t)),RC(1),'m*')
plot(squeeze(-PRECIP(xx,yy,1:33,t)),RC(1:33),'c-')
plot(squeeze(-BIO(xx,yy,1:33,t)),RC(1:33),'r-')
plot(squeeze(RES(xx,yy,1:33)),RC(1:33),'k--')
legend('dCdt', '-ADV_{corr}', '-DIFF', '-SURF','-PRECIP','-BIO','RES')


figure(2);
subplot(241)
imagesc(dCdt(:,:,1).'); colorbar;
title('dCdt');
subplot(242)
imagesc(ADV_corr(:,:,1,t).'); colorbar;
title('ADV');
subplot(243)
imagesc(BIO(:,:,1,t).'); colorbar;
title('BIO');
subplot(244)
imagesc(DIFF(:,:,1,t).'); colorbar;
title('DIFF');
subplot(245)
imagesc(SURF(:,:,1,t).'); colorbar;
title('SURF');
subplot(246)
imagesc(PRECIP(:,:,1,t).'); colorbar;
title('PRECIP');
subplot(247)
imagesc(DIV_t(:,:,1,t).'); colorbar;
title('DIV');
subplot(248)
imagesc(RES(:,:,1,t).'); colorbar;
title('RES');


end
