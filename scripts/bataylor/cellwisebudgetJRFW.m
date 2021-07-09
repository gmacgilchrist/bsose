% B-SOSE DIC budget for JRFW
% Ariane Verdy, May 2018
% edited Ben Taylor, Jan 2 2020

% I'm rewriting this with no pretenses of flexibility. 
% This script will only work for the specific files Matt sent me on Dec 28
% or so.

clear all
close all

tracername = 'DIC';
tracernum = 1;

% file name
%diag_budget_file = 'diag_dic_budget';
budgetfolder = '../../../data/bSOSE/iter129jrfw/monthly/';

% for 2008-2012 monthly output
% whichrun = '2008';
% %dts = 730.2; tmax=43812;
% cd /data/averdy/soccom/bgc_budgets/diags/
% diag_state_file = 'diag_state';
% dt = 2628900; 




% % for rdmds
% addpath ~/scripts_m
% 
% % for plotting
% addpath ~/scripts_m/
% addpath ~/scripts_m/m_map/
% load ~/scripts_m/redbluecolormap.mat

% time stepping
ts = 1 %tmax); % I don't think I use this ... 
nt = length(ts);% or this? 
monstart = 2; % starting month
monend =2; % ending month

%right now I only want one time step, at month 5

% select area
x = 1:2160; % all longitudes
ylen = 570;
y = 1:ylen; % everywhere poleward of pretty far north; almost everywhere
% don't let y= 588


z = 1:52; % top 1000 m


load grid.mat 
% XC = ncread('../../../data/bSOSE/iter122/setup/grid.nc', 'XC');
% YC = ncread('../../../data/bSOSE/iter122/setup/grid.nc', 'YC');
% RC = ncread('../../../data/bSOSE/iter122/setup/grid.nc', 'RC');
% hFacC = ncread('../../../data/bSOSE/iter122/setup/grid.nc', 'hFacC');
% hFacS = ncread('../../../data/bSOSE/iter122/setup/grid.nc', 'hFacS');
% hFacW = ncread('../../../data/bSOSE/iter122/setup/grid.nc', 'hFacW');
% RF = ncread('../../../data/bSOSE/iter122/setup/grid.nc', 'RF');
% DRC = ncread('../../../data/bSOSE/iter122/setup/grid.nc', 'DRC');
% RAC = ncread('../../../data/bSOSE/iter122/setup/grid.nc', 'RAC');
% 
% DRF = ncread('../../../data/bSOSE/iter122/setup/grid.nc', 'DRF');
% DXG = ncread('../../../data/bSOSE/iter122/setup/grid.nc', 'DXG');
% DYG = ncread('../../../data/bSOSE/iter122/setup/grid.nc', 'DYG');
[nx,ny,nz] = size(hFacC);
dz = permute(repmat(DRF(z),[1,nx,ylen]),[2,3,1]).*hFacC(x,y,z); % what is this?

% What is WTRAC01? 
% Why are we dividing it by dz?
% What is dz? 
% What do permute and repmat do? 

%save ('grid.mat' ,'hFacS', 'hFacC', 'hFacW', 'XC', 'YC', 'RC', 'RF', 'DRC','RAC','DRF','DXG','DYG')
%%
% cell volume, face areas (for flux calculations)

volume = zeros(nx,ylen,nz);
areaWest = zeros(nx+1,ylen,nz);
areaSouth = zeros(nx,ylen+1,nz);
areaTop = zeros(nx,ylen,nz+1);
for k=1:nz
 volume(:,y,k) = hFacC(x,y,k).*RAC(x,y)*DRF(k);
 areaTop(:,y,k) = RAC(x,y);
 if x(end)==nx
  areaWest(:,:,k)  = DYG([x 1],y).*DRF(k).*hFacW([x 1],y,k);
 else
  areaWest(:,:,k)  = DYG([x x(end)+1],y).*DRF(k).*hFacW([x x(end)+1],y,k);
 end
 areaSouth(:,:,k) = DXG(x,[y ylen+1]).*DRF(k).*hFacS(x,[y ylen+1],k); 
end
areaTop(:,:,nz+1) = RAC(x,y);
area = RAC(x,y);


% initialize
surf = single(zeros(nx,ylen,nz,nt));
dilut = single(zeros(nx,ylen,nz,nt));
bio = single(zeros(nx,ylen,nz,nt));
mix = single(zeros(nx,ylen,nz,nt));
adv = single(zeros(nx,ylen,nz,nt));
adv_h = single(zeros(nx,ylen,nz,nt));
adv_v = single(zeros(nx,ylen,nz,nt));
corr = single(zeros(nx,ylen,nz,nt));
div = single(zeros(nx,ylen,nz,nt));
tend = single(zeros(nx,ylen,nz,nt));

%%
% read diagnostics
% calculate tendencies in mol/m3/s
for t=1:monend-monstart+1

% tendency due to air-sea flux 
% diagnostic: BLGCFLX
% surface flux (mol/m3/s)
% rdmds reads these MIT GCM files. Why do I need it? as opposed to just
% ncread? I truly don't know... 

% flux = rdmds('diag_airsea',ts(t),'rec',3); 
flux = ncread(strcat(budgetfolder,'bsose_i129_2013to2018_monthly_surfCO2flx.nc'), 'BLGCFLX', [1 1 t+monstart-1], [Inf ylen 1]);
surf(:,:,1,t) = flux(x,y)./(DRF(1)*squeeze(hFacC(x,y,1)));
%%
% tendency due to dilution
% diagnostic: FORCTR01
% forcing tendency (mol/m3/s) includes effects of E-P-R and sponge layer contributions

dilut(:,:,1,t) = ncread(strcat(budgetfolder,'bsose_i129_2013to2018_monthly_ForcTr01.nc'), 'ForcTr01', [1 1 1 t+monstart-1], [Inf ylen 1 1]);

% tendency due to biology
% diagnostic: BLGBIOC
% - uptake + remin + carbonate system (mol/m3/s)
% if strcmp(whichrun,'2013_5d')
%  tmp = rdmds(diag_budget_file,ts(t),'rec',1000);
% elseif strcmp(whichrun,'2013')
%  tmp = rdmds(diag_budget_file,ts(t),'rec',8);
% elseif strcmp(whichrun,'2008')
%  tmp = rdmds(diag_budget_file,ts(t),'rec',8);
% else
%  display('error')
% end
bio(:,:,:,t) = ncread(strcat(budgetfolder,'bsose_i129_2013to2018_monthly_BLGBIOC.nc'), 'BLGBIOC', [1 1 1 t+monstart-1], [Inf ylen Inf 1]);

%%
% advection
% diagnostics: ADVxTr01, ADVyTr01, ADVrTr01
% advective flux (mol/s)
% advflux = rdmds(diag_budget_file,ts(t),'rec',1:3);
advflux_x(x,:,:) = ncread(strcat(budgetfolder,'bsose_i129_2013to2018_monthly_ADVx_Tr01.nc'), 'ADVxTr01', [1 1 1 t+monstart-1], [Inf ylen Inf 1]);
advflux_y(:,:,:) = ncread(strcat(budgetfolder,'bsose_i129_2013to2018_monthly_ADVy_Tr01.nc'), 'ADVyTr01', [1 1 1 t+monstart-1], [Inf ylen+1 Inf 1]);
% advrDIC is still probably not complete - so should delete that dude!
advflux_z(:,:,z) = ncread(strcat(budgetfolder,'bsose_i129_2013to2018_monthly_ADVr_Tr01.nc'), 'ADVrTr01', [1 1 1 t+monstart-1], [Inf ylen Inf 1]);
% why are they using x,y,z ? 
if x(end)==nx
 advflux_x = advflux_x([x 1],y,z);
else
 advflux_x = advflux_x([x x(end)+1],y,z);
end

% if y(end)==ny
%     advflux_y = advflux_y(x,[y y(end)],z);
% else 
% advflux_y = advflux_y(x,[y y(end)+1],z);
% end
if z(end)==nz
 advflux_z = advflux_z(x,y,[z z(end)]); 
 advflux_z(:,:,end) = 0*advflux_z(:,:,end);
else
 advflux_z = advflux_z(x,y,[z z(end)+1]); 
end
% what does diff do to them? Just takes differences between adjacent
% values.
adv_x = diff(advflux_x,1,1)./volume;
adv_y = diff(advflux_y,1,2)./volume;
adv_z = -diff(advflux_z,1,3)./volume;

% minus sign because advective flux is on lhs, moving to rhs
adv(:,:,:,t) = -(adv_x + adv_y + adv_z);

 

% % THIS STUFF IS WRONG. tryna write a divergence term using the UTRAC01 terms.
% utrans = ncread(strcat(budgetfolder,'bsose_i129_2013to2018_monthly_UTRAC01.nc'), 'UTRAC01', [1 1 1 t+monstart-1], [Inf ylen Inf 1]);
% vtrans = ncread(strcat(budgetfolder,'bsose_i129_2013to2018_monthly_VTRAC01.nc'), 'VTRAC01', [1 1 1 t+monstart-1], [Inf ylen+1 Inf 1]);
% wtrans = ncread(strcat(budgetfolder,'bsose_i129_2013to2018_monthly_WTRAC01.nc'), 'WTRAC01', [1 1 1 t+monstart-1], [Inf ylen Inf 1]);
% 
% % calculate div from diff's
% div(1:2159, :, :) = div(1:2159, :, :) + diff(utrans,1);
% div(2160,:, : ) = div(2160, : ,:) + (utrans(2160,:,:) - utrans(1,:,:));
% div(:,:,:) = div(:,:,:) + diff(vtrans,1,2);
% div(:,:,1:51) = div(:,:,1:51) + diff(wtrans,1, 3);
% div(:,:,52) = div(:,:,52) + wtrans(:,:,52);

% divergence, correctly written. 
% diagnostics: UVEL, VVEL, WVEL
%vel = rdmds(diag_state_file,ts(t),'rec',3:5);
% vel = zeros(nx,ylen+1,nz,3);
% 
% vel(:,1:ylen,:,1) = ncread(strcat(budgetfolder,'../bsose_i129_2013to2018_monthly_Uvel.nc'), 'UVEL', [1 1 1 t+monstart-1], [Inf ylen Inf 1]);
% vel(:,:,:,2) = ncread(strcat(budgetfolder,'../bsose_i129_2013to2018_monthly_Vvel.nc'), 'VVEL', [1 1 1 t+monstart-1], [Inf ylen+1 Inf 1]);
% vel(:,1:ylen,:,3) = ncread(strcat(budgetfolder,'../bsose_i129_2013to2018_monthly_Wvel.nc'), 'WVEL', [1 1 1 t+monstart-1], [Inf ylen Inf 1]);
% if x(end)==nx
%  U = vel([x 1],y,z,1).*areaWest;
% else
%  U = vel([x x(end)+1],y,z,1).*areaWest;
% end
% V= vel(x,[y ylen+1],z,2).*areaSouth;
% 
% if z(end)==nz
%  W = vel(x,y,[z z(end)],3).*areaTop; 
%  W(:,:,end) = 0*W(:,:,end);
% else
%  W = vel(x,y,[z z(end)+1],3).*areaTop; 
% end
% div_x=diff(U,1,1)./volume;
% div_y=diff(V,1,2)./volume;
% div_z=-diff(W,1,3)./volume;
% 
% % tracer field (mol/m3)
% tracer = ncread(strcat(budgetfolder,'../bsose_i122_2013to2017_monthly_DIC.nc'), 'TRAC01', [1 1 1 t+monstart-1], [Inf ylen Inf 1]);
% 
% % tmp = rdmds('diag_bgc',ts(t),'rec',1);
% % tracer = tmp(x,y,z);
% div(:,:,:,t) = tracer.*(div_x + div_y + div_z);
% 
% % advection components
% adv_h(:,:,:,t) = -(adv_x + adv_y)+tracer.*(div_x + div_y);
% adv_v(:,:,:,t) = -(adv_z)+tracer.*(div_z);

% correction to vertical advection at z=0
% diagnostic: WTRAC01
tmp = ncread(strcat(budgetfolder,'bsose_i129_2013to2018_monthly_WTRAC01.nc'), 'WTRAC01', [1 1 1 t+monstart-1], [Inf ylen 1 1]);
corr(:,:,1,t) = tmp(x,y,1)./dz(x,y,1);
%%
% mixing
% diagnostics: DFxETr01, DFyETr01, DFrITr01
% diffusive flux (mol/s)
%diffflux = rdmds(diag_budget_file,ts(t),'rec',4:6);
diffflux = zeros(nx, ylen+1, nz, 3);
diffflux(:,1:ylen,:,1) = ncread(strcat(budgetfolder,'bsose_i129_2013to2018_monthly_DFxE_Tr01.nc'), 'DFxETr01', [1 1 1 t+monstart-1], [Inf ylen Inf 1]);
diffflux(:,:,:,2) = ncread(strcat(budgetfolder,'bsose_i129_2013to2018_monthly_DFyE_Tr01.nc'), 'DFyETr01', [1 1 1 t+monstart-1], [Inf ylen+1 Inf 1]);
diffflux(:,1:ylen,:,3) = ncread(strcat(budgetfolder,'bsose_i129_2013to2018_monthly_DFrI_Tr01.nc'), 'DFrITr01', [1 1 1 t+monstart-1], [Inf ylen Inf 1]);
if x(end)==nx
 diffflux_x = diffflux([x 1],y,z,1);
else
 diffflux_x = diffflux([x x(end)+1],y,z,1);
end
diffflux_y = diffflux(x,[y ylen+1],z,2);
if z(end)==nz
 diffflux_z = diffflux(x,y,[z z(end)],3); 
 diffflux_z(:,:,end) = 0*diffflux_z(:,:,end);
else
 diffflux_z = diffflux(x,y,[z z(end)+1],3); 
end
mix_x = diff(diffflux_x,1,1)./volume;
mix_y = diff(diffflux_y,1,2)./volume;
mix_z = -diff(diffflux_z,1,3)./volume;

% minus sign because diffusive flux is on lhs, moving to rhs
mix(:,:,:,t) = -(mix_x + mix_y + mix_z);
%%
% total tendency
dt = 3600*24*30.4306; %approximate number of seconds in a month
snap = ncread(strcat(budgetfolder,'bsose_i129_2013to2018_MonthlySnapshots_TRAC01.nc'), 'TRAC01', [1 1 1 t+monstart-1], [Inf ylen Inf 2]);
tend(:,:,:,t) = diff(snap(x,y,z,:),1,4)/dt;

end % for t

clear tmp snap flux diffflux* div_* vel U V W advflux* %mix_* \%



% fix divisions by hFacC=0
surf(isnan(surf)) = 0;
adv(isnan(adv)) = 0;
adv_h(isnan(adv_h)) = 0;
adv_v(isnan(adv_v)) = 0;
mix(isnan(mix)) = 0;
corr(isnan(corr)) = 0;

%%
% remove correction from advection
adv = adv-corr;
adv_v = adv_v-corr;
%%

% residual
res = adv+div+mix+bio+surf+dilut-tend;

%adv=adv+div;
%clear div
%save('carbBudgetmon' + string(monstart) + 'thru'+string(monend)+'.mat', 'adv', 'adv_h', 'adv_v', 'adv_z', 'adv_y', 'adv_x', 'mix', 'mix_x', 'mix_y', 'mix_z', 'bio', 'surf', 'dilut', 'tend', 'res','div');
%What would I do if I got this to work, tho?? This is wild ...
% Also notice that this is entirely in terms of concentration, but it's
% likely that we'd rather work directly in mol/time... 
% I am very capable of making that change though! 

%% 
% load('carbBudgetmon2thru2.mat', 'adv', 'mix', 'tend','res', 'dilut', 'surf', 'bio');
% 
% check that the terms balance locally
% plot a single time, single location
%%
t=1; x1=520; y1=470;
load('grid.mat', 'RC');
figure(5);
zfig = 1:1:36;
RC = RC(zfig);
hold on
plot(squeeze(tend(x1,y1,zfig,t)),RC); hold on
plot(squeeze(surf(x1,y1,zfig,t)),RC);
plot(squeeze(dilut(x1,y1,zfig,t)),RC);
plot(squeeze(mix(x1,y1,zfig,t)),RC);
plot(squeeze(bio(x1,y1,zfig,t)),RC);
plot(squeeze(adv(x1,y1,zfig,t)),RC);
plot(squeeze(res(x1,y1,zfig,t)),RC,'--k');
ylabel('depth'); xlabel('mol/m3/s');
xlim([-5e-9 5e-9])
legend('tend','surf','dilut','mix','bio','adv','res') % tend and res removed
%%
title('Carbon Profile at Lat = 100 W, 45 S')

%% residual
figure(111)
dep=1;
xlon=1830;
plot(res(xlon, :, dep))
hold on 
plot(adv(xlon,: , dep))
plot(tend(xlon,: , dep)*1)
plot(surf(xlon,: , dep)*1)
plot(dilut(xlon,: , dep)*1)
plot(mix(xlon,:,dep)*1)
%%
legend('res', 'adv', 'tend','surf','dilut','mix')
%%
title('Surface terms thru Rio outflow, unscaled')

