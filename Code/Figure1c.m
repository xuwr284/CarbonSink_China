clear
clc
folder = ['E:\CarbonF\CarbonSink_China\'];
folderdir='Output\AGB_age_curves';
if ~exist([folder,folderdir], 'dir')
mkdir([folder,folderdir]);
end
agb_int=readmatrix([folder,'Data\agb_20_zonal.txt']);% statistic AGB in each 1 by 1 grid; 2id 6 maxAGB 15 90% AGB 16 95% AGB 17 100% AGB 
agb_int(:,6:17)=agb_int(:,6:17)/10;% convert to Mg/ha
agb_int(:,13)=agb_int(:,16);
grids=readgeoraster([folder,'Data\Grids.tif']);% Unique ID for each 1 by 1 grid
dis10=readmatrix([folder,'Data\Disturbance.txt']);%%1ID 2Year 3BA 4grids 5 agb 6CCI AGB 7treecover  8Aforest 9 deforest 10-33plant 11 dis
dis10(:,5)=dis10(:,5)/10;
dis10(dis10(:,5)<0,:)=[]; % remove agb<0 pixels
dis10(dis10(:,11)==0,:)=[];
dis10(dis10(:,2)==2.02,:)=[];% remove pixels disturbed in 2020 year
dis_na=dis10(dis10(:,10)==2,:);
GridsN=unique(dis10(:,4));
dis_natural_curve=readmatrix([folder,'Data\disturbance_curves.csv']); % the original fitting curve
F=@(p,xdata)p(1)*power((1-exp(-xdata/p(3))),p(2));

%% figure
igrid=495; %example grid
ID=dis_natural_curve(igrid,1);
[x,y]=find(grids==ID);
data=dis_na(dis_na(:,4)==ID,:);data(data(:,5)==0,:)=[];
i=0;
while isempty(data)||length(data(:,1))<1000
    i=i+1;x1=x-i;x2=x+i;y1=y-i;y2=y+i;
    if x2>36
       x2=36;
    end
    if x1<1
       x1=1;
    end
    if y2>62
       y2=62;
    end
    if y1<1
       y1=1;
    end
    focalnum=grids(x1:x2,y1:y2); B = focalnum'; B = B(:)';
    data=dis_na(ismember(dis_na(:,4),B),:);data(data(:,5)==0,:)=[];
end
data(:,2)=round(data(:,2)*1000,0);
data(:,1)=2020-data(:,2);
years=unique(data(:,2));
% remove errors
Bio_t=nan(length(years(:,1)),3);
Bio_t(:,1)=years;
for iy=1:length(years(:,1))
   Bio_t(iy,2)=mean(data(data(:,2)==years(iy),5));
   Bio_t(iy,3)=median(data(data(:,2)==years(iy),5));
end
Bio_t(:,4)=2020-Bio_t(:,1);
Bio_t(Bio_t(:,1)==2020,:)=[];
p1=dis_natural_curve(igrid,4:6);
squares=((F(p1,data(:,1)) - data(:,5))./mean(data(:,5))).^2;
data(squares>=1.5,:)=[];
%% end remove errors
years=unique(data(:,2));
Bio_t=nan(length(years(:,1)),3);
Bio_t(:,1)=years;
for iy=1:length(years(:,1))
   Bio_t(iy,2)=mean(data(data(:,2)==years(iy),5));
   Bio_t(iy,3)=median(data(data(:,2)==years(iy),5));
end
Bio_t(:,4)=2020-Bio_t(:,1);
Bio_t(Bio_t(:,1)==2020,:)=[];
if i>0
Bio_t(Bio_t(:,3)>agb_int(agb_int(:,2)==ID,13),:)=[];
end
bio200=agb_int(agb_int(:,2)==ID,13);Bio_t_o=[Bio_t;[0,bio200,bio200,100];[0,bio200,bio200,150]];
t_o=Bio_t_o(1:end,4);y_o=Bio_t_o(1:end,3);
c_init = max(Bio_t(:,4))/((max(Bio_t(:,3)) - min(Bio_t(:,3))) / (max(Bio_t(:,4)) - min(Bio_t(:,4))));
p0 = dis_natural_curve(igrid,4:6);%[agb_int(agb_int(:,2)==ID,16) 3 c_init];
Fsumsquares = @(p)sum((F(p,t_o) - y_o).^2);
options = optimset('Display','iter','FunValCheck','on', ...
    'MaxFunEvals',2000,'MaxIter',2000, ...
    'TolFun',1e-6,'TolX',1e-6);
paramslb = [0 1.5 c_init];  % lower bound
paramsub = [agb_int(agb_int(:,2)==ID,13)*1.1 Inf Inf];
[p,resnorm,residual,existflag,output,lambda,jacobian] = ...
    lsqcurvefit(F,p0,t_o,y_o,paramslb,paramsub,options); 
%rsq=corr(Bio_t(1:end,3),F(p,Bio_t(1:end,4)))^2;
%predfun = @(beta)p(1).*(1-exp(-1/beta(2).*t_o)).^beta(1);
predfun = @(param)param(1).*(1-exp(-1/param(3).*t_o)).^param(2);
rsq=corr(y_o,F(p,t_o))^2;
    %rsq2=1-sum((y-F(p,t)).^2)/sum((y-mean(y)).^2);
sse=nansum((predfun(p)-y_o).^2);
rmse1=sqrt(sse/(length(t_o)-1));
rsquare=1-sse/sum((Bio_t(1:end,3)-mean(Bio_t(1:end,3))).^2);
rmse=sqrt(mean((Bio_t(1:end,3)-F(p,Bio_t(1:end,4))).^2));
 sdr = sqrt(sum((y_o - predfun(p)).^2)/(length(t_o)-1));
 J = f_jacobianest(predfun,p);
 Sigma = sdr^2*inv(J'*J);
 se = sqrt(diag(Sigma))';
 se(p-se<0)=p(p-se<0);
xmax=100;ymax=agb_int(agb_int(:,2)==ID,13)*1.2;
xall=[0:100]';
[~,q25,q75] = f_montecarlo(p,se,xall);
p11(:,1)=q25;
p11(:,2)=q75;
p11(:,3)=0:xmax;
curve(igrid,1)=ID;curve(igrid,2:4)=p;curve(igrid,5)=rsq;curve(igrid,6)=rmse;
f = figure('visible','on');
 hold on
plot(0:100, F(p, 0:100), 'linewidth',1,Color=[0 0.4470 0.7410]);
hold on
X_plot  = [xall', fliplr((xall)')];
Y_plot  = [p11(:,1)', fliplr(p11(:,2)')];
fill(X_plot, Y_plot , 1,....
       'facecolor',[0 0.4470 0.7410], ...
       'edgecolor','none', ...
       'facealpha', 0.3);
hold on
plot([0;Bio_t(1:end,4)],[0;Bio_t(1:end,3)],'DisplayName','AGB', 'MarkerSize',5,'Marker','o',...
   'Color', [0 0.4470 0.7410], 'LineStyle','none');
hold on
plot(50,agb_int(agb_int(:,2)==ID,13),'DisplayName','AGB', 'MarkerSize',20,'Marker','.',...
        'LineStyle','none',...
        'Color',[0,90/255,0]);
hold on
line([0, xmax], [p(1), p(1)], 'Color', [0 0.4470 0.7410], 'LineStyle','-');
legend( 'fitted curve','Uncertainty range','AGB-age pairs','AGBmax','AGB_eq', 'Location', 'southeast' ); % 'GlobBiomass2010',
 ylim([0 ymax])
%text(30,ymax*0.15,sprintf('R^2 = %3.2g\nRMSE = %2.4g',rsq,rmse),'FontSize',9)


%%%% afforestation 
disaff=readmatrix([folder,'Data\Aff.txt']);%%1ID 2Year 3BA 4grids 5agb 6Aforest 7 deforest 8-33plant 9 dis
disaff(:,5)=disaff(:,5)/10;
disaff(disaff(:,5)<0,:)=[];disaff(disaff(:,6)<0,:)=[];disaff(disaff(:,7)>0,:)=[];disaff(disaff(:,9)>0,:)=[];
disaff(disaff(:,2)==2.02,:)=[];disaff(:,2)=disaff(:,2)-3/1000; % detect 3y early
dis_aff_curve=readmatrix([folder,'Data\Afforestation_curves.csv']); 
data=disaff(disaff(:,4)==ID,:);data(data(:,5)==0,:)=[];
i=0;
years=unique(data(:,2));
Bio_t=nan(length(years(:,1)),3);
Bio_t(:,1)=years;
for iy=1:length(years(:,1))
   Bio_t(iy,2)=mean(data(data(:,2)==years(iy),5));
   Bio_t(iy,3)=median(data(data(:,2)==years(iy),5));
end
bio300=min(agb_int(agb_int(:,2)==ID,13),max(Bio_t(:,3)));
%bio300=min(agb_int(agb_int(:,2)==ID,16),max(data(:,5)));
while isempty(data)||length(data(:,1))<100
    i=i+1;x1=x-i;x2=x+i;y1=y-i;y2=y+i;
    if x2>36
       x2=36;
    end
    if x1<1
       x1=1;
    end
    if y2>62
       y2=62;
    end
    if y1<1
       y1=1;
    end
    focalnum=grids(x1:x2,y1:y2); B = focalnum'; B = B(:)';
    data=disaff(ismember(disaff(:,4),B),:);data(data(:,5)==0,:)=[];
end
data(:,2)=round(data(:,2)*1000,0);
data(:,1)=2020-data(:,2);
years=unique(data(:,2));
% remove errors
Bio_t=nan(length(years(:,1)),3);
Bio_t(:,1)=years;
for iy=1:length(years(:,1))
   Bio_t(iy,2)=mean(data(data(:,2)==years(iy),5));
   Bio_t(iy,3)=median(data(data(:,2)==years(iy),5));
end
Bio_t(:,4)=2020-Bio_t(:,1);
Bio_t(Bio_t(:,1)==2020,:)=[];
p1=dis_aff_curve(dis_aff_curve(:,1)==ID,4:6);
squares=((F(p1,data(:,1)) - data(:,5))./mean(data(:,5))).^2;
data(squares>=1.5,:)=[];
%% end remove errors
years=unique(data(:,2));
Bio_t=nan(length(years(:,1)),3);
Bio_t(:,1)=years;
for iy=1:length(years(:,1))
   Bio_t(iy,2)=mean(data(data(:,2)==years(iy),5));
   Bio_t(iy,3)=median(data(data(:,2)==years(iy),5));
end
Bio_t(:,4)=2020-Bio_t(:,1);
Bio_t(Bio_t(:,1)==2020,:)=[];
if i>0
Bio_t(Bio_t(:,3)>agb_int(agb_int(:,2)==ID,13),:)=[];
end
if isempty(bio300)
   bio300=min(agb_int(agb_int(:,2)==ID,13),max(Bio_t(:,3)));
end
bio200=agb_int(agb_int(:,2)==ID,13);%Bio_t_o=Bio_t;%[Bio_t;[0,bio200,bio200,200]];
Bio_t_o=[Bio_t;[0,bio300,bio300,50];[0,bio300,bio300,100]];
t_o=Bio_t_o(1:end,4);y_o=Bio_t_o(1:end,3);
%c_init = max(Bio_t(:,4))/((max(Bio_t(:,3)) - min(Bio_t(:,3))) / (max(Bio_t(:,4)) - min(Bio_t(:,4))));
p0 = p1;%[agb_int(agb_int(:,2)==ID,16) 3 c_init];
Fsumsquares = @(p)sum((F(p,t_o) - y_o).^2);
options = optimset('Display','iter','FunValCheck','on', ...
    'MaxFunEvals',2000,'MaxIter',2000, ...
    'TolFun',1e-6,'TolX',1e-6);
paramslb = [0 0.5 20];  % lower bound
paramsub = [agb_int(agb_int(:,2)==ID,13)*1.1 Inf Inf];
[p,resnorm,residual,existflag,output,lambda,jacobian] = ...
    lsqcurvefit(F,p0,t_o,y_o,paramslb,paramsub,options); 
%rsq=corr(Bio_t(1:end,3),F(p,Bio_t(1:end,4)))^2;
%predfun = @(beta)p(1).*(1-exp(-1/beta(2).*t_o)).^beta(1);
predfun = @(param)param(1).*(1-exp(-1/param(3).*t_o)).^param(2);
rsq=corr(y_o,F(p,t_o))^2;
    %rsq2=1-sum((y-F(p,t)).^2)/sum((y-mean(y)).^2);
sse=nansum((predfun(p)-y_o).^2);
rmse1=sqrt(sse/(length(t_o)-1));
rsquare=1-sse/sum((Bio_t(1:end,3)-mean(Bio_t(1:end,3))).^2);
rmse=sqrt(mean((Bio_t(1:end,3)-F(p,Bio_t(1:end,4))).^2));
 sdr = sqrt(sum((y_o - predfun(p)).^2)/(length(t_o)-1));
 J = f_jacobianest(predfun,p);
 Sigma = sdr^2*inv(J'*J);
 se = sqrt(diag(Sigma))';
 se(p-se<0)=p(p-se<0);
xmax=100;ymax=agb_int(agb_int(:,2)==ID,13)*1.2;
xall=[0:100]';
[~,q25,q75] = f_montecarlo(p,se,xall);
p11(:,1)=q25;
p11(:,2)=q75;
p11(:,3)=0:xmax;
curve(igrid,1)=ID;curve(igrid,2:4)=p;curve(igrid,5)=rsq;curve(igrid,6)=rmse;
plot(0:100, F(p, 0:100), 'linewidth',1,Color=[0.9290 0.6940 0.1250]);
hold on
X_plot  = [xall', fliplr((xall)')];
Y_plot  = [p11(:,1)', fliplr(p11(:,2)')];
fill(X_plot, Y_plot , 1,....
       'facecolor',[0.9290 0.6940 0.1250], ...
       'edgecolor','none', ...
       'facealpha', 0.3);
hold on
plot([0;Bio_t(1:end,4)],[0;Bio_t(1:end,3)],'DisplayName','AGB', 'MarkerSize',5,'Marker','o',...
   'Color', [0.9290 0.6940 0.1250], 'LineStyle','none');
line([0, xmax], [p(1), p(1)], 'Color', [0.9290 0.6940 0.1250], 'LineStyle','--');
%legend( 'fitted curve','Uncertainty range','AGB points','95% AGB (region)','AGB_eq', 'Location', 'southeast' ); % 'GlobBiomass2010',
%text(20,ymax*0.25,sprintf('R^2 = %3.2g\nRMSE = %2.4g',rsq,rmse),'FontSize',9)
lgd = findobj('type', 'legend');
%delete(lgd)
hold off
saveas(f,[folder,folderdir,'\curve_',num2str(ID),'.jpg'])
close(f)










