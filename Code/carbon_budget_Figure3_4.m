%%%%%% back cast forest biomass for each year
clear
clc
folder = 'I:\CarbonF\CarbonSink_China\Data\';
folderdir='I:\CarbonF\CarbonSink_China\Output';
grids=readgeoraster([folder,'Grids.tif']);
F=@(p,xdata)p(1)*power((1-exp(-xdata/p(3))),p(2));
agb_int=readmatrix([folder,'agb_20_zonal.txt']);agb_int(:,15:17)=agb_int(:,15:17)/10; %2id 6 max 15 90% 16 95% 17 100%
dis_natural_curve=readmatrix([folder,'REF_natural_curve.csv']); 
dis_plant_curve=readmatrix([folder,'REF_plantation_curve.csv']);
Aff_curve=readmatrix([folder,'AFF_curve.csv']);
CC_NEP=readmatrix([folder,'CcNep.txt']);CC_NEP=CC_NEP(:,[6 11:12]);%ID NEP CC
CC_NEP(:,2)=CC_NEP(:,2)/10;
carbon_total=nan(length(dis_natural_curve(:,1)),30);
%com_2020_biomass=nan(length(dis_natural_curve(:,1)),3);
%com_2020_biomass(:,1)=dis_natural_curve(:,1);
wood=readmatrix([folder,'woodproductf.xlsx']); 
woodp=readmatrix([folder,'woodproduct_INf.xlsx']);
ig=491; % example grid
ID=dis_natural_curve(ig,1);
%fiberboard	Other	paper	plywood/veneer panels	sawnwood
k=[6.557 4.124 3.196 4.161 2.151];o=[3.509 6.242 0.683 9.334 29.982];
k1=[0.042 0.042 0.084 0.042 0.042];k2=[0.025 0.025 0.05 0.025 0.025];
syms x

%%%% disturbance fire + other
[dis,disr]=readgeoraster([folder,'dis',num2str(ID-1),'.tif']);
fire=dis; fire(fire>100)=0;fire(fire>0)=fire(fire>0)+1984;fire=double(fire);
other_dis=dis;other_dis(other_dis<100|other_dis>236)=0;other_dis(other_dis>200)=other_dis(other_dis>200)-200;
other_dis(other_dis>100)=other_dis(other_dis>100)-100;
other_dis(other_dis>0)=other_dis(other_dis>0)+1984; other_dis=double(other_dis);
disT=fire+other_dis;disT=double(disT);
%%%% grid area
wgs84 = wgs84Ellipsoid("km");
lats=repmat(disr.LatitudeLimits(2):-0.000269494585235856:disr.LatitudeLimits(1)+0.000269494585235856,length(dis(1,:)), 1);
lons=repmat([disr.LongitudeLimits(1):0.000269494585235856:disr.LongitudeLimits(2)-0.000269494585235856]',1,length(dis(:,1)));
S=@(x,y)areaquad(x,y,x-0.000269494585235856,y+0.000269494585235856,wgs84);
area_s=arrayfun(S, lats, lons).*100;% ha
%%%%% forest
[F1985,F1985r]=readgeoraster([folder,'F1985_',num2str(ID-1),'.tif']);
F1985=imresize(F1985,[length(dis(:,1)) length(dis(1,:))],'nearest');
%%%%% forest age at 2010
Fage=readgeoraster([folder,'Fage2010_',num2str(ID-1),'.tif']);
Fage=imresize(Fage,[length(dis(:,1)) length(dis(1,:))],'nearest');
%%%%% forest age at 2019
Fage2=readgeoraster([folder,'Age2019_',num2str(ID-1),'.tif']);
Fage2=imresize(Fage2,[length(dis(:,1)) length(dis(1,:))],'nearest');
Fage2(Fage2==65536)=0;
%%%% added disturbance
other_dis2=fire;other_dis2(:,:)=0;%other_dis2(Fage2<35&Fage2>0&disT==0&F1985==2)=2019-Fage2(Fage2<35&Fage2>0&disT==0&F1985==2);other_dis2=double(other_dis2);
length(other_dis2(other_dis2>0))
disTF=disT+other_dis2;
log=other_dis+other_dis2;
%%% fire percent
if isempty(fire(fire>0))
Fire_other_rate=0.064;
else
Fire_other_rate=length(fire(fire>0))/length(disTF(disTF>0)); 
end
%%%%% plantation
if isfile([folder,'PNF2020_',num2str(ID-1),'.tif'])
   plant=readgeoraster([folder,'PNF2020_',num2str(ID-1),'.tif']);
   plant=imresize(plant,[length(dis(:,1)) length(dis(1,:))],'nearest');
   plant(plant==128)=0;plant(plant==2)=0;%plant(plant>0)=plant(plant>0)-3;
else
   plant=dis; plant(:,:)=0; 
end
%%%%% afforestation
if isfile([folder,'AFF',num2str(ID-1),'.tif'])
   afforest=readgeoraster([folder,'AFF',num2str(ID-1),'.tif']);
   afforest=imresize(afforest,[length(dis(:,1)) length(dis(1,:))],'nearest');
   afforest(afforest>10000)=0;afforest(afforest>0)=afforest(afforest>0)-3;
else
   afforest=dis;afforest(:,:)=0; 
end
afforest(disTF>0)=0;afforest(afforest<0)=0;%reforest(isnan(agb30))=nan;
%%%%% deforestation
deforest=readgeoraster([folder,'Deforest',num2str(ID-1),'.tif']);
deforest=imresize(deforest,[length(dis(:,1)) length(dis(1,:))],'nearest');
deforest(deforest>10000)=0;deforest(deforest>1985)=deforest(deforest>1985)+1;deforest(deforest==1985)=1990;% data problem year+1
deforest=double(deforest);
%%% final forest area
forest=F1985;forest(:,:)=0;
forest(F1985==2|disTF>0|deforest>0|afforest>0)=2;
%%%%% agb
agb=readgeoraster([folder,'agb',num2str(ID-1),'.tif'],"OutputType","double");
%agb30=imresize(agb,[length(dis(:,1)) length(dis(1,:))],'nearest');
agb30=agb/10;
if agb30(4,1)==-3276.8
agb30(:,1:3)=agb30(:,4:6);
end
if agb30(1,4)==-3276.8
agb30(1:3,:)=agb30(4:6,:);
end
if agb30(4,end)==-3276.8
agb30(:,end-2:end)=agb30(:,end-5:end-3);
end
if agb30(end,4)==-3276.8
agb30(end-2:end,:)=agb30(end-5:end-3,:);
end
agb30(agb30<0)=nan;
agb30(forest==0)=nan;
%%%% biomass of deforested area
agb30_0=agb30; agb30_0(isnan(agb30_0))=0;
avg = conv2(agb30_0, ones(7)/49, 'same');
Biomass_deforest=tsnanmean(avg(deforest>0&agb30>0&forest==2));
%%%%% biomass of no disturbance
biomass_no_dis=tsnanmean(agb30(disTF==0 & reforest==0 & agb30>0&forest==2));
icurve=1;
%% forest age at 2020 seperate post-dis and expansion (reforestation)
 icurve = 1; % using fitted curve as example, curve 2-5 means the 75%,95% upper and lower range from montecarlo
   if icurve ==1
       dis_natural_p=dis_natural_curve(dis_natural_curve(:,1)==ID,2:4); dis_plant_p=dis_plant_curve(dis_plant_curve(:,1)==ID,2:4);
       Aff_p=Aff_curve(Aff_curve(:,1)==ID,2:4); 
   end
   if icurve ==2
       dis_natural_p=dis_natural_curve(dis_natural_curve(:,1)==ID,23:25); dis_plant_p=dis_plant_curve(dis_plant_curve(:,1)==ID,23:25);
       expant_natural_p=expant_natural_curve(expant_natural_curve(:,1)==ID,23:25); expant_plant_p=expant_plant_curve(expant_plant_curve(:,1)==ID,23:25);
   end
   if icurve ==3
       dis_natural_p=dis_natural_curve(dis_natural_curve(:,1)==ID,26:28); dis_plant_p=dis_plant_curve(dis_plant_curve(:,1)==ID,26:28);
       expant_natural_p=expant_natural_curve(expant_natural_curve(:,1)==ID,26:28); expant_plant_p=expant_plant_curve(expant_plant_curve(:,1)==ID,26:28);
   end
   if icurve ==4 %sd 11-13 75% 17:19 20:22
       dis_natural_p=dis_natural_curve(dis_natural_curve(:,1)==ID,17:19); dis_plant_p=dis_plant_curve(dis_plant_curve(:,1)==ID,17:19);
       expant_natural_p=expant_natural_curve(expant_natural_curve(:,1)==ID,17:19); expant_plant_p=expant_plant_curve(expant_plant_curve(:,1)==ID,17:19);
   end
   if icurve ==5
       dis_natural_p=dis_natural_curve(dis_natural_curve(:,1)==ID,20:22); dis_plant_p=dis_plant_curve(dis_plant_curve(:,1)==ID,20:22);
       expant_natural_p=expant_natural_curve(expant_natural_curve(:,1)==ID,20:22); expant_plant_p=expant_plant_curve(expant_plant_curve(:,1)==ID,20:22);
   end
     
    if isempty(dis_plant_p)
        dis_plant_p=dis_natural_p;
    end
    MbioD=agb_int(agb_int(:,2)==ID,16);
   %MbioD=F(dis_natural_p,300); MbioE=F(expant_natural_p,300);
%    f=figure('visible','on');
%    hold on
%    plot(0:300,F(dis_natural_p,0:300),'linewidth',3);
   %%%%% modify legacy recovery curve
   Fage_t=Fage(agb30<MbioD&reforest==0&disTF==0&forest==2&agb30>0);   Fage_t(isnan(Fage_t))=[];
   Fage_t2=Fage2(agb30<MbioD&reforest==0&disTF==0&forest==2&agb30>0&other_dis2==0);   Fage_t2(isnan(Fage_t2))=[];
   if icurve==1
   Fage2020=mean(Fage_t)+10;
   Fage2020_2=mean(Fage_t2)+1;
   end
   if icurve==2
       Fage2020=prctile(Fage_t,95)+10;
   end
   if icurve==3
       Fage2020=prctile(Fage_t,5)+10;
   end
   if icurve==4
       Fage2020=prctile(Fage_t,75)+10;
   end
   if icurve==5
       Fage2020=prctile(Fage_t,25)+10;
   end
   curvef=0;dis_legacy_p=dis_natural_p; if dis_legacy_p(1)<MbioD; dis_legacy_p(1)=MbioD;end
 %  if dis_legacy_p(2)<2; dis_legacy_p(2)=2; end
  % plot(0:300,F(dis_legacy_p,0:300),'linewidth',3);
   while F(dis_legacy_p,Fage2020_2)>mean(agb30(agb30<MbioD&reforest==0&disTF==0&forest==2&agb30>0))
   curvef=curvef+1;
   dis_legacy_p(3)=dis_legacy_p(3)*1.03;
   end
%    syms x;
%    equation = dis_legacy_p(1) * (1 - exp(-x / dis_legacy_p(3)))^dis_legacy_p(2) == floor(agb30(F1985==2&reforest==0&disT==0));
%    Bage(F1985==2&reforest==0&disT==0) = floor(solve(equation, x, 'Real', true));      
  if floor(MbioD)==MbioD
      MageD_bio=floor(MbioD)-1;
  else
      MageD_bio=floor(MbioD);
  end
   MageD=floor(solve(dis_legacy_p(1)*power((1-exp(-x/dis_legacy_p(3))),dis_legacy_p(2)) == MageD_bio,x,'Real',true)); 
   MageD=double(MageD);MageD=MageD(1);
   growthrateD=[0:MageD+1;F(dis_legacy_p,0:MageD+1)];growthrateD=growthrateD';
   Bage=agb30;Bage(:,:)=0;
   for iage= 1:MageD
      Bage(agb30>=growthrateD(iage,2)&agb30<growthrateD(iage+1,2))=iage;
   end
    Bage(agb30>=floor(MbioD)&reforest==0&disTF==0)=300;
    Bage_dis=Bage;
    Bage_dis(disTF>0&agb30>0)=2020-disTF(disTF>0&agb30>0);
    Bage_dis(reforest>0)=2020-reforest(reforest>0);%the reforest has minus 3 years
    clear("Bage")
%%  carbon dynamic storage

    carbon_dis=nan(length(1986:2020),28);
    carbon_dis(:,1)=1986:2020;
    %Areas:
    carbon_dis(1,2)=sum(sum(area_s(Bage_dis>=300)));%length(Bage_dis(Bage_dis>=300))*9/100 % oldgrowth area
    carbon_dis(2,2)=sum(sum(area_s(Bage_dis<300&disTF==0&reforest==0&F1985==2&agb30>0))); % legacy area
    carbon_dis(3,2)=sum(sum(area_s(fire>0&F1985==2))); % fire area
    carbon_dis(4,2)=sum(sum(area_s(other_dis>0&F1985==2))); % harvest area
    carbon_dis(5,2)=sum(sum(area_s(other_dis2>0&F1985==2&agb30>0))); % other-harvest area
    carbon_dis(6,2)=sum(sum(area_s(reforest>0))); % reforesation area
    carbon_dis(7,2)=sum(sum(area_s(deforest>0))); % deforesation area
    % old growth
    carbon_old_growth=mean(agb30(Bage_dis>=300))*0.5;%%%% mean old growth carbon in 2020
    Biomass_no_dis=mean(agb30(disTF==0&agb30>0));
    carbon_dis(8,2)=carbon_old_growth; % oldgrowth carbon    carbon_dis(8,2)=carbon_old_growth;
    carbon_dis(9,2)=Biomass_no_dis;%*CC_NEP(CC_NEP(:,1)==ID,2)/100;
    Sink_oldgrowth=agb30;Sink_oldgrowth(:,:)=0;
  %  Sink_oldgrowth(Bage_dis>=300)=1;%CC_NEP(CC_NEP(:,1)==ID,2)/100;%.*area_s(Bage_dis>=300).*(2019-1985)./(2019-1985)./area_s(Bage_dis>=300); % Tg C/pixel
%     if ~exist([folderdir,'\Carbon30m\oldgrowth\sinkold',num2str(icurve)], 'dir')
%     mkdir([folderdir,'\Carbon30m\oldgrowth\sinkold',num2str(icurve)]);
%     end
%     geotiffwrite([folderdir,'\Carbon30m\oldgrowth\sinkold',num2str(icurve),'\sinkold_',num2str(ID),'.tif'],Sink_oldgrowth,disr);
   Age1985=Bage_dis-(2020-1985); Age1985(Age1985<0)=0;%floor(tsnanmean(Bage_dis(disTF==0 & reforest==0 & agb30>0&forest==2))-(2020-1985));Age1985(forest==0)=0;
    if tsnanmean(Bage_dis(disTF==0 & reforest==0 & agb30>0&forest==2))-(2020-1985)<0
        carbon_dis(10,2)=1;
    end
   E_fire=agb30;E_fire(:,:)=0; E_fuel=agb30;E_fuel(:,:)=0;E_deforest=agb30;E_deforest(:,:)=0;E_CWD_leg=agb30;E_CWD_leg(:,:)=0;E_CWD=agb30;E_CWD(:,:)=0;E_HWP=agb30;E_HWP(:,:)=0;
    Sink_leg=agb30;Sink_leg(:,:)=0; Sink_dis=agb30;Sink_dis(:,:)=0; Sink_aff=agb30;Sink_aff(:,:)=0;

    
    for iy=1986:2019
      Age=Bage_dis-(2020-iy); Age(forest==0)=0;Age(Age<0)=0;Age_pre=Age;Age_pre(forest==2)=Age(forest==2)-1;Age_pre(Age_pre<0)=0;
              Bio=agb30; 
        Bio(Bage_dis<300&(disTF==0|disTF>iy)&reforest==0&forest==2&agb30>0)=F(dis_legacy_p,Age(Bage_dis<300&(disTF==0|disTF>iy)&reforest==0&forest==2&agb30>0));%legacy
        Bio(reforest<=iy&reforest>0)=F(Aff_p,Age(reforest<=iy&reforest>0));
        Bio(disTF<=iy&disTF>0&plant==0)=F(dis_natural_p,Age(disTF<=iy&disTF>0&plant==0));Bio(disTF<=iy&disTF>0&plant>0)=F(dis_plant_p,Age(disTF<=iy&disTF>0&plant>0));
        Bio(deforest<=iy&deforest>0)=0;Bio(Bage_dis>=300)=agb30(Bage_dis>=300);
        Bio2=Bio; Bio2(disTF>iy)=tsnanmean(avg(disTF>0));
%          if iy==2019
%              if ~exist([folderdir,'\Carbon30m\biomass',num2str(iy),'\bio',num2str(icurve)], 'dir')
%                  mkdir([folderdir,'\Carbon30m\biomass',num2str(iy),'\bio',num2str(icurve)]);
%              end
%            geotiffwrite([folderdir,'\Carbon30m\biomass',num2str(iy),'\bio',num2str(icurve),'\Bio_',num2str(ID),'.tif'],Bio2,disr);
%         end
        Pre_dis_biomass=tsnanmean(Bio2(forest==2&disTF==0&Bio2>0&reforest==0));%tsnanmean(avg(disT>0));
        carbon_dis(iy-1985,4)=Pre_dis_biomass; 
 %  annual area
        carbon_dis(iy-1985,5)=sum(sum(area_s(fire==iy)));
        carbon_dis(iy-1985,6)=sum(sum(area_s(other_dis==iy)));
        carbon_dis(iy-1985,7)=sum(sum(area_s(other_dis2==iy)));
        carbon_dis(iy-1985,8)=sum(sum(area_s(disTF==iy&plant==0)));carbon_dis(iy-1985,9)=sum(sum(area_s(disTF==iy&plant>0)));
 %fire emission
        carbon_dis(iy-1985,10)=-biomass_no_dis*0.5*CC_NEP(CC_NEP(:,1)==ID,3)*carbon_dis(iy-1985,5);%-Pre_dis_biomass*0.5*CC_NEP(CC_NEP(:,1)==ID,3)*carbon_dis(iy-1985,5);
        E_fire_y=agb30;E_fire_y(:,:)=0;E_fire_y(fire==iy)=-Pre_dis_biomass*0.5*CC_NEP(CC_NEP(:,1)==ID,3).*area_s(fire==iy);
        E_fire=E_fire+E_fire_y; 
 %fuelwood  15%residuals other disturbance emission
        carbon_dis(iy-1985,11)=-Pre_dis_biomass*(wood(wood(:,1)==iy,4))*(1-0.15)*0.5*(carbon_dis(iy-1985,6)+carbon_dis(iy-1985,7));
        E_fuel_y=agb30;E_fuel_y(:,:)=0; E_fuel_y(log==iy)=-Pre_dis_biomass*(wood(wood(:,1)==iy,4))*(1-0.15)*0.5.*area_s(log==iy);
        E_fuel=E_fuel+E_fuel_y;
 %legacy CWD carbon emission
        carbon_dis(iy-1985,12)=MbioD*0.5*(Fire_other_rate*(1-CC_NEP(CC_NEP(:,1)==ID,3))+(1-Fire_other_rate)*0.15)*...
            tsnansum(tsnansum((exp(-0.05*(Age(Bage_dis<300&disTF==0&reforest==0)))-exp(-0.05*(Age_pre(Bage_dis<300&disTF==0&reforest==0)))).*area_s(Bage_dis<300&disTF==0&reforest==0)))*0.5;
        E_CWD_y=agb30;E_CWD_y(:,:)=0;
        E_CWD_y(Bage_dis<300&disTF==0&reforest==0)=MbioD*0.5*(Fire_other_rate*(1-CC_NEP(CC_NEP(:,1)==ID,3))+(1-Fire_other_rate)*0.15)*...
            (exp(-0.05*(Age(Bage_dis<300&disTF==0&reforest==0)))-exp(-0.05*(Age_pre(Bage_dis<300&disTF==0&reforest==0))))*0.5.*area_s(Bage_dis<300&disTF==0&reforest==0);%% cwd leg
        E_CWD_leg=E_CWD_leg+E_CWD_y; 
 %new CWD carbon emission
       % CWD_initial=MbioD*0.5*(Fire_other_rate*(1-CC_NEP(CC_NEP(:,1)==ID,3))+(1-Fire_other_rate)*0.15)*tsnanmean(tsnanmean(exp(-0.05*(Age_pre(Bage_dis<300&disTF==0&reforest==0)))));
        carbon_dis(iy-1985,13)=Biomass_deforest*0.5*sum(sum(area_s(deforest==iy&disTF==0)));% deforestation CWD
        carbon_dis(iy-1985,14)=Pre_dis_biomass*0.5*(1-CC_NEP(CC_NEP(:,1)==ID,3))*sum(sum(area_s(fire==iy)));%fire CWD
        carbon_dis(iy-1985,15)=Pre_dis_biomass*0.5*0.15*sum(sum(area_s(log==iy)));%logging CWD
        E_CWD_de=agb30;E_CWD_de(:,:)=0; E_CWD_fire=agb30;E_CWD_fire(:,:)=0;E_CWD_log=agb30;E_CWD_log(:,:)=0;
        if iy>1986
        E_CWD_de(deforest<iy&disTF==0&deforest>0)=Biomass_deforest*0.5.*area_s(deforest<iy&disTF==0&deforest>0).*...
             (exp(-0.05*(iy-deforest(deforest<iy&deforest>0&disTF==0)))-exp(-0.05*(iy-1-deforest(deforest<iy&disTF==0&deforest>0))))*0.5;
        E_CWD_fire(fire<iy&fire>0)=Pre_dis_biomass*0.5*(1-CC_NEP(CC_NEP(:,1)==ID,3)).*area_s(fire<iy&fire>0).*...
             (exp(-0.05*(iy-fire(fire<iy&fire>0)))-exp(-0.05*(iy-1-fire(fire<iy&fire>0))))*0.5;
        E_CWD_log(log<iy&log>0)=Pre_dis_biomass*0.5*0.15.*area_s(log<iy&log>0).*...
             (exp(-0.05*(iy-log(log<iy&log>0)))-exp(-0.05*(iy-1-log(log<iy&log>0))))*0.5;
        carbon_dis(iy-1985,16)=tsnansum(tsnansum(E_CWD_de));% deforestation CWD emission
        carbon_dis(iy-1985,17)=tsnansum(tsnansum(E_CWD_fire));%fire CWD emission
        carbon_dis(iy-1985,18)=tsnansum(tsnansum(E_CWD_log));%logging CWD emission
        end
        E_CWD=E_CWD+E_CWD_de+E_CWD_fire+E_CWD_log;
% HWP %fiberboard	Other	paper	plywood/veneer panels	sawnwood
        HWP_Carbon=(Pre_dis_biomass*carbon_dis(iy-1985,5))*(1-wood(wood(:,1)==iy,4))*(1-0.2)*0.5;
        carbon_dis(iy-1985,19)=tsnansum(tsnansum(Pre_dis_biomass*(1-(wood(wood(:,1)==iy,4)))*(1-0.15)*0.5.*area_s(log==iy)));%HWP carbon pool     
        for iwod=1:5
            HWP=agb30;HWP(:,:)=0;
              f=@(t) t.^(k(iwod)-1)/(gamma(k(iwod)).*o(iwod).^k(iwod)).*exp(-t/o(iwod));
            if iy>1986
                f1=@(t) exp(-t.*k1(iwod));f2=@(t) exp(-t.*k2(iwod)); %aeroric 0.4 anaerobic 0.3
              HWP_y=agb30;HWP_y(:,:)=0;
%               HWP_y(log>0&log<iy)=Pre_dis_biomass*(1-(wood(wood(:,1)==iy,4)))*(1-0.15)*0.5.*area_s(log>0&log<iy)*...
%                   (integral(f,0,iy-1-log(log>0&log<iy))-integral(f,0,iy-log(log>0&log<iy)))
              for idecay = 1:iy-1986
                 % HWP_y=agb30;HWP_y(:,:)=0;
                  if iwod==3
                  HWP_y(log==1986+idecay-1)=-Pre_dis_biomass*(1-(wood(wood(:,1)==1986+idecay-1,4)))*0.5*0.5.*area_s(log==1986+idecay-1)*...
                  woodp(woodp(:,1)==1986+idecay-1,iwod+1)*(integral(f,0,iy-1985-idecay)-integral(f,0,iy-1985-idecay-1));
                  else
                  HWP_y(log==1986+idecay-1)=Pre_dis_biomass*(1-(wood(wood(:,1)==1986+idecay-1,4)))*0.5*0.5.*area_s(log==1986+idecay-1)*...
                      woodp(woodp(:,1)==1986+idecay-1,iwod+1)*...
                      (integral(f,0,iy-1985-idecay)-integral(f,0,iy-1985-idecay-1))*((f1(iy-1985-idecay)-f1(iy-1986-idecay))*0.4+(f2(iy-1985-idecay)-f2(iy-1986-idecay))*0.3);
                  end
              end
              HWP=HWP+HWP_y;
                 carbon_dis(iy-1985,19+iwod)=tsnansum(tsnansum(HWP));
            end
        end 
        E_HWP=E_HWP+HWP;
          
%% regrowth
        %legacy
        Sink_leg_y=agb30;Sink_leg_y(:,:)=0;
        Sink_leg_y(Bage_dis<300&(disTF==0|disTF>iy)&reforest==0&(deforest>iy|deforest==0))=1.116*(F(dis_legacy_p,Age(Bage_dis<300&(disTF==0|disTF>iy)&reforest==0&(deforest>iy|deforest==0)))-...
            F(dis_legacy_p,Age_pre(Bage_dis<300&(disTF==0|disTF>iy)&reforest==0&(deforest>iy|deforest==0))))*0.5.*area_s(Bage_dis<300&(disTF==0|disTF>iy)&reforest==0 & (deforest>iy|deforest==0));% legacy recovery
        carbon_dis(iy-1985,25)=tsnansum(tsnansum(Sink_leg_y));
        Sink_leg=Sink_leg+Sink_leg_y;
        %% disturbance
        Sink_disn_y=agb30;Sink_disn_y(:,:)=0;Sink_disp_y=agb30;Sink_disp_y(:,:)=0;
        Sink_disn_y(disTF>0&disTF<iy&plant==0)=1.12*(F(dis_natural_p,Age(disTF>0&disTF<iy&plant==0))-F(dis_natural_p,Age_pre(disTF>0&disTF<iy&plant==0)))*0.5.*area_s(disTF>0&disTF<iy&plant==0);% disturb product natural recovery
        Sink_disp_y(disTF>0&disTF<iy&plant>0)=1.111*(F(dis_plant_p,Age(disTF>0&disTF<iy&plant>0))-F(dis_plant_p,Age_pre(disTF>0&disTF<iy&plant>0)))*0.5.*area_s(disTF>0&disTF<iy&plant>0);% disturb product plantation recovery
        carbon_dis(iy-1985,26)=tsnansum(tsnansum(Sink_disn_y));
        carbon_dis(iy-1985,27)=tsnansum(tsnansum(Sink_disp_y));
        Sink_dis=Sink_dis+Sink_disp_y+Sink_disn_y; 
        %%% afforestation   
        Sink_aff_y=agb30;Sink_aff_y(:,:)=0;
        Sink_aff_y(reforest>0&reforest<iy)=1.116*(F(Aff_p,Age(reforest>0&reforest<iy))-F(Aff_p,Age_pre(reforest>0&reforest<iy)))*0.5.*area_s(reforest>0&reforest<iy);% expansion plantation recovery
        %Sink_enatural_y(reforest>0&plant==0&reforest<iy)=1.12*(F(expant_natural_p,Age(reforest>0&plant==0&reforest<iy))-F(expant_natural_p,Age_pre(reforest>0&plant==0&reforest<iy)))*0.5.*area_s(reforest>0&plant==0&reforest<iy);% expansion plantation recovery
        carbon_dis(iy-1985,28)=tsnansum(tsnansum(Sink_aff_y));
        Sink_aff=Sink_aff+Sink_aff_y;    
     end
 %%%% end 1986-2019   
     if ~exist([folderdir,'\Carbon30m\E_AGC',num2str(icurve)], 'dir')
              mkdir([folderdir,'\Carbon30m\E_AGC',num2str(icurve)]);
     end
    geotiffwrite([folderdir,'\Carbon30m\E_AGC',num2str(icurve),'\E_AGC_',num2str(ID),'.tif'],(E_fire+E_fuel)/(2019-1986)./area_s,disr);
     if ~exist([folderdir,'\Carbon30m\E_CWD',num2str(icurve)], 'dir')
              mkdir([folderdir,'\Carbon30m\E_CWD',num2str(icurve)]);
     end
    geotiffwrite([folderdir,'\Carbon30m\E_CWD',num2str(icurve),'\E_CWD_',num2str(ID),'.tif'],(E_CWD_leg+E_CWD)/(2019-1986)./area_s,disr);

     if ~exist([folderdir,'\Carbon30m\Sink_leg_newcurve',num2str(icurve)], 'dir')
              mkdir([folderdir,'\Carbon30m\Sink_leg_newcurve',num2str(icurve)]);
     end
    geotiffwrite([folderdir,'\Carbon30m\Sink_leg_newcurve',num2str(icurve),'\Sink_leg_',num2str(ID),'.tif'],Sink_leg/(2019-1986)./area_s,disr);
    if ~exist([folderdir,'\Carbon30m\Sink_dis_newcurve',num2str(icurve)], 'dir')
              mkdir([folderdir,'\Carbon30m\Sink_dis_newcurve',num2str(icurve)]);
     end
    geotiffwrite([folderdir,'\Carbon30m\Sink_dis_newcurve',num2str(icurve),'\Sink_dis_',num2str(ID),'.tif'],Sink_dis/(2019-1986)./area_s,disr);
    if ~exist([folderdir,'\Carbon30m\Sink_aff_newcurve',num2str(icurve)], 'dir')
              mkdir([folderdir,'\Carbon30m\Sink_aff_newcurve',num2str(icurve)]);
     end
    geotiffwrite([folderdir,'\Carbon30m\Sink_aff_newcurve',num2str(icurve),'\Sink_aff_',num2str(ID),'.tif'],Sink_aff/(2019-1986)./area_s,disr);
varNames={'Year','Area:old,leg,fire,log,log2,refo,defo','oldgrowthCarbon','predis-biomass','fireArea','logArea','log2Area','dis_naArea','dis_plantArea',...
        'E_fire','E_fuel','E_CWD_leg','carbon_CWD_def','carbon_CWD_fire','carbon_CWD_log','E_CWD_def','E_CWD_fire','E_CWD_log','carbon_HWP','E_HWP1','E_HWP2',...
        'E_HWP3','E_HWP4','E_HWP5','S_leg','S_disn','S_disp','S_aff'};
T = array2table(carbon_dis);T.Properties.VariableNames = varNames;

if ~exist([folderdir,'\carbon_out',num2str(icurve)], 'dir')
              mkdir([folderdir,'\carbon_out',num2str(icurve)]);
end
writetable(T, [folderdir,'\carbon_out',num2str(icurve),'\carbon_',num2str(ID),'.csv']);
%writematrix(carbon_dis_F,['C:\Project\Carbon\Outputs\carbon_dynamic_HWP_AFF_A\carbon_dis_',num2str(ID),'_',num2str(icurve),'.csv']);
%carbon_total(ig,:)=tsnansum(carbon_dis_F);
