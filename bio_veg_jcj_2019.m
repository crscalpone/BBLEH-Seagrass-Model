% Estuarine SAV Model adapted from Straub et al. (2015)

% Cara Scalpone, September 2019
% Code structure from Jeremy Testa, 2015, Chesapeake Biological Laboratory

% Model Inputs:

% z : depth (m)
% temp_mult : slope derived from Air Temp ~ Water Temp
%             regressions
% temp_offset : off-set derived from Air Temp ~ Water Temp
%               regressions
% ATEMP : air temperature (C)
% PARs : photosynthetically active radiation at the surface (uE m-2 s-1)
% DINsed : dissolved inorganic nitrogen in sediment (uM)
% DINwc : dissolved inorganic nitrogen in water column (uM)
% kmAG1 : aboveground mortality constant before JDm (d-1)
% km1 : belowground mortality constant before JDm (d-1)

function [AGB,BGB,PPz,Pmaxz,LLIMz,NLIMz]=bio_veg_jcj_2019(z,temp_mult,...
    temp_offset,ATEMP,PARs,DINsed,DINwc,kmAG1,km1,AGBr,BGBr) 
%% Interpolate forcing files to daily values
time=1:1:365;

% Define JD for May 1 model initiation
JD = zeros(1,365);     
JD(:,1:245) = 121:365; 
JD(:,246:365) = 1:120; 

WTEMP=ATEMP.*temp_mult + temp_offset;

%% Environmental Parameters

kd = 1.9;      % Light attenuation coefficient 
sr = 0.1;      % Surface reflectance of light

%% SAV Parameters

knsed  = 28.558;  %Half saturation coefficient for sediment N uptake (uM)
knwcz  = 7.319;   %Half saturation coefficient for water-column N uptake (uM)
Toptz  = 22.5;    %Optimum temperature for zostera growth (deg C)
kiz    = 57.5;    %Half-saturation for light 
rc     = 0.00005; %Root respiration constant and rate when ATEMP = Topt (1/d)
rrc    = 1.25;    %Root respiration constant (1/d)
Downt  = 0.4;     %Downward translocation cofficient (1/d)
trns   = 0.02;    %Upward translocation cofficient (1/d)
tcrit  = 10;      %Critical temperature for development of aboveground biomass (deg C)  
km2    = 0.031;   %Belowground mortality constant after JDm (1/d)
kmAG2  = 0.034;   %Aboveground mortality constant after JDm (1/d)
JDm    = 166;     %Julien Day of mortality increase (June 15)

%% Initial Conditions and set arrays

PARd    = ones(length(time),1);
LLIMz   = ones(length(time),1);
knt     = ones(length(time),1);
NLIMz   = ones(length(time),1);
Pmaxz   = ones(length(time),1);
PhotoPD = ones(length(time),1);
PRz     = ones(length(time),1);
PPz     = ones(length(time),1);
AGM     = ones(length(time),1);
AGR     = ones(length(time),1);
AGBG    = ones(length(time),1);
BGAG    = ones(length(time),1);
BGR     = ones(length(time),1);
BGM     = ones(length(time),1);
AGB     = ones(length(time),1);
BGB     = ones(length(time),1);

AGB(1,1) = AGBr;
BGB(1,1) = BGBr;

%% Run SAV model for a year
for i = 2:length(time)
    
    if WTEMP(i) > tcrit %Prevents upward translocation until temp is high enough for germination
        if AGB(i-1,1) < 0.44 %If density is greater than 0.44 g C m-2, no germination occurs
            thresh = 1;
        else
            thresh = 0;
        end
    else
        thresh = 0;
    end
    
    
    if JD(i) < JDm %Changes mortality constant
            kmAG = kmAG1;
            km = km1;
        else
            kmAG = kmAG2;
            km = km2;
    end
 
    
    %Light, nutrients for Zostera
    
    PARd(i,1)   = (PARs(i,1)*(1-sr))*exp(-kd*z);                                         %PAR at depth of meadow
    LLIMz(i,1)  = PARd(i,1)/(PARd(i,1)+kiz);
    knt(i,1)    = knwcz/knsed;                                                           %water-column/sediment DIN factor
    NLIMz(i,1)  = (DINwc(i,1)+(knt(i,1)*DINsed(i,1)))/(knwcz+(DINwc(i,1)+(knt(i,1)*DINsed(i,1)))); %WC+SED Nitrogen limitation factor for Zostera
    
    
    %Maximum production rate for Zostera
    Pmaxz(i,1) = (0.0948+0.0309.*exp(-0.5.*(((WTEMP(i,1)-Toptz)./3.2964)^2)));           %Pmax Zostera
    PhotoPD(i,1) = (12-(2.5*cos(2*pi*(JD(i)-354)/365)))/24;
    
    %Run SAV model
    PRz(i,1)     = Pmaxz(i,1)*PhotoPD(i,1)*min(LLIMz(i,1),NLIMz(i,1));                   %Primary Production Rate for Zostera
    PPz(i,1)     = AGB(i-1,1)*PRz(i,1);                                                  %Addition from Primary Production
    AGM(i,1)     = AGB(i-1,1)*kmAG;
    AGR(i,1)     = AGB(i-1,1)*(PRz(i,1)*(0.00317*(WTEMP(i,1)+0.105)+exp(0.135*WTEMP(i,1)-10.1))); %Zostera Above ground respiration
    AGBG(i,1)    = PPz(i,1)*Downt;                                                       %Translocation of above ground biomass to below ground
    BGAG(i,1)    = BGB(i-1,1)*trns*thresh;                                               %Translocation of below ground biomass to above ground (because no reproduction formulation)
    BGR(i,1)     = BGB(i-1,1)*rc*rrc^(WTEMP(i,1)-Toptz);                                 %Below ground biomass respiration
    BGM(i,1)     = km*BGB(i-1,1);                                                        %Below ground biomass mortality
    AGB(i,1)     = AGB(i-1,1)+PPz(i,1)+BGAG(i,1)-AGM(i,1)-AGR(i,1)-AGBG(i,1);            %Compute new AGB Biomass (g C m-2)
    BGB(i,1)     = BGB(i-1,1)+AGBG(i,1)-BGAG(i,1)-BGM(i,1)-BGR(i,1);                     %Compute new BGB Biomass  (g C m-2)
    
    
end