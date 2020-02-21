% Code for Running Spatial Biomass Model

% Cara Scalpone, January 2020
% Code structure from Jeremy Testa and Neil Ganju

%% Required input variables (created from raw data files)
% ATEMP                -- daily air temperature grid
% mgrid, bgrid         -- slode grid and offset grid from regression analysis
% dailyPARgrid         -- daily photosynthetically active radiation grid
% baymask              -- masks area outside of BB-LEH
% DINsed_uM, DINwc_uM  -- daily DIN in sediment and water column
% h                    -- depth grid

%% Set calibrated parameters
kmAG1 = 0.003; % aboveground mortality constant before JDm (d-1)
km1 = 0.01;    % aboveground mortality constant before JDm (d-1)
AGBr = 12;     % initial AGB on May 1
BGBr = 18;     % initial BGB on May 1

%% Set arrays for output
AGBmap=nan(365,800,160); %aboveground biomass
BGBmap=nan(365,800,160); %belowground biomass
        
h(h<0.1)=NaN;    %eliminate intertidal areas, land
h(mgrid==0)=NaN;     %elimate all land areas if not already eliminated with bay mask
        		% i.e. locations where slope == 0, where there is no water
                
%% Run SAV model with spatial inputs
        for j=1:800
            for k=1:160
                
                if isnan(h(j,k)) || isnan(baymask(j,k)) %mask Great Bay, outgoing channels, and land
                    AGBmap(:,j,k)=NaN;
                    BGBmap(:,j,k)=NaN;
                    
                else
                    [AGB,BGB]=bio_veg_jcj_2019(h(j,k), mgrid(j,k),bgrid(j,k),sq(ATEMP(:,j,k)), ...
                        sq(dailyPARgrid(:,j,k)), DINsed_uM, DINwc_uM, kmAG1, km1, AGBr, BGBr);
                    
                    AGBmap(:,j,k)=AGB; % aboveground biomass map
                    BGBmap(:,j,k)=BGB; % belowground biomass map

                end
                
            end
            
        end

