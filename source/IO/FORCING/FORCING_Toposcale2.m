

classdef FORCING_Toposcale2 < matlab.mixin.Copyable
    
    properties
        forcing_index
        DATA            % forcing data time series
        TEMP            % forcing data interpolated to a timestep
        PARA            % parameters
        STATUS         
        CONST
    end
    
    
    methods
        
        
        function forcing = provide_PARA(forcing)         

            forcing.PARA.filename = [];   %filename of Matlab file containing forcing data
            forcing.PARA.forcing_path = [];
            forcing.PARA.start_time = []; % start time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.end_time = [];   % end time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.forcing_timestep = []; %in hours
            
            forcing.PARA.rain_fraction = [];  %rainfall fraction assumed in sumulations (rainfall from the forcing data file is multiplied by this parameter)
            forcing.PARA.snow_fraction = [];  %snowfall fraction assumed in sumulations (snowfall from the forcing data file is multiplied by this parameter)
            forcing.PARA.precipitation_lapse_rate = [];
            
            forcing.PARA.all_rain_T = [];
            forcing.PARA.all_snow_T = [];
            forcing.PARA.slope_angle = []; %slope angle in degrees
            forcing.PARA.aspect = []; %aspect of the slope in degrees
            forcing.PARA.albedo_surrounding_terrain = [];
            forcing.PARA.sky_view_factor = []; %sky view factor (0.5 for vertical rock walls)
            forcing.PARA.horizon_angle_bins = [];
            forcing.PARA.horizon_angle = [];
            
            forcing.PARA.heatFlux_lb = [];  % heat flux at the lower boundary [W/m2] - positive values correspond to energy gain
            forcing.PARA.min_height_pl_above_ground = [];
            forcing.PARA.number_of_grid_cells = [];
        end
        
        
        
        function forcing = provide_CONST(forcing)
            forcing.CONST.sigma = [];
            forcing.CONST.R_spec = [];
            forcing.CONST.g = [];
            forcing.CONST.rho_w = [];
        end
        
        function forcing = provide_STATVAR(forcing)
            
        end

        
        function forcing = finalize_init(forcing, tile)
            disp('loading forcing file')
            load([forcing.PARA.forcing_path forcing.PARA.filename]);
            
            disp('processing forcing data')
            
            kNN = forcing.PARA.number_of_grid_cells;
            
            [lat, lon] = meshgrid(era.lat, era.lon);
            target_lat = tile.PARA.latitude;
            target_lon = tile.PARA.longitude;
            target_altitude = tile.PARA.altitude + forcing.PARA.min_height_pl_above_ground;
            surface_altitude = tile.PARA.altitude;
            slope = forcing.PARA.slope_angle;
            aspect = forcing.PARA.aspect;
            horizon_angle_bins = deg2rad(forcing.PARA.horizon_angle_bins); % Not important, but these are typical horizon elevation bins (south=0) from TopoPAR.
            horizon_angle = deg2rad(forcing.PARA.horizon_angle); % Elevation angles in radians (just 0s for flat horizon)
            sky_view_factor = forcing.PARA.sky_view_factor;
            precipitation_lapse_rate = forcing.PARA.precipitation_lapse_rate;
            
            R_spec = forcing.CONST.R_spec;  % Gas constant for dry air [JK^-1kg^-1]
            g = forcing.CONST.g; % Acceleration of gravity [ms^-1]
            eps0=0.622; % Ratio of molecular weight of water and dry air [-]
            S0=1370; % Solar constat (total TOA solar irradiance) [Wm^-2] used in ECMWF's IFS
            sbc = forcing.CONST.sigma;
            
            d=sqrt((target_lat - lat).^2 + (target_lon-lon).^2);
            idw=1./(d.^2);
            [~,rankis]=sort(d(:),1,'ascend');
            idw(rankis>kNN)=0;
            idw=idw./sum(idw(:),1);
            idw(isnan(idw)) = 1;
            
            %ccordinates of the closest points in xy
            coords =[mod((rankis(1:kNN,1)+1), size(era.lat,1))+1  floor(((rankis(1:kNN,1))+1)./size(era.lon,1))];
            
            %find correct pressure level for T, q and wind u and v
            
            %0. for each pressure level and ERA grid cell,  eliminate all
            %values of pressure levels below the orography
            %make array 2x2xtime with first level above
            
            layer_above_orography = zeros(size(era.Z,1),size(era.Z,2),size(era.Z,4));
            for i=1:size(era.Z,3)
                layer_above_orography = layer_above_orography + double(squeeze(era.Z(:,:,i,:)) > repmat(era.Zs,1,1, size(era.Z,4)));
            end
            layer_above_orography = double(layer_above_orography);
            %1. do the same for the true altitude plus the target offset, find layers below and above
            
            
            layer_above_target_altitude = zeros(size(era.Z,1),size(era.Z,2),size(era.Z,4));
            layer_above_surface_altitude = zeros(size(era.Z,1),size(era.Z,2),size(era.Z,4));
            for i=1:size(era.Z,3)
                layer_above_target_altitude = layer_above_target_altitude + double(squeeze(era.Z(:,:,i,:)) > repmat(target_altitude,size(era.Z,1),size(era.Z,2),size(era.Z,4)));
                layer_above_surface_altitude = layer_above_surface_altitude + double(squeeze(era.Z(:,:,i,:)) > repmat(surface_altitude,size(era.Z,1),size(era.Z,2),size(era.Z,4)));
            end
            layer_above_target_altitude = double(layer_above_target_altitude);
            layer_above_surface_altitude = double(layer_above_surface_altitude);
            
            %layer above target altitude must also be above orography
            layer_above_target_altitude = min(layer_above_target_altitude, layer_above_orography);
            
            %interpolation is possible
            layer_below_is_above_orography = (layer_above_target_altitude + 1) <= layer_above_orography  & (layer_above_target_altitude + 1) <= size(era.Z,3);
            
            
            %condition: if layer_below_is_above_orography, then interpolate
            %to correct altitude, otherwise use layer_above_target_altitude
            altitude = 0;
            temperature = 0;
            wind = 0;
            q = 0;
            p = 0;
            p_forcing_height = 0;
            
            %loop over relevant coordinates
            for i=1:size(coords,1)
                la_ab_ta_alt = squeeze(layer_above_target_altitude(coords(i,1), coords(i,2),:));
                la_be_ta_alt = min((la_ab_ta_alt + 1), size(era.Z,3));
                la_ab_su_alt = squeeze(layer_above_surface_altitude(coords(i,1), coords(i,2),:));
                
                la_be_is_ab_oro = squeeze(layer_below_is_above_orography(coords(i,1), coords(i,2),:));
                %loop over pressure levels
                la_ab_ta_alt_matrix = zeros(size(era.Z,3), size(era.Z,4));
                la_be_ta_alt_matrix = zeros(size(era.Z,3), size(era.Z,4));
                la_be_ta_alt_matrix_p = zeros(size(era.Z,3), size(era.Z,4));
                for j=1:size(era.Z,3)
                    la_ab_ta_alt_matrix(j,:) = double(la_ab_ta_alt == j);
                    la_be_ta_alt_matrix(j,:) = double(la_be_ta_alt == j);
                    la_ab_ta_alt_matrix_p(j,:) = double(min(la_ab_su_alt, size(era.Z,3)-1) == j);
                    la_be_ta_alt_matrix_p(j,:) = double(min(la_ab_su_alt+1, size(era.Z,3)) == j);
                end
                z_above_target_altitude = double(sum(squeeze(era.Z(coords(i,1), coords(i,2),:, :)) .* la_ab_ta_alt_matrix, 1));
                z_below_target_altitude = double(sum(squeeze(era.Z(coords(i,1), coords(i,2),:, :)) .* la_be_ta_alt_matrix, 1));
                T_above_target_altitude = double(sum(squeeze(era.T(coords(i,1), coords(i,2),:, :)) .* la_ab_ta_alt_matrix, 1));
                T_below_target_altitude = double(sum(squeeze(era.T(coords(i,1), coords(i,2),:, :)) .* la_be_ta_alt_matrix, 1));
                q_above_target_altitude = sum(squeeze(era.q(coords(i,1), coords(i,2),:, :)) .* la_ab_ta_alt_matrix, 1);
                q_below_target_altitude = sum(squeeze(era.q(coords(i,1), coords(i,2),:, :)) .* la_be_ta_alt_matrix, 1);
                p_above_target_altitude = sum(repmat((era.p)', 1, size(era.Z,4)) .* la_ab_ta_alt_matrix, 1);
                p_below_target_altitude = sum(repmat((era.p)', 1, size(era.Z,4)) .* la_be_ta_alt_matrix, 1);
                wind_above_target_altitude = sqrt((sum(squeeze(era.u(coords(i,1), coords(i,2),:, :)) .* la_ab_ta_alt_matrix, 1)).^2 +...
                    (sum(squeeze(era.v(coords(i,1), coords(i,2),:,:)) .* la_ab_ta_alt_matrix, 1)).^2);
                wind_below_target_altitude = sqrt((sum(squeeze(era.u(coords(i,1), coords(i,2),:, :)) .* la_be_ta_alt_matrix, 1)).^2 +...
                    (sum(squeeze(era.v(coords(i,1), coords(i,2),:, :)) .* la_be_ta_alt_matrix, 1)).^2);
                
                z_above_surface_altitude = sum(squeeze(era.Z(coords(i,1), coords(i,2),:, :)) .* la_ab_ta_alt_matrix_p, 1);
                z_below_surface_altitude = sum(squeeze(era.Z(coords(i,1), coords(i,2),:, :)) .* la_be_ta_alt_matrix_p, 1);
                p_above_surface_altitude = sum(repmat((era.p)', 1, size(era.Z,4)).* la_ab_ta_alt_matrix_p, 1);
                p_below_surface_altitude = sum(repmat((era.p)', 1, size(era.Z,4)).* la_be_ta_alt_matrix_p, 1);
                
                fraction_above = (target_altitude - z_below_target_altitude)./(z_above_target_altitude - z_below_target_altitude);
                fraction_above(isnan(fraction_above)) = 0;
                fraction_above(fraction_above == Inf) = 0;
                fraction_above(fraction_above == -Inf) = 0;
                
                fraction_above_p = (surface_altitude - z_below_surface_altitude)./(z_above_surface_altitude - z_below_surface_altitude);
                fraction_above_p(isnan(fraction_above_p)) = 0;
                fraction_above(fraction_above_p == Inf) = 0;
                fraction_above(fraction_above_p == -Inf) = 0;
                
                temperature_interp = (1-fraction_above) .* T_below_target_altitude + fraction_above .* T_above_target_altitude;
                wind_interp = (1-fraction_above) .* wind_below_target_altitude + fraction_above .* wind_above_target_altitude;
                q_interp = (1-fraction_above) .* q_below_target_altitude + fraction_above .* q_above_target_altitude;
                p_interp = (1-fraction_above) .* p_below_target_altitude + fraction_above .* p_above_target_altitude;
                p_interp_surface = (1-fraction_above_p) .* p_below_surface_altitude + fraction_above_p .* p_above_surface_altitude;
                altitude = altitude + idw(coords(i,1), coords(i,2)) .* (double(la_be_is_ab_oro) .* target_altitude + double(~la_be_is_ab_oro) .* z_above_target_altitude');
                temperature = temperature + idw(coords(i,1), coords(i,2)) .* (double(la_be_is_ab_oro) .* temperature_interp' + double(~la_be_is_ab_oro) .* T_above_target_altitude');
                wind = wind + idw(coords(i,1), coords(i,2)) .* (double(la_be_is_ab_oro) .* wind_interp' + double(~la_be_is_ab_oro) .* wind_above_target_altitude');
                q = q + idw(coords(i,1), coords(i,2)) .* (double(la_be_is_ab_oro) .* q_interp' + double(~la_be_is_ab_oro) .* q_above_target_altitude');
                p_forcing_height = p_forcing_height + idw(coords(i,1), coords(i,2)) .* (double(la_be_is_ab_oro) .* p_interp' + double(~la_be_is_ab_oro) .* p_above_target_altitude');
                p = p + idw(coords(i,1), coords(i,2)) .* p_interp_surface';
            end
            
            temperature = temperature - 273.15;
            
            %% precipitation
            precipitation = 0;
            altitude_ERA_orography = 0;
            
            %loop over relevant coordinates
            for i=1:size(coords,1)
                precipitation = precipitation + idw(coords(i,1), coords(i,2)) .* squeeze(era.P(coords(i,1), coords(i,2),:));
                altitude_ERA_orography = altitude_ERA_orography + idw(coords(i,1), coords(i,2)) .* era.Zs(coords(i,1), coords(i,2));
            end
            
            precipitation  = 24.* precipitation .* precipitation_lapse_rate.^((surface_altitude-altitude_ERA_orography)./100); %convert from mm/h to mm/day
            precipitation(precipitation<0) = 0;
            
            
            %% Longwave routine.
            

            % Useful in-line functions.
            K2C=@(Tk) Tk-273.15;
            q2w=@(q) 0.5.*(1-sqrt(1-4.*q)); % Mixing ratio from specific humidity based on 2nd order Taylor series expansion.
            wp2e=@(w,p) 0.5.*p.*(-1+sqrt(1+4.*w./eps0)); % Vapor pressure from mixing ratio based on 2nd order Taylor series expansion.
            
            % AERK from Alduchov and Eskridge (1996).
            A1=17.625; B1=243.04; C1=610.94;
            Magnus=@(tc) C1.*exp(A1.*tc./(B1+tc)); % A version of the Magnus formula with the AERK parameters.
            % Note, e=Magnus(tdc) and es=Magnus(tc)
            
            
            T_surface_level = 0;
            LW_surface_level = 0;
            T_dew_surface_level = 0;
            p_surface_level = 0;
            
            %loop over relevant coordinates
            for i=1:size(coords,1)
                T_surface_level = T_surface_level + idw(coords(i,1), coords(i,2)) .* squeeze(era.T2(coords(i,1), coords(i,2),:));
                LW_surface_level = LW_surface_level + idw(coords(i,1), coords(i,2)) .* squeeze(era.LW(coords(i,1), coords(i,2),:));
                T_dew_surface_level = T_dew_surface_level + idw(coords(i,1), coords(i,2)) .* squeeze(era.Td2(coords(i,1), coords(i,2),:));
                p_surface_level = p_surface_level + idw(coords(i,1), coords(i,2)) .* squeeze(era.ps(coords(i,1), coords(i,2),:));
            end
            
            vapor_pressure_surface_level = Magnus(K2C(T_dew_surface_level));
            q_surface_level = 0.622 .* vapor_pressure_surface_level ./ p_surface_level;
            
            Lin_pressure_level = compute_Lin(forcing, LW_surface_level, T_surface_level-273.15, temperature, q_surface_level, q, p_surface_level, p); 
            
            % Diagnose the all sky emissivity at grid.
            
%             all_sky_emissivity_surface_level = LW_surface_level./(sbc .* T_surface_level.^4);
%             
%             % Use the vapor pressure and temperature to calculate clear sky
%             % emssivity at grid and subgrid. [also function]
%             x1=0.43; x2=5.7;
%             %cef=0.23+x1.*(vpf./Tout).^(1/x2);
%             clear_sky_emissivity_surface_level = 0.23+x1.*(vapor_pressure_surface_level./T_surface_level).^(1/x2);
%             
%             cloud_emissivity_surface_level = all_sky_emissivity_surface_level - clear_sky_emissivity_surface_level;
            
            %to get clear-sky emissivity at the true altitude, cef, one needs q, p and
            %T  at the surface from previous timestep
            %     wf=q2w(qout); % Convert to mixing ratio at fine grid.
            %     vpf=wp2e(wf,pout);
            %     clear_sky_emissivity_surface_altitude = 0.23+x1.*(vpf./Tout).^(1/x2);
            %     % Use the former cloud emissivity to compute the all sky emissivity at subgrid.
            %     aef=clear_sky_emissivity_surface_altitude + all_sky_emissivity_surface_level - clear_sky_emissivity_surface_level;
            %     LWf=aef.*sbc.*Tout.^4;
            
            %     % Scale LW with terrain configuration, considering both occlusion by
            %     % and emissions from the surrounding terrain. From Dozier & Frew 1990.
            %     fout.LW(:,n)=svf.*LWf;%+0.5.*(1+cos(slp)).*(1-svf).*0.99.*5.67e-8.*(273.15.^4);
            
            
            %% Shortwave radiation
            
            
            % Get the solar geometry, which is assumed to be approximately invariant
            % for the entire area of interest. This is relatively unproblematic
            % even for a domain with an extent of 0.5 deg (500 km) or so.
            solar_azimuth = era.t'.*0;
            solar_zenith = era.t'.*0;
            for i = 1:size(era.t,2)
                [sazn, szen] = solargeom(forcing, era.t(1,i), target_lat, target_lon);
                solar_azimuth(i,1) = sazn;
                solar_zenith(i,1) = szen;
            end
            
            % Estimate the average solar zenith and azimuth over the last time step
            % (if available). This is important, otherwise the (average) shortwave
            % radiation for a timestep (which is what we are estimating) can be quite biased.
            solar_azimuth = [solar_azimuth(1,1); (solar_azimuth(1:end-1,1) + solar_azimuth(2:end,1))./2];
            solar_zenith = [solar_zenith(1,1); (solar_zenith(1:end-1,1) + solar_zenith(2:end,1))./2];
            
            
            % Compute downwelling TOA SW irradiance (i.e. the incoming shortwave
            % incident on a horizontal plane at the TOA), by accounting for the
            % solar zentih angle.
            mu0=max(cos(solar_zenith),0); % Trunacte negative values.
            % May also want to consider treating values mu0<0 for prominent topography
            % when the horizon  angles are less than 0.
            sunset=mu0<cosd(89);%(mu0==0); % Sunset switch.  sunset=(mu0==0); % Sunset switch.
            % Note, it would be better to use the true average ((1/tau) integral_t^(t+tau) mu0 dt)
            % but this approximation should be ok.
            SW_toa = S0.*mu0;
            
            % Get the surface shortwave and pressure from ERA5 and interpolate to fine grid.
            SW_surface_level = 0;
            p_surface_level = 0;
            %loop over relevant coordinates
            for i=1:size(coords,1)
                SW_surface_level = SW_surface_level + idw(coords(i,1), coords(i,2)) .* squeeze(era.SW(coords(i,1), coords(i,2),:));
                p_surface_level = p_surface_level + idw(coords(i,1), coords(i,2)) .* squeeze(era.ps(coords(i,1), coords(i,2),:));
            end
            %correct negative values
            SW_surface_level = max(SW_surface_level, 0);
            
            clearness_index = SW_surface_level./SW_toa;
            diffuse_fraction = 0.952 - 1.041 .* exp(-1.*exp(2.3-4.702.*clearness_index));
            diffuse_fraction = max(diffuse_fraction, 0);
            
            % Diffuse component.
            SW_surface_level_diffuse = diffuse_fraction .* SW_surface_level; % Diffuse shortwave radiation.
            SW_surface_level_direct = SW_surface_level - SW_surface_level_diffuse; %Direct shortwave radiation.
            
            SW_surface_level_diffuse = sky_view_factor .* SW_surface_level_diffuse; % Scale diffuse part with the sky-view factor. ERROR IN ORIGINAL SCRIPT, CHECK?
            
            % Note that even if it is physically correct, the above sky view scaling
            % can bias models (by removing a compensating error source) that
            % don't account for reflections from the surrounding terrain. -> corrected
            % in Juditha's forcing class, needs to be implemented!
            
            % Scale direct shortwave using Beer's law (see Aalstad 2019, Appendix A)
            SW_true_altitude_direct = SW_toa.*((SW_surface_level_direct./(max(SW_toa,1e-12))).^(p./p_surface_level)); % Combined but
            
            
            % Illumination angles.
            cosif = mu0.*cos(slope) + sin(solar_zenith).*sin(slope).*cos(solar_azimuth - aspect); % Cosine of illumination angle at fine (target) grid.
            selfshadow = cosif<0; % Self shadowing, occurs when |saz-asp|>90
            cosif(selfshadow)=0;
            cosic = mu0; % Cosine of illumination angle at coarse grid (ERA5 assumes slope=0 deg)
            
            % Binary shadow masks.
            % absdiff = abs(horizon_angle_bins - solar_azimuth); % To identify the relevant (minimum difference) horizon azimuth bin.
            % [~,binis] = min(absdiff); % Index of the relevant bin, returns first optimum index if there are multiple.
            % horizons = horizon_angle(:,binis); % Identify the relevant horizon angles for each target grid cell.
            
            horizons = solar_azimuth .*0;
            diff = solar_azimuth .*0 .*1e9;
            for i=1:size(horizon_angle_bins,1)
                update = abs(horizon_angle_bins(i,1) - solar_azimuth) < diff;
                diff(update) = abs(horizon_angle_bins(i,1) - solar_azimuth(update));
                horizons(update) = horizon_angle(i,1);
            end
            
            seln = max(pi/2 - solar_zenith, 0); % Solar elevation angle (in degrees for consistency with horizon angles).
            shade = horizons>seln; % True if a grid cell is in shadow.
            
            % Terrain corrected direct shortwave.
            SW_true_altitude_direct = SW_true_altitude_direct.*(cosif./cosic).*(1 - shade); % Terrain corrected direct shortwave at subgrid.
            
            SW_true_altitude = SW_true_altitude_direct + SW_surface_level_diffuse;
            SW_true_altitude(sunset) = 0;


            forcing.DATA.precipitation = precipitation;
            forcing.DATA.Tair = temperature;
            %forcing.DATA.sky_emissivity = cloud_emissivity_surface_level;
            forcing.DATA.Lin_pressure_level = Lin_pressure_level;
            forcing.DATA.Sin = SW_true_altitude;
            forcing.DATA.Lin_pressure_level = Lin_pressure_level;
            forcing.DATA.L_in = LW_surface_level;
            forcing.DATA.q = q;
            forcing.DATA.p = p;
            forcing.DATA.p_forcing_height = p_forcing_height;
            forcing.DATA.z = altitude - tile.PARA.altitude;
            forcing.DATA.wind = wind;
            forcing.DATA.timeForcing = (era.t)';
            
            
            forcing.DATA.wind(forcing.DATA.wind<0.5)=0.5; %set min wind speed to 0.5 m/sec to avoid breakdown of turbulence            
            
            
            if isempty(forcing.PARA.start_time) || isnan(forcing.PARA.start_time(1,1)) %|| ~ischar(forcing.PARA.start_time)
                forcing.PARA.start_time = forcing.DATA.timeForcing(1,1);
            else
                forcing.PARA.start_time = datenum(forcing.PARA.start_time(1,1), forcing.PARA.start_time(2,1), forcing.PARA.start_time(3,1));
            end
             
            if isempty(forcing.PARA.end_time) || isnan(forcing.PARA.end_time(1,1)) %~ischar(forcing.PARA.end_time)
                forcing.PARA.end_time = floor(forcing.DATA.timeForcing(end,1));
            else
                forcing.PARA.end_time = datenum(forcing.PARA.end_time(1,1), forcing.PARA.end_time(2,1),forcing.PARA.end_time(3,1));
            end
            
            %initialize TEMP
            forcing.TEMP.snowfall = 0;
            forcing.TEMP.rainfall = 0;
            forcing.TEMP.Tair = 0;
            forcing.TEMP.Lin = 0;
            forcing.TEMP.Sin = 0;
            forcing.TEMP.q = 0;
            forcing.TEMP.p = 0;
            forcing.TEMP.airT_height = 0;
            forcing.TEMP.wind_height = 0;
            forcing.TEMP.wind = 0;
            forcing.TEMP.timeForcing = 0;
            
            forcing.PARA.airT_height = forcing.PARA.min_height_pl_above_ground; % can be removed later     
            forcing.PARA.wind_height = forcing.PARA.min_height_pl_above_ground;
        end
        
        function forcing = interpolate_forcing(forcing, tile)
            t = tile.t;

            posit=floor((t-forcing.DATA.timeForcing(1,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1)))+1;
            
            forcing.TEMP.Sin = forcing.DATA.Sin(posit,1)+(forcing.DATA.Sin(posit+1,1)-forcing.DATA.Sin(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
            forcing.TEMP.Sin = double(forcing.TEMP.Sin);
            
            forcing.TEMP.Tair=forcing.DATA.Tair(posit,1)+(forcing.DATA.Tair(posit+1,1)-forcing.DATA.Tair(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
            forcing.TEMP.Tair = double(forcing.TEMP.Tair);
            
            forcing.TEMP.wind=forcing.DATA.wind(posit,1)+(forcing.DATA.wind(posit+1,1)-forcing.DATA.wind(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
            forcing.TEMP.wind = double(forcing.TEMP.wind);
            
            forcing.TEMP.q = forcing.DATA.q(posit,1)+(forcing.DATA.q(posit+1,1)-forcing.DATA.q(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
            forcing.TEMP.q = double(forcing.TEMP.q);
            
            forcing.TEMP.p = forcing.DATA.p(posit,1)+(forcing.DATA.p(posit+1,1)-forcing.DATA.p(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
            forcing.TEMP.p = double(forcing.TEMP.p);           
            
            forcing.TEMP.airT_height = forcing.DATA.z(posit,1)+(forcing.DATA.z(posit+1,1)-forcing.DATA.z(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
            forcing.TEMP.airT_height = double(forcing.TEMP.airT_height);
            
            forcing.TEMP.wind_height = forcing.TEMP.airT_height;
            
            precipitation = forcing.DATA.precipitation(posit,1)+(forcing.DATA.precipitation(posit+1,1)-forcing.DATA.precipitation(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
            
            forcing.TEMP.snowfall = precipitation .* forcing.PARA.snow_fraction .* (double(tile.TOP.NEXT.STATVAR.T2m <= forcing.PARA.all_snow_T)  + ...
                double(tile.TOP.NEXT.STATVAR.T2m > forcing.PARA.all_snow_T & tile.TOP.NEXT.STATVAR.T2m <= forcing.PARA.all_rain_T) .* ...
                (tile.TOP.NEXT.STATVAR.T2m - forcing.PARA.all_snow_T) ./ max(1e-12, (forcing.PARA.all_rain_T - forcing.PARA.all_snow_T)));
            forcing.TEMP.snowfall = double(forcing.TEMP.snowfall);
            
            forcing.TEMP.rainfall = precipitation .* forcing.PARA.rain_fraction .* (double(tile.TOP.NEXT.STATVAR.T2m >= forcing.PARA.all_rain_T)  + ...
                double(tile.TOP.NEXT.STATVAR.T2m > forcing.PARA.all_snow_T & tile.TOP.NEXT.STATVAR.T2m < forcing.PARA.all_rain_T) .* ...
                (1 - (tile.TOP.NEXT.STATVAR.T2m - forcing.PARA.all_snow_T) ./ max(1e-12, (forcing.PARA.all_rain_T - forcing.PARA.all_snow_T))));
            forcing.TEMP.rainfall = double(forcing.TEMP.rainfall);
            
            %sky_emissivity = forcing.DATA.sky_emissivity(posit,1)+(forcing.DATA.sky_emissivity(posit+1,1)-forcing.DATA.sky_emissivity(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
            
            Lin_pressure_level = forcing.DATA.Lin_pressure_level(posit,1)+(forcing.DATA.Lin_pressure_level(posit+1,1)-forcing.DATA.Lin_pressure_level(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
            p_forcing_height = forcing.DATA.p_forcing_height(posit,1)+(forcing.DATA.p_forcing_height(posit+1,1)-forcing.DATA.p_forcing_height(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
  
            Lin_surface = compute_Lin(forcing, Lin_pressure_level, forcing.TEMP.Tair, tile.TOP.NEXT.STATVAR.T2m, forcing.TEMP.q, tile.TOP.NEXT.STATVAR.q2m, p_forcing_height, forcing.TEMP.p);

            %disp([Lin_pressure_level Lin_surface])
            %disp([p_forcing_height forcing.TEMP.p])
%             eps0=0.622;
%             q2w=@(q) 0.5.*(1-sqrt(1-4.*q)); % Mixing ratio from specific humidity based on 2nd order Taylor series expansion.
%             wp2e=@(w,p) 0.5.*p.*(-1+sqrt(1+4.*w./eps0)); % Vapor pressure from mixing ratio based on 2nd order Taylor series expansion.
% 
%             %wf=q2w(tile.TOP.NEXT.STATVAR.q2m); % Convert to mixing ratio at fine grid.
%             wf=q2w((tile.TOP.NEXT.STATVAR.q2m + forcing.TEMP.q)./2);
%             %wf=q2w(forcing.TEMP.q);
%             vpf=wp2e(wf, forcing.TEMP.p);
%             
%             x1=0.43; x2=5.7;
%             %clear_sky_emissivity_surface_altitude = real(0.23+x1.*(vpf./(tile.TOP.NEXT.STATVAR.T2m + 273.15)).^(1/x2));
%             clear_sky_emissivity_surface_altitude = real(0.23+x1.*(vpf./((tile.TOP.NEXT.STATVAR.T2m + forcing.TEMP.Tair)./2 + 273.15)).^(1/x2));
%             aef = max(0, min(1, clear_sky_emissivity_surface_altitude + sky_emissivity));
%             %LW = aef.* forcing.CONST.sigma .* (tile.TOP.NEXT.STATVAR.T2m+273.15).^4;
%             LW = aef.* forcing.CONST.sigma .* ((tile.TOP.NEXT.STATVAR.T2m + forcing.TEMP.Tair)./2 + 273.15).^4;
                     
            forcing.TEMP.Lin = forcing.PARA.sky_view_factor .* Lin_surface + (1 - forcing.PARA.sky_view_factor) .* forcing.CONST.sigma .* (tile.TOP.NEXT.STATVAR.T2m+273.15).^4;
            forcing.TEMP.Lin = double(forcing.TEMP.Lin);

            forcing.TEMP.t = t;
        end
        

        
        
        
        function [solar_azimuth,solar_zenith] = solargeom(forcing, Time,Latitude,Longitude)
            %% [saz,szen]=solargeom(time,latitude,longitude)
            % Adopted from the Sandia National Labs PVL Toolbox ephemeris routine.
            % Inputs:
            %   time = Time stamp vector (matlab datenum format) assumed to be in UTC
            %   latitude = Latitude
            %   longitude = Longitude
            % Outputs:
            %   saz = Solar azimuth angle [radians, anticlockwise from south]
            %   szen = Solar zentih angle [radians].
            % Link to the original toolbox:
            % https://pvpmc.sandia.gov/applications/pv_lib-toolbox/
            % References:
            % Stein et al. (2012), doi:10.1109/PVSC.2012.6318225 [MATLAB version]
            % Holmgren et al. (2018), doi:10.21105/joss.00884 [Python version]
            TZone=0;
            Longitude=-Longitude;
            tv=datevec(Time);
            Year=tv(:,1);
            v0=zeros(size(Year)); v1=ones(size(Year));
            DayOfYear=Time-datenum([Year v1 v1 v0 v0 v0])+1;
            DecHours=(Time - floor(Time)) .* 24;
            RadtoDeg=180/pi;
            DegtoRad=pi/180;
            Abber = 20/3600;
            LatR = Latitude * DegtoRad;
            UnivDate = DayOfYear + floor((DecHours + TZone)/24);
            UnivHr = mod((DecHours + TZone), 24);
            Yr = Year-1900;
            YrBegin = 365 * Yr + floor((Yr-1)/4)-0.5;
            Ezero = YrBegin + UnivDate;
            T = Ezero / 36525;
            GMST0 = 6/24 +38/1440 + (45.836 + 8640184.542 * T + 0.0929 * T.^2)/86400;
            GMST0 = 360 * (GMST0 - floor(GMST0));
            GMSTi = mod(GMST0 + 360*(1.0027379093 * UnivHr / 24),360);
            LocAST = mod((360 + GMSTi - Longitude), 360);
            EpochDate = Ezero + UnivHr / 24;
            T1 = EpochDate / 36525;
            ObliquityR = DegtoRad * (23.452294 - 0.0130125 * T1 - 0.00000164 * T1.^2 ...
                + 0.000000503 * T1.^3);
            MlPerigee = 281.22083 + 0.0000470684 * EpochDate + 0.000453 * T1 .^ 2 + ...
                0.000003 * T1 .^ 3;
            MeanAnom = mod((358.47583 + 0.985600267 * EpochDate - 0.00015 * T1 .^ 2 - ...
                0.000003 * T1 .^ 3), 360);
            Eccen = 0.01675104 - 0.0000418 * T1 - 0.000000126 * T1 .^ 2;
            EccenAnom = MeanAnom;
            E=0;
            while max(abs(EccenAnom - E)) > 0.0001
                E = EccenAnom;
                EccenAnom = MeanAnom + RadtoDeg .* Eccen .* sin(DegtoRad .* E);
            end
            TrueAnom = 2 * mod(RadtoDeg * atan2(((1 + Eccen) ./ (1 - Eccen)).^ 0.5 .* tan(DegtoRad * EccenAnom / 2), 1), 360) ;
            EcLon = mod(MlPerigee + TrueAnom, 360) - Abber ;
            EcLonR = DegtoRad * EcLon;
            DecR = asin(sin(ObliquityR) .* sin(EcLonR));
            %Dec = RadtoDeg * DecR;
            RtAscen = RadtoDeg * atan2(cos(ObliquityR).*(sin(EcLonR)),cos(EcLonR));
            HrAngle = LocAST - RtAscen ;
            HrAngleR = DegtoRad .* HrAngle ;
            %HrAngle = HrAngle - (360 .* sign(HrAngle) .* (abs(HrAngle) > 180));
            SunAz = RadtoDeg .* atan2(-1 * sin(HrAngleR), cos(LatR) .* tan(DecR) - sin(LatR) .* cos(HrAngleR));
            SunAz = SunAz + (SunAz < 0) * 360; %shift from range of [-180,180] to [0,360]
            SunEl = asind(cos(LatR) .* cos(DecR) .* cos(HrAngleR) + sin(LatR) .* sin(DecR));
            
            % Convert solar azimuth angle from [N,E,S,W]=[0,90,180,270] to [180, 90, 0
            % -90], i.e. the same as the aspect and horizon angle system.
            solar_azimuth=deg2rad(SunAz);
            solar_azimuth=(5*pi/2)-solar_azimuth;
            solar_azimuth=solar_azimuth-2.*pi.*(solar_azimuth>2.*pi);
            solar_azimuth=solar_azimuth+pi/2;
            solar_azimuth=solar_azimuth-2.*pi.*(solar_azimuth>pi);
            
            % Calculate solar zenith angle from solar elevation angle
            SunEl=deg2rad(SunEl);
            solar_zenith=(pi/2)-SunEl;
        end

        function Lin_pressure_level = compute_Lin(forcing, L_in, T_1, T_2, q_1, q_2, p_1, p_2)
            %based on Kaye L. Brubaker  Dara Entekhabi: An Analytic Approach to Modeling Land?Atmosphere Interaction: 1. Construct and Equilibrium Behavior
            % Wilfried Brutsaert: On a derivable formula for long?wave radiation from clear skies
            
            A=0.75;
            m= 1/7;
            ps = 1.05e5; %reference pressure level at sea level
            g= forcing.CONST.g;
            sigma = forcing.CONST.sigma;
            
            rho_w = forcing.CONST.rho_w;
            
            inverted = p_1>p_2;
            T_top = double(~inverted) .* T_1 + double(inverted) .* T_2;
            T_bottom = double(inverted) .* T_1 + double(~inverted) .* T_2;
            q_top = double(~inverted) .* q_1 + double(inverted) .* q_2;
            q_bottom = double(inverted) .* q_1 + double(~inverted) .* q_2;
            p_top = double(~inverted) .* p_1 + double(inverted) .* p_2;
            p_bottom = double(inverted) .* p_1 + double(~inverted) .* p_2;
            
            increments = 20;
            delta_p = -(p_top-p_bottom)./increments;
            epsilon =0;
            LW = 0;
            p_bottom = min(p_bottom, 1.05e5);
            p_top = min(p_top, p_bottom);
            
            for fraction = 1/(increments-1):1/(increments-1):1
                T = fraction .* T_top + (1-fraction).* T_bottom;
                q = fraction .* q_top + (1-fraction).* q_bottom;
                p = fraction .* p_top + (1-fraction).* p_bottom;
                %a = 2./3 .* q .*ps ./ rho_w ./g .* ( (p_bottom./ps).^(3/2) - (p_top./ps).^(3/2)) .* 100;
                a = 2./3 .* q .*ps ./ rho_w ./g .* ( 1 - (p./ps).^(3/2)) .* 100;
                dEps_da = A.* m.* (a).^(m-1); %in per cm
                da_dp = q ./ rho_w ./g .* (p./ps).^(1/2) .* 100;
                epsilon = epsilon + delta_p .* dEps_da .* da_dp;
                LW = LW + sigma.*(T+273.15).^4 .* delta_p .* dEps_da .* da_dp;
            end
            
            L_in_1 = L_in .* (1-epsilon) + LW; %projecting downwards
            L_in_2 = (L_in - LW) ./ (1-epsilon); %projecting upwards

            Lin_pressure_level = double(~inverted) .* L_in_1 + double(inverted) .* L_in_2;
            Lin_pressure_level = real(double(Lin_pressure_level));
        end
 
                
    end
end