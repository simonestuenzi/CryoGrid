%made for snow extent DA in ESA CCI, runs for multiple grid cells simultaneously 

%good to the point that PBS produces weights , bst-fitting ensemble member
%can be identified and written in outoput file

classdef OBS_SCFG_ESA_CCI < matlab.mixin.Copyable

    
    properties
        PARA
        CONST
        STATVAR
        TEMP
    end
    
    methods
        function obs  = provide_PARA(obs )
            
            obs.PARA.obs_filename = [];
            obs.PARA.obs_folder = [];
            obs.PARA.std_observations = [];
        end
        
        function obs = provide_CONST(obs)
            
        end
        
        function obs = provide_STATVAR(obs)
            
        end
        
        function obs = finalize_init(obs, tile)
            if exist([obs.PARA.obs_folder obs.PARA.obs_filename '_' datestr(tile.t,'yyyy') '.nc'])==2
                obs.STATVAR.observations = double(ncread([obs.PARA.obs_folder obs.PARA.obs_filename '_' datestr(tile.t,'yyyy') '.nc'], 'scfg', [tile.PARA.range(1) 1], [tile.PARA.range(end)-tile.PARA.range(1)+1 Inf], [1 1]));
                obs.STATVAR.observations = obs.STATVAR.observations';
                obs.STATVAR.observations(obs.STATVAR.observations==200) = NaN;
                obs.STATVAR.observations = obs.STATVAR.observations ./ 100;
                obs.STATVAR.time = [datenum(year(tile.t),9,1)+0.5:datenum(year(tile.t),9,1)+size(obs.STATVAR.observations,1)-0.5]';
                obs.STATVAR.obs_variance = zeros(size(obs.STATVAR.observations,1), size(obs.STATVAR.observations,2)) + obs.PARA.std_observations.^2;
            else
                obs.STATVAR.observations = [];
                obs.STATVAR.time = [];
                obs.STATVAR.obs_variance = [];
            end
        end

    end
end

