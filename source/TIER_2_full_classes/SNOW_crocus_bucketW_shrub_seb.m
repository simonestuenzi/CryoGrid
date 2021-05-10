%========================================================================
% CryoGrid GROUND class SNOW_crocus_bucketW_seb
% CROCUS snow model Vionnet et al., 2012, but with simpler layer splitting and regridding scheme compared to CROCUS 
% temperature and windspeed-dependent initial snow density, snow microstructure (dendricity, sphericity, grain size), 
% compaction, sublimation, water flow, refreezing, variable albedo.
% Meltwater exceeding the available pore space within the snow cover is automatically removed.
% R. Zweigel, S. Westermann, October 2020
%========================================================================

classdef SNOW_crocus_bucketW_shrub_seb < SNOW_crocus_bucketW_seb

    properties
        
    end
    
    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------

       
        function snow = compute_diagnostic(snow, tile)
            snow = compute_diagnostic@SNOW_crocus_bucketW_seb(snow, tile);
            snow.STATVAR.albedo = snow.STATVAR.albedo.*
        end
        
        function snow = compute_diagnostic_CHILD(snow, tile)
            
            snow = compute_diagnostic_CHILD@SNOW_crocus_bucketW_seb(snow, tile);
            snow.STATVAR.albedo = snow.STATVAR.albedo.* (1-snow.PARA.stem_area); 
            
        end
        
      
    
end
