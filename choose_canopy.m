function ground = choose_canopy(ground, t) % Simone LAI

% t_1 = datenum(t);

dateComponents = datevec(t);
% years = dateComponents(:, 1); % Extract the years
months = dateComponents(:, 2); % Extract the months
snowCoveredMonths = (months == 9 | months == 10 | months == 11 | months == 12 | months == 1 | months == 2 | months == 3 | months == 4 | months == 5);
snowCoveredData = t(snowCoveredMonths);

    ground.STATVAR.canopy.pai = [];
    ground.STATVAR.canopy.lai = [];
    ground.STATVAR.canopy.dlai = [];
    ground.STATVAR.canopy.sumlai = [];
    ground.STATVAR.mlcanopyinst.zw = [];
    ground.STATVAR.mlcanopyinst.zs = [];
    ground.STATVAR.canopy.dpai = [];
    ground.STATVAR.mlcanopyinst.sumpai = [];
    ground.STATVAR.flux.albsoib = [];
    ground.STATVAR.flux.albsoid = [];

if snowCoveredData   
    % Perform actions specific to snow-covered months
%     disp('Current timestep is in a snow-covered month.');
    ground.STATVAR.canopy.pai = ground.STATVAR.canopy.pai_winter;
    ground.STATVAR.canopy.lai = ground.STATVAR.canopy.lai_winter;
    ground.STATVAR.canopy.dlai = ground.STATVAR.canopy.dlai_winter(:,:);
    ground.STATVAR.canopy.sumlai = ground.STATVAR.canopy.sumlai_winter(:,:);
    ground.STATVAR.mlcanopyinst.zw = ground.STATVAR.mlcanopyinst.zw_winter(:,:);
    ground.STATVAR.mlcanopyinst.zs = ground.STATVAR.mlcanopyinst.zs_winter(:,:);
    ground.STATVAR.canopy.dpai = ground.STATVAR.canopy.dpai_winter(:,:);
    ground.STATVAR.mlcanopyinst.sumpai = ground.STATVAR.mlcanopyinst.sumpai_winter(:,:);
    ground.STATVAR.flux.albsoib = ground.STATVAR.flux.albsoib_winter(:,:);
    ground.STATVAR.flux.albsoid = ground.STATVAR.flux.albsoid_winter(:,:);
else
    % Perform actions for other months
%     disp('Current timestep is not in a snow-covered month.');
    ground.STATVAR.canopy.pai = ground.STATVAR.canopy.pai_summer;
    ground.STATVAR.canopy.lai = ground.STATVAR.canopy.lai_summer;
    ground.STATVAR.canopy.dlai = ground.STATVAR.canopy.dlai_summer(:,:);
    ground.STATVAR.canopy.sumlai = ground.STATVAR.canopy.sumlai_summer(:,:);
    ground.STATVAR.mlcanopyinst.zw = ground.STATVAR.mlcanopyinst.zw_summer(:,:);
    ground.STATVAR.mlcanopyinst.zs = ground.STATVAR.mlcanopyinst.zs_summer(:,:);
    ground.STATVAR.canopy.dpai = ground.STATVAR.canopy.dpai_summer(:,:);
    ground.STATVAR.mlcanopyinst.sumpai = ground.STATVAR.mlcanopyinst.sumpai_summer(:,:);
    ground.STATVAR.flux.albsoib = ground.STATVAR.flux.albsoib_summer(:,:);
    ground.STATVAR.flux.albsoid = ground.STATVAR.flux.albsoid_summer(:,:);
end

% currentMonth = month(t_1); % Extract the months from the datetime array
% % currentYear  = year(t);   % Extract the years from the datetime array
% 
% % snowCoveredMonths = (months == 10 | months == 11 | months == 12 | months == 1 | months == 2 | months == 3 | months == 4);
% 
% isSnowCoveredMonth = ismember(currentMonth, [10, 11, 12, 1, 2, 3, 4]); % && ismember(currentYear, unique(years(snowCoveredMonths)));
% 
% if isSnowCoveredMonth
%     % Perform actions specific to snow-covered months
%     disp('Current timestep is in a snow-covered month.');
%     ground.STATVAR.canopy.pai = ground.STATVAR.canopy.pai_winter;
%     ground.STATVAR.canopy.lai = ground.STATVAR.canopy.lai_winter;
%     ground.STATVAR.canopy.dlai = ground.STATVAR.canopy.dlai_winter(:,:);
%     ground.STATVAR.canopy.sumlai = ground.STATVAR.canopy.sumlai_winter(:,:);
%     ground.STATVAR.mlcanopyinst.zw = ground.STATVAR.mlcanopyinst.zw_winter(:,:);
%     ground.STATVAR.mlcanopyinst.zs = ground.STATVAR.mlcanopyinst.zs_winter(:,:);
%     ground.STATVAR.canopy.dpai = ground.STATVAR.canopy.dpai_winter(:,:);
%     ground.STATVAR.mlcanopyinst.sumpai = ground.STATVAR.mlcanopyinst.sumpai_winter(:,:);
%     ground.STATVAR.flux.albsoib = ground.STATVAR.flux.albsoib_winter(:,:);
%     ground.STATVAR.flux.albsoid = ground.STATVAR.flux.albsoid_winter(:,:);
% else
%     % Perform actions for other months
%     disp('Current timestep is not in a snow-covered month.');
%     ground.STATVAR.canopy.pai = ground.STATVAR.canopy.pai_summer;
%     ground.STATVAR.canopy.lai = ground.STATVAR.canopy.lai_summer;
%     ground.STATVAR.canopy.dlai = ground.STATVAR.canopy.dlai_summer(:,:);
%     ground.STATVAR.canopy.sumlai = ground.STATVAR.canopy.sumlai_summer(:,:);
%     ground.STATVAR.mlcanopyinst.zw = ground.STATVAR.mlcanopyinst.zw_summer(:,:);
%     ground.STATVAR.mlcanopyinst.zs = ground.STATVAR.mlcanopyinst.zs_summer(:,:);
%     ground.STATVAR.canopy.dpai = ground.STATVAR.canopy.dpai_summer(:,:);
%     ground.STATVAR.mlcanopyinst.sumpai = ground.STATVAR.mlcanopyinst.sumpai_summer(:,:);
%     ground.STATVAR.flux.albsoib = ground.STATVAR.flux.albsoib_summer(:,:);
%     ground.STATVAR.flux.albsoid = ground.STATVAR.flux.albsoid_summer(:,:);
% end

% snowCoveredDataIndices = snowCoveredMonths & ismember(years, unique(years(snowCoveredMonths)));
% 
% snowCoveredTimestamps = timestamps(snowCoveredDataIndices);

% if t == snowCoveredTimestamps
%     ground.STATVAR.canopy.pai = ground.STATVAR.canopy.pai_winter;
%     ground.STATVAR.canopy.lai = ground.STATVAR.canopy.lai_winter;
%     ground.STATVAR.canopy.dlai = ground.STATVAR.canopy.dlai_winter(:,:);
%     ground.STATVAR.canopy.sumlai = ground.STATVAR.canopy.sumlai_winter(:,:);
%     ground.STATVAR.mlcanopyinst.zw = ground.STATVAR.mlcanopyinst.zw_winter(:,:);
%     ground.STATVAR.mlcanopyinst.zs = ground.STATVAR.mlcanopyinst.zs_winter(:,:);
%     ground.STATVAR.canopy.dpai = ground.STATVAR.canopy.dpai_winter(:,:);
%     ground.STATVAR.mlcanopyinst.sumpai = ground.STATVAR.mlcanopyinst.sumpai_winter(:,:);
%     ground.STATVAR.flux.albsoib = ground.STATVAR.flux.albsoib_winter(:,:);
%     ground.STATVAR.flux.albsoid = ground.STATVAR.flux.albsoid_winter(:,:);
% else % Summer canopy
%     ground.STATVAR.canopy.pai = ground.STATVAR.canopy.pai_summer;
%     ground.STATVAR.canopy.lai = ground.STATVAR.canopy.lai_summer;
%     ground.STATVAR.canopy.dlai = ground.STATVAR.canopy.dlai_summer(:,:);
%     ground.STATVAR.canopy.sumlai = ground.STATVAR.canopy.sumlai_summer(:,:);
%     ground.STATVAR.mlcanopyinst.zw = ground.STATVAR.mlcanopyinst.zw_summer(:,:);
%     ground.STATVAR.mlcanopyinst.zs = ground.STATVAR.mlcanopyinst.zs_summer(:,:);
%     ground.STATVAR.canopy.dpai = ground.STATVAR.canopy.dpai_summer(:,:);
%     ground.STATVAR.mlcanopyinst.sumpai = ground.STATVAR.mlcanopyinst.sumpai_summer(:,:);
%     ground.STATVAR.flux.albsoib = ground.STATVAR.flux.albsoib_summer(:,:);
%     ground.STATVAR.flux.albsoid = ground.STATVAR.flux.albsoid_summer(:,:);
% end

% if ground.STATVAR.canopy.lai_winter < (0.5*ground.STATVAR.canopy.lai_summer)+(0.5*0.1)
%     ground.STATVAR.pftcon.slatop = ground.STATVAR.pftcon.slatop_deciduous;
% end  

end

% 
% 
% if t>735517 && t<735749 || t>735882 && t<736114 || t>736247 && t<736480 || t>736613 && t<736845  || t>736978 && t<737210 || t>737343 && t<737575        
%     ground.STATVAR.canopy.pai = ground.STATVAR.canopy.pai_winter;
%     ground.STATVAR.canopy.lai = ground.STATVAR.canopy.lai_winter;
%     ground.STATVAR.canopy.dlai = ground.STATVAR.canopy.dlai_winter(:,:);
%     ground.STATVAR.canopy.sumlai = ground.STATVAR.canopy.sumlai_winter(:,:);
%     ground.STATVAR.mlcanopyinst.zw = ground.STATVAR.mlcanopyinst.zw_winter(:,:);
%     ground.STATVAR.mlcanopyinst.zs = ground.STATVAR.mlcanopyinst.zs_winter(:,:);
%     ground.STATVAR.canopy.dpai = ground.STATVAR.canopy.dpai_winter(:,:);
%     ground.STATVAR.mlcanopyinst.sumpai = ground.STATVAR.mlcanopyinst.sumpai_winter(:,:);
%     ground.STATVAR.flux.albsoib = ground.STATVAR.flux.albsoib_winter(:,:);
%     ground.STATVAR.flux.albsoid = ground.STATVAR.flux.albsoid_winter(:,:);
% else % Summer canopy
%     ground.STATVAR.canopy.pai = ground.STATVAR.canopy.pai_summer;
%     ground.STATVAR.canopy.lai = ground.STATVAR.canopy.lai_summer;
%     ground.STATVAR.canopy.dlai = ground.STATVAR.canopy.dlai_summer(:,:);
%     ground.STATVAR.canopy.sumlai = ground.STATVAR.canopy.sumlai_summer(:,:);
%     ground.STATVAR.mlcanopyinst.zw = ground.STATVAR.mlcanopyinst.zw_summer(:,:);
%     ground.STATVAR.mlcanopyinst.zs = ground.STATVAR.mlcanopyinst.zs_summer(:,:);
%     ground.STATVAR.canopy.dpai = ground.STATVAR.canopy.dpai_summer(:,:);
%     ground.STATVAR.mlcanopyinst.sumpai = ground.STATVAR.mlcanopyinst.sumpai_summer(:,:);
%     ground.STATVAR.flux.albsoib = ground.STATVAR.flux.albsoib_summer(:,:);
%     ground.STATVAR.flux.albsoid = ground.STATVAR.flux.albsoid_summer(:,:);
% end
% 
% if ground.STATVAR.canopy.lai_winter < (0.5*ground.STATVAR.canopy.lai_summer)+(0.5*0.1)
%     ground.STATVAR.pftcon.slatop = ground.STATVAR.pftcon.slatop_deciduous;
% end  
% 
% end
% 
