# Matlab - Ice-retreat model with hydrated salt

A Matlab version of the model described in Landis et al. (2025), “Role of natron in delaying retreat of buried ice tables on Ceres” is contained here. 

## salt_retreat_natron.m

Script to calculate the salt vapor pressure with depth within a set of subsurface layers assuming only salt vapor, not the presence of an ice table. This script represents step 1 to estimate vapor pressures with depth in the second script contained in this archive, described below.  

## Salt_reatreat_time_calculation_salt_plus_ice_table.m

Script to calculate the salt and water ice table model described in the paper text. This script uses “salt_retreat_natron.m” for an initialization stage, and the user defines the depths in the subsurface at which the pressures and ice table retreats will be calculated. 
