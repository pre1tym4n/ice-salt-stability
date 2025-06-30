%Salt_retreat_exposure_time_iteration1.m

clc
clear all 

%% Pre-calculate salt pressure with depth without ice table
wt_frac_salt=0.10;    % weight fraction of salts in the material being exposed
depths=[0,25e-6,7.5e-5,0.0001,0.0004,0.002,0.003,0.005,0.007,0.01,0.02,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5]; %sets depths of interest for the model, Iteration 0 must be pre-run at these depths and ice contents
T=155;      %temperature of interest (pick from a variety, use an external thermal model for Ceres-realistic parameters)

for n=2:max(size(depths))
    n
    salt_retreat_natron(wt_frac_salt,T,depths(n));
end


%% Load no-ice-table salt retreat model results  
% Model results from iteration 0 (salt_retreat_natron.m) which gives the salt dehydration pressures with depths without the presence of the ice table

final_pressure(1)=0; %sets final surface pressure (assumes  surface in contact with vacuum) 

for n=2:max(size(depths))
    load(['salt_dehydration_iteration0_' num2str(depths(n)) 'depth' num2str(wt_frac_salt) 'wt_fraction_salt.mat'])
    final_pressure(n)=P_depth(end); 
end 

%% initializing parameters

%heat/vapor constants

k=8.6173e-5;        %Boltzmann's constant in electron volts (McCord+ equations all use eV!!)
R=8.31;             % Jules/mol/kelvin, universal gas constant
Av=6.02214076e23;   %Avagadro's number to convert number density to mol density
delta_H=51058;      %enthalpy of sublimation of water 
Pt=611;             %reference pressure for C-C equation (in Schorghofer, 2005)
Tt=273.16; 
rho_ice=925;        %density of water ice in kg/m^3  %Tom uses 930 kg/m^3 vs this (source of error??)

tau=2;     %tortuoisty 
phi=0.2;   %porosity
r=(10e-6)/2;   %particle radius (NOT diameter!!)
kb=1.380649e-23;   % J/K        %Boltzman's constant in MKS units (THIS ONE MUST BE IN J/K FOR SCHORGHOFER EQS)
m_molecular=2.9915e-26; %kg molecular mass of water

% time scaling
eyr=3.155e7;    %earth year in seconds

%array setup for depth and layer thickness

z0=0;            %initial dessicated overlayer 
dz_spacing=25e-6;
z=z0+(0:dz_spacing:0.1);

for n=1:max(size(z))-1
   z_midpoints(n)=0.5*(z(n)+z(n+1));
end

dz_edges=z-circshift(z,1);
dz_edges(1)=[]; 

dz_center=[diff(z_midpoints) z_midpoints(1)-z(1)]; 

% Salt properties, including volume fraction, doing natron to start with 

Ea=0.7;   %escape energy in eV
nu=2e12;  %escape attempt frequency in 1/s
n_waters=10; %number of waters per salt grain

lithic_density=2000;    %kg/m^3 density of Ceres regolith
rho_bulk=lithic_density*(1-phi);  %density of regolith lag layer
molar_mass_natron=0.28614; %kg/mol of natron 

n0=(wt_frac_salt./molar_mass_natron).*Av.*rho_bulk*n_waters; %source flux of water molecules from regolith material of a certain density
n0_array=n0.*ones(1,max(size(dz_center)));  %assign source flux to each layer

alpha=nu*exp(-Ea/(k*T)); %calculate alpha (parameter in McCord+) for water vapor loss from salts directly

Psat=Pt.*exp((-(delta_H/R))*((1/T)-(1/Tt))); %calculate saturation vapor pressure for this temperature

Dk=(pi/(8+pi))*(phi/(1-phi))*(r/tau)*sqrt(8*kb*T/(pi*m_molecular)); % Knudsen diffusion equation prefactor D

%% Iteration 0: calculate pressure curve with depth w/o any affect from salts to set initial pressure conditions 

% retreat rate and exposure time exposure time

J_bar=((4*pi)/(8+pi)).*(phi/(1-phi))*(1/tau)*sqrt((m_molecular)/(2*pi*kb*T)).*Psat.*(r./cumsum(dz_center)); 

Ri=J_bar./rho_ice; %ice table recession rate

recession_time=cumsum(dz_center./Ri);

time_since_exposure=recession_time; 

n_remain_flux=alpha.*n0_array.*(exp(-alpha.*time_since_exposure)); %calculate source term for each layer at the elapsed time and save in a time by z array                    
n_density=zeros(1, max(size(dz_edges))); 

for n=1:max(size(dz_edges))
     for m=1:n
         n_density(n)=n_density(n)+dz_edges(m)*sum(dz_edges(m:end).*n_remain_flux(m:end)); % need the source at the center, integrate at top and bottom of layers
     end
end

n_density_div_Dk=n_density./Dk; 
n_density_div_Av=n_density_div_Dk./Av; 

%calculate P based on ideal gas law from n_density above

P_depth=(n_density_div_Av*R*T);       %N is really number density which includes the volume, so this is fine 

%normalize P with depth to the water saturation vapor pressure

P_depth_normal=P_depth/Psat; 

%% Iteration 1: calculate the retreat rates with with input from salt pressure at each layer

% Cacluate partial pressure that occurs when you get to each layer, the interpolation from the imported density profile with depth

P_iteration0=interp1(depths,final_pressure,cumsum(dz_center));    %% use z_midpoints for midpoints of layers instead of layer edges

P_effective=Psat-P_iteration0; 

J_bar_effective=((4*pi)/(8+pi)).*(phi/(1-phi))*(1/tau)*sqrt((m_molecular)/(2*pi*kb*T)).*P_effective.*(r./cumsum(dz_center)); 

Ri_effective=J_bar_effective./rho_ice; 

recession_time_effective=cumsum(dz_center./Ri_effective);

%% calculate time difference

time_difference=recession_time_effective-time_since_exposure; 

time_difference=recession_time_effective-recession_time; 
time_ice_only=recession_time/eyr/1e6;   %recession rate in Myrs
time_salt_and_ice=recession_time_effective/eyr/1e6; %recession rate in Myrs
time_diff_Myr=time_difference/eyr/1e6;  %recession rate in Myrs

%% make plots

figure(1)
loglog(cumsum(dz_center), P_iteration0, 'b*');
hold on 
loglog(depths, final_pressure, 'r*-')
set(gca, 'FontSize', 15)
legend('interpolation', 'values from iteration 0')
xlabel('Depth (m)')
ylabel('Pressure (Pa)')

zs_plot=cumsum(dz_center); 

figure(2)
hold on 
plot(zs_plot, time_diff_Myr, 'r.-')
plot(zs_plot, time_ice_only, 'b--')
plot(zs_plot, time_salt_and_ice, 'g--')
xlabel('depth')
ylabel('time (Myr)')
legend('time difference', 'ice table retreat', 'iteration 1 modified retreat')
set(gca, 'FontSize', 15)