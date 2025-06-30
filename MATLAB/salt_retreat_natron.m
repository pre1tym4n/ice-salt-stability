function [P_depth,dz_center] = salt_retreat_natron(wt_frac_salt,T,max_depth)

%% initializing parameters
%constants

k=8.6173e-5;        %Boltzmann's constant in electron volts (THIS ONE MUST BE IN eV FOR MCCORD ET AL EQS)
R=8.31;             % Jules/mol/kelvin
Av=6.02214076e23;   %Avagadro's number to convert number density to mol density
delta_H=51058;
Pt=611; 
Tt=273.16; 
rho_ice=925;        %density of water ice in kg/m^3

tau=2;     %tortuoisty 
phi=0.2;   %porosity
r=(10e-6)/2;   %particle radius (NOT diameter!!)
kb=1.380649e-23;   % J/K        %Boltzman's constant in MKS units (THIS ONE MUST BE IN J/K FOR SCHORGHOFER EQS)
m_molecular=2.9915e-26; %kg molecular mass of water

% time scaling
eyr=3.155e7;    %earth year in seconds


for zz=1:max(size(max_depth))
%array setup for depth and layer thickness

z0=0;            %initial dessicated overlayer 
dz_spacing=25e-6;
z=z0+(0:dz_spacing:max_depth(zz));

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

lithic_density=2000;    %kg/m^3 density of Ceres regolith if 2700, 2000 to match Tom's result
rho_bulk=lithic_density*(1-phi);  %density of overburden layer
molar_mass_natron=0.28614; %kg/mol of natron 

n0=(wt_frac_salt./molar_mass_natron).*Av.*rho_bulk*n_waters; %source flux of water molecules from regolith material of a certain density
n0_array=n0.*ones(1,max(size(dz_edges)));  %assign source flux to each layer

%calculate alpha for water vapor loss from salts directly

alpha=nu*exp(-Ea/(k*T));

%calculate saturation vapor pressure for this temperature
Psat=Pt.*exp((-(delta_H/R))*((1/T)-(1/Tt))); 

%Knudsen diffusion equation prefactor D

Dk=(pi/(8+pi))*(phi/(1-phi))*(r/tau)*sqrt(8*kb*T/(pi*m_molecular)); 


%% retreat rate and exposure time exposure time

J_bar=((4*pi)/(8+pi)).*(phi/(1-phi))*(1/tau)*sqrt((m_molecular)/(2*pi*kb*T)).*Psat.*(r./cumsum(dz_center)); 

Ri=J_bar./rho_ice; 

recession_time=flip(cumsum(flip(dz_center./Ri)));

%% Calculate loss rate based on exposure time 

tic

time_since_exposure=recession_time; 

n_remain_flux=alpha.*n0_array.*(exp(-alpha.*time_since_exposure)); %calculate source term for each layer at the elapsed time and save in a time by z array                    
toc

%% Calculate pressure due to salts being released at each depth 

n_density=zeros(1, max(size(dz_edges))); 

for n=1:max(size(dz_edges))
     for m=1:n
         n_density(n)=n_density(n)+dz_edges(m)*sum(dz_edges(m:end).*n_remain_flux(m:end)); %% NEED THE SOURCE AT THE CENTER, INTEGRATE AT TOP AND OBTTOM LAYERS
     end
end

n_density_div_Dk=n_density./Dk; 
n_density_div_Av=n_density_div_Dk./Av; 

%calculate P based on ideal gas law from n_density above

P_depth=(n_density_div_Av*R*T);       %%N is really number density which includes the volume, so this is fine 

%normalize P with depth to the water saturation vapor pressure

P_depth_normal=P_depth/Psat; 


%% save file 

name=convertCharsToStrings(['salt_dehydration_iteration0_' num2str(max(z)) 'depth' num2str(wt_frac_salt) 'wt_fraction_salt' num2str(T) 'K_temp.mat']); 


save(name, 'T', 'z', 'dz_center', 'dz_edges', 'Dk', 'Psat', 'n_density', 'Ri', 'n_remain_flux', 'P_depth', 'recession_time');  

end

end