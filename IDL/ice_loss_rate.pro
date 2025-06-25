function ice_loss_rate, temperature=temperature, $
                        depth=depth,             $
                        porosity=porosity,       $
                        tortuosity=tortuosity,   $
                        diameter=diameter,       $ 
                        Pp=Pp  
;
; Purpose: Calculate the instantaneous ice retreat rate
;          for a porous surface based on Shorghofer (2008).
;          
; Inputs: (keywords, defaults in parens)
;       temperature - temperature of the lag deposit  - Kelvin (293K)         
;       depth       - depth of the ice table in meters (1 mm)
;       porosity    - porosity of the regolith (0.5)
;       tortuosity  - tortuosiy of the pore network (2)
;       diameter    - grain size or pore diameter in meters (100 microns)
;       Pp          - partial pressure of water vapor at the ice surface
;                     from sources other than ice (e.g. hydrated salt). 
;                     The mass flux is reduced when Pp <> 0 and Pp < Ps. The 
;                     default is Pp = 0 (only sublimation of ice contributes
;                     to the vapor pressure at the surface of the ice 
;                     table).
;                     NOTE: If Pp > Ps, then the returned mass flux will
;                     be negative, indicating a growing ice table.
;                     
; Outputs: (structure)
;       .Ps         - equilibrium vapor pressure (Pa)
;       .Pp         - partial vapor pressure from other water vapor sources (Pa)
;       .vbar       - mean thermal velocity (m/s)
;       .Dk         - diffusion coefficient (m2/s)
;       .MFP        - mean free path of molecules (m) == 2 Dk / vbar  ***
;       .DT         - mean time between collisions (s) == .MFP / vbar
;       .Ltrk       - transport mean free path 
;       .Jbar       - mass flux of water from the surface (kg/m2/s)
;       .Ri         - instantaneous ice loss rate (m/s) (recession velocity)
; 
; *** From charts by Xiaolei Chen (diffusion_equation.pdf)
; 
; Example:  From Schorghofer (2008) - particle diameter of 100 microns,
;           tortuosity of 2, porosity of 0.5, depth of 1 m, and a 
;           temperature of 150K should give a retreat rate of 9 m per Gy.
;           
;           r=ice_loss_rate(temperature=150.,depth=1)
;           print, r.ri*1.e9*3.154e+7   ; note: 3.154e+7 s/year
;           ; this gives 8.5 m (per billion years)
; 
; Version:
;    Subroutine written by THPrettyman in 2017
;    1-May-2024 Updated by THPrettyman to include Pp to reduce ice table retreat rate.

compile_opt idl3

; defaults
T=293.d0  ; 20degC (about room temperature)
deltaZ=0.001 ; meters (default is 1 mm)
phi=0.5   ; porosity
tau=2d0   ; tortuosity 
Ppart=0d  ; partial vapor pressure from sources other than sublimation of ice
rr=100./2.d0*1.e-6 ; particle radius in meters (100 micron diameter particle default)
if keyword_set(temperature) then T=temperature
if keyword_set(depth) then deltaZ=depth
if keyword_set(porosity) then phi=porosity
if keyword_set(tortuosity) then tau=tortuosity
if keyword_set(diameter) then rr=diameter/2.d0*1.e-6 ; diameter is specified in microns
if keyword_set(Pp) then Ppart=Pp

; parameters
deltaH=51.058 ; water ice sublimation enthalphy, assumed constant (J/kg)
Pt=611. ; Triple point pressure (Pa)
Tt=273.16 ; Triple point temperature (K)
R=8.314459848 ; Universal gas constant (MJ/Mmol-K)
kB=1.3806485279d-23 ; Boltzmann constant (J/K)
;NA=6.02214085774d23 ; Avagadro's number (/mol)
; mass of a water molecule
;   c=parse_chemical_formula('H2O')
;   m=c.MW/NA/1000.d0 ; mass of a molecule in kg
m=2.9914699d-26 ; kg (from the above two lines of code)
rho_bulk_ice=930.d0 ; kg/m3

; equilibrium vapor pressure (Pa)
Ps=Pt*exp(-1000.*deltaH/R*(1./T-1./Tt)) ; factor of 1000 fits http://www.kayelaby.npl.co.uk/chemistry/3_4/3_4_1.html
if Ppart gt Ps then print, 'ice_loss_rate: Warning - the input partial pressure exceeds the saturation vapor pressure.'

; mean thermal velocity
vbar=sqrt(8.d0*kB*T/(!dpi*m))  ; m/s

; diffusion coefficient (m2/s)
Dk=!dpi/(8.d0+!dpi)*phi/(1.d0-phi)*vbar*rr/tau

; mass flux of water from the surface (kg/m2/s)
Jbar=2.d0*!dpi/(8.d0+!dpi)*phi/(1.d0-phi)*(1./tau)*sqrt(2.d0*m/(!dpi*kB*T))*(Ps-Ppart)*rr/deltaZ

; instantaneous loss rate (m/s)
Ri=Jbar/rho_bulk_ice

return, {Ps:Ps, Pp:Ppart, vbar:vbar, Dk:Dk, MFP:2.d0*Dk/vbar, DT:2.d0*Dk/vbar^2, Jbar:Jbar, Ri:Ri, temperature:T}
end