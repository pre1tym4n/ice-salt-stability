function salty_ice_loss1, zf, wsalt=wsalt, pp_profile=pp_profile, verbose=verbose
;
; Purpose: Model ice retreat assuming the partial pressure from 
;          water vapor sourced from salt is known.
;          
; Inputs:
;    zf         - Ice table depth (m)
;    wsalt      - Weight fraction of salt (g/g) - default is 0.05 g/g
;    pp_profile - Structure containing the assumed depth profile of
;                 partial pressure (Pa), e.g. to ignore partial pressure
;                 effects, use
;                 
;                 pp_profile={z:[0d,1d],pp:[0d,0d]}
;                 
;                 This is the default if the keyword is not set.
;    verbose    - If set, then generate plots and diagnostics
;                 
; Outputs:
;    See output structure. The most important parameters are:
;    .Pp        - Estimated partial pressure contribution from dehydration 
;                 of salt at the final ice table depth (Pa)
;    .tt        - Vector of recession times given the input partial pressure profile (s).
;                 Note that tt[-1] is the time required to recede to a depth of zf.
; 
; THP - 11-Apr-2024 
;
compile_opt idl2

pp_prof={z:[0d,1d],pp:[0d,0d]}
if keyword_set(pp_profile) then pp_prof=pp_profile 

; Regolith parameters
;zf=0.15d0        ; Ice table depth (m)
dz=0.000025       ; Step size (m)
ncl=50L           ; Number of steps close to the ice table where vapor is being released from salt
temp=155d0        ; 155K temperature
porosity=0.2      ; porosity
diameter=10       ; grain size (um)
rho_grain=2000d0  ; grain density (kg/m3)
rho_bulk=rho_grain*(1-porosity)

; Salt parameters
wfrac=0.05           ; Weight fraction
if keyword_set(wsalt) then wfrac=wsalt
print, wfrac
form='Na2CO3(H2O)10' ; Emperical formula
nwat=10              ; Number of bound water modecules
Ea=0.7               ; Activation energy in eV
nu=2d12              ; Lattice frequency (s-1)

; Constants
spjyr=31557600d0      ; seconds per Julian year
nA=6.02214076d23      ; Avogadro's number (#/mol)
kB=1.380649d-23       ; Boltzman's constant (J/K)
eV=1.602176634e-19    ; J/eV
kB_eV=kB/eV           ; Boltzman's constant (eV/K)
Rg=8.31446261815324d0 ; Ideal gas constant (m3-Pa/mol/K)

; Calculate the volumetric number density of bound water
pp=parse_chemical_formula(form)
n0=wfrac/(pp.mw/1000)*nA*rho_bulk*nwat       ; water molecules per m3

; Ice table retreat sans salt
tic
nn=long(zf/dz)+1  ; number of steps
zz=dindgen(nn)*dz
vv=dblarr(nn)
ppp=interpol(pp_prof.pp,pp_prof.z,zz)
for i=1L,nn-1 do begin & $
  r=ice_loss_rate(temperature=temp,depth=zz[i],pp=ppp[i],porosity=porosity,diameter=diameter) & $
  vv[i]=r.ri & $
endfor
tt=dblarr(nn)
for i=1L,nn-1 do tt[i]=tt[i-1]+dz/vv[i]

; Diffusion coefficent and saturation vapor pressure given by last call to ice_loss_rate
dk=r.dk  ; m2/s
Ps=r.ps  ; Pa

; exposure time for layers
zmid=(zz[1:-1]+zz[0:-2])/2
delz=zz[1:-1]-zz[0:-2]
texp=interpol(tt[-1]-tt,zz,zmid)

; Calculate the instantaneous vapor release from salt
; at the last time step
alpha=nu*exp(-Ea/(kB_eV*temp))
source=alpha*n0*exp(-alpha*texp)  

; Number density with depth
n=dblarr(source.length)   ; number of molecules per m3
if ncl eq 0 then begin
  i0=0
endif else begin
  i0=source.length-1-ncl & if i0 lt 0 then i0=0
endelse
for i=i0,source.length-1 do begin
  for j=0L,i do begin
      n[i] += delz[j]*total(delz[j:-1]*source[j:-1])
  endfor
endfor
n=n/dk
toc

; Calculate pressure with depth
Psalt=(n/nA)*Rg*temp

if keyword_set(verbose) then begin

  print, 'Partial pressure at the ice table (fraction of the saturation pressure) = ',Psalt[-1]/Ps

  ; Source as a function of depth
  w=window(dimensions=[800,800])
  p=objarr(10)
  p[0]=plot(zmid*100,source/nA,/current)
  p[0].symbol='o'
  p[0].sym_filled=1
  p[0].sym_size=0.5
  p[0].xtickdir=1 & p[0].xticklen=0.02
  p[0].xthick=2
  p[0].ytickdir=1 & p[0].yticklen=0.02
  p[0].ythick=2
  p[0].xtitle='Depth (cm)'
  p[0].ytitle='Source (mol/m$^{3}$/s)'
  p[0].xtickfont_name='Arial' & p[0].xtickfont_size=24
  p[0].ytickfont_name='Arial' & p[0].ytickfont_size=24
  p[0].position=[0.255,0.13,0.95,0.95]

  ; Pressure as a function of depth
  w=window(dimensions=[800,800])
  p=objarr(10)
  p[0]=plot(zmid*100,Psalt/ps,/current)
  p[0].symbol='o'
  p[0].sym_filled=1
  p[0].sym_size=0.5
  p[0].xtickdir=1 & p[0].xticklen=0.02
  p[0].xthick=2
  p[0].ytickdir=1 & p[0].yticklen=0.02
  p[0].ythick=2
  p[0].xtitle='Depth (cm)'
  p[0].ytitle='Pressure from salt normalized to P$_{S}$ (unitless)'
  p[0].xtickfont_name='Arial' & p[0].xtickfont_size=24
  p[0].ytickfont_name='Arial' & p[0].ytickfont_size=24
  p[0].position=[0.21,0.13,0.95,0.95]

endif

return, {temp:temp, porosity:porosity, diameter:diameter, wfrac:wfrac, form:form, nwat:nwat, $
         zz:zz, tt:tt, Pp:Psalt[-1], Ps:Ps, Dk:dk, zmid:zmid, alpha:alpha, texp:texp, source:source, n:n, Psalt:Psalt}
end