# Energy conservation during precipitation in numerical models 

**Deepak Waman, KIT Nov./Dec. 2025**

## Abstract
Sedimentation of precipitation transports not only hydrometeor mass but also sensible and latent heat vertically through the atmospheric column. 
During sedimentation, hydrometeors carry internal energy (heat content) and latent energy (associated with phase) with them. 
This vertical energy transport should be properly accounted for in numerical simulations such that the total energy of an atmospheric column changes only through energy removal by precipitation reaching the surface. 
This energy conservation needs to be accounted for in numerical models to account for the temperature change from hydrometeor sedimentation.
Here, in the 2-moment NWP set-up of the ICON model (https://gitlab.dkrz.de/icon/icon-nwp.git, kitDW_2mom_energy_conservation branch), we implemented a formulation that accounts for a temperature change from energy conservation during precipitation fallout. 
The implementation is tested by simulating an idealised case, called Weisman-Klemp 1982, with convection initialized with warm bubbles.

## Implementation
```
alv = Latent heat of vaporization (2.5008E6 J kg-1)
als = Latent heat of sublimation (2.8345E6 J kg-1)
rv  = Gas constant of moist air (461 J k-1 kg-1)
rd  = Gas constant of dry air (287 J k-1 kg-1)
cpv = Specific heat of dry air at constant pressure (1867.46 J k-1 kg-1)
cvd = Specific heat of dry air at constant volume (717.6 J k-1 kg-1)

clw = Specific heat of liquid water (4192.6641 J k-1 kg-1)
cvv = Specific heat of moist air at constant volume (1407.95 J k-1 kg-1)
tmelt = melting temperature of ice/snow (273.15 K)
ci = Specific heat of ice (2108 J k-1 kg-1)
lvc = alv - tmelt x (cpv - clw)
lsc = als - tmelt x (cpv - ci) 

Hydrometeors mass before sedimentation
qv_old = qv; qc_old = qc; qr_old = qr; qi_old = qi;
qs_old = qs; qg_old = qg; qh_old = qh;
tk_old = tk   (Temperature after microphysical changes in qx due to latent heat change)
```
There are three steps:
1. **Internal energy that all hydrometeors hold at a certain vertical level (Before sedimentation)**
   
   e_int = $\rho$ $\times$ dz $\times$ $\times$ (cv $\times$ tk_old $-$ qliq $\times$ lvc $-$ qice $\times$ lsc)   ... [J m-2]
   
   Where, cv = cvd $\times$ (1.0 $-$ qtot) + cvv $\times$ qv + clw $\times$ qliq + ci $\times$ qice            ... [J k-1 kg-1]

   qliq = qc + qr; qice = qi + qs + qg + qh; qtot = qv + qliq + qice;

2. **Internal energy due to rain and ice precipitation (After sedimentation, which updates qx)**
   
   rhodzbydt = $\rho$ $\times$ $\frac{dz}{dt}$
   
   mass_flux_r = (qr_old - qr) $\times$ rhodzbydt   ... [kg m-2 s-1]

   mass_flux_ci = (qi_old - qi)  $\times$ rhodzbydt

   mass_flux_s = (qs_old - qs)  $\times$ rhodzbydt

   mass_flux_g = (qg_old - qg) $\times$ rhodzbydt

   mass_flux_h = (qh_old - qh)  $\times$ rhodzbydt

   mass_flux_ice = mass_flux_r + mass_flux_ci + mass_flux_s +mass_flux_g + mass_flux_h
   
   Then, the energy flux in rain and total ice is:

   tk_old1 - tempk at level k

   tk_old2 - tempk at level k-1 (one level below k)
   
   eflx_rain = mass_flux_rain $\times$ (clw $\times$ tk_old1 $-$ cvd $\times$ tk_old2 $-$ lvc)   ... [J s-1 m-2 or W m-2]

   term1 == clw $\times$ tk_old is the enthalpy leaving level k
   
   term2 == cvd $\times$ tk_old is the enthalpy of air moving up (from level k-1) & replacing enthalpy at level k

   Terms 1 and 2 in combination define the rain enthalpy (latent heat terms here account for phase change)
   
   'lvc' is the energy that would be released if rain evaporates
   
   eflx_ice = mass_flux_ice $\times$ (ci $\times$ tk_old $-$ cvd $\times$ tk_old $-$ lsc)

   ! Update total mass flux
   mass_flux_ice = mass_flux_ice + mass_flux_rain

   eflx = dt $\times$ [eflx_rain + eflx_ice]

   qliq = qc + qr

   qice = qi + qs + qg + qh

   qtot = qliq + qice

3. **Change in internal energy after sedimentation:**
   
   e_int = e_int - eflx    ... [J m-2]


   cv = ( cvd $\times$ (1. - qtot) + cvv $\times$ qv + clw $\times$ qliq + ci $\times$ qice ) $\times$ rho $\times$ dz

   
   tk = $\frac{e_{int} + \rho \times dz \times (qliq \times lvc + qice \times lsc)}{cv}$





4. **Idealized simulation with Weismen-Klamp warm bubble, the temperature with and without Energy Conservation (EC):**
   <img width="1389" height="540" alt="image" src="https://github.com/user-attachments/assets/7819cac3-6eab-4774-9d15-42b550cfe1b9" />

