pro f_vth_ch,engs,spec,tmk=tmk,volem=volem,engmin=engmin,engmax=engmax,engres=engres

  ; To try and get a f_vth like isothermal out of CHIANTI
  ;
  ;
  ; 22-Apr-2020 IGH
  ;
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ;defsysv,'!abund_file',!xuvtop+'/abundance/sun_coronal_1992_feldman_ext.abund';,/test
  ;defsysv,'!ioneq_file',!xuvtop+'/ioneq/chianti.ioneq';,/test

  if (n_elements(engmin) ne 1) then engmin=0.5
  if (n_elements(engmax) ne 1) then engmax=50
  if (n_elements(engres) ne 1) then engres=0.1
  if (n_elements(tmk) ne 1) then tmk=1
  if (n_elements(volem) ne 1) then volem=1d47

  hc=1.98644568e-15
  kevinj=1.602176565e-16

  wavemin=hc/engmax/kevinj
  wavemax=hc/engmin/kevinj

  waveres=(wavemax-wavemin)/((engmax-engmin)/engres)
  isotemp=tmk*1d6
  discarea=3.14*6.955e10^2
  ; Giving it cm^-5 not cm^-3 => output spectrum will be in cm^-2
  colem=volem/discarea

  ch_synthetic,wavemin,wavemax, output=trans , pressure=1e15, $
    ioneq_name=!xuvtop+'/ioneq/chianti.ioneq', /photons,$
    logt_isothermal=alog10(isotemp), $
    logem_isothermal=alog10(colem),/all

  make_chianti_spec, trans, engs, ss,/cont, $
    bin_size=engres,  wrange=[engmin,engmax],$
    abund_name=!xuvtop+'/abundance/sun_coronal_1992_feldman.abund',/photons,/kev

  ;   To get rid of the final 1/sr
  spec=ss.spectrum*6.87e-5

end