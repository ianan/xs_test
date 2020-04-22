pro fvth_vs_chianti

  ; Comparison of what f_vth gives vs CHIANTI's synthetic spectrum for a chosen T, EM
  ;
  ; f_vth is from:
  ;   https://hesperia.gsfc.nasa.gov/ssw/packages/xray/idl/f_vth.pro
  ; which is calling chianti_kev.pro
  ;   https://hesperia.gsfc.nasa.gov/ssw/packages/xray/idl/chianti_kev.pro
  ; which loads in a previously calculated database from CHIANTI which is in
  ;   sswdb_xray_chianti = concat_dir('$SSWDB_XRAY','chianti') => ssw/packages/xray/dbase/chianti/
  ; and for default of line+cont we have
  ;   getenv('CHIANTI_LINES_FILE') => chianti_lines_1_10_v71.sav
  ; or can also get this via
  ;   print,chianti_kev_version()
  ; default abundances used are Feldman coronal from where given as absolute value in alog10()
  ;   ssw/packages/xray/dbase/chianti/sun_coronal_abund.txt
  ; which is the same as http://articles.adsabs.harvard.edu/pdf/1992ApJS...81..387F and
  ;   ssw/packages/chianti/dbase/abundance/sun_coronal_1992_feldman.abund
  ; which can also be seen in xr_rd_abundance() where given as actual ratio relative to H
  ; 
  ; CHIANTI: Calculating it from scratch doing it two different ways:
  ;   * ch_synthetic + make_chianti_spec
  ;   * isothermal
  ; In both cases doing the full lines + continuum (freebound, freefree, two_photon) calc
  ;
  ;
  ; 22-Apr-2020 IGH
  ;
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ; Setup the model parameters
  mk2kev=0.08617
  em=1d47
  tmk=10

  ; Energy bin setup
  mine=2.5
  maxe=25
  ebin=0.04
  nengs=floor((maxe-mine)/ebin)

  edges=dblarr(2,nengs)
  edges[0,*]=mine+findgen(nengs)*ebin
  edges[1,*]=mine+(1.+findgen(nengs))*ebin
  de=edges[1,0]-edges[0,0]
  engs=get_edges(edges,/mean)
  
  ; Get the f_vth model spectrum, photons/s/keV/cm^2

  modf=f_vth( edges, [em*1d-49, tmk*mk2kev, 1.] )

  ; Now CHIANTI - my wrapper onto ch_synthetic.pro and make_chianti_spec.pro
  ; This takes a while to run as doing the full CHIANTI thing
  ; Work things out at the energy bin mid-points (though might be slightly offset.....)
  f_vth_ch,engsc,spec,tmk=tmk,volem=em,engmin=mine+0.5*ebin,engmax=maxe-0.5*ebin,engres=ebin

  ; Now try with CHIANTI's isothermal.pro procedure
  hc=1.98644568e-15
  kevinj=1.602176565e-16
  wavemin=hc/max(engs)/kevinj
  wavemax=hc/min(engs)/kevinj
  ; This returns the spectrum in photons/s/cm^2/sr if em in is in cm^-5 (1/area of Sun visible surface)
  isothermal, wavemin,wavemax,(wavemax-wavemin)/n_elements(engs),tmk*1d6,lmbda,speciso,$
    em=em/(3.14*6.955e10^2),/photons,/cont,pressure=1e15,$
    abund_name=!xuvtop+'/abundance/sun_coronal_1992_feldman.abund',$
    ioneq_name=!xuvtop+'/ioneq/chianti.ioneq'

  ; Might be some oddities with the lines as fixed bin width in lambda not energy
  dle=get_edges(hc/lmbda/kevinj,/width)
  ; To go from photons/s/cm^2/sr to photons/s/cm^2/keV
  speci=speciso*6.87e-5/dle
  speci=speci[0:nengs-1]
  engsi=hc/lmbda[0:nengs-1]/kevinj
  ; Just put them in to ascending energy order, to avoid confusion
  speci=reverse(speci)
  engsi=reverse(engsi)

  ;; Make a nice plot of everything
  @post_outset
  tlc_igh

  fid='ex1'

  figname='fvth_chianti_'+fid+'.eps'
  set_plot,'ps'
  device, /encapsulated, /color, /isolatin1,/inches, $
    bits=8, xsize=5, ysize=8,file=figname

  !p.multi=[0,1,3]
  !p.thick=4
  !p.charsize=2
  yr=[5d-4,5d4]
  xr=[2,16]
  xtit='Energy [keV]'
  ytit='photons/s/keV/cm!U2!N'
  
  plot,engs,modf,xtit=xtit,ytit=ytit,yrange=yr,xrange=xr,/ylog,/nodata
  oplot,engs,modf,color=8
  xyouts,0.9*xr[1],0.1*yr[1],'f_vth',color=8,/data,align=1,chars=1
  xyouts,0.9*xr[1],0.01*yr[1],string(tmk,format='(f4.1)')+' MK, 10^'+$
    string(alog10(em),format='(f4.1)')+' cm!U-3!N',color=0,/data,align=1,chars=1
  plot,engsc,spec,xtit=xtit,ytit=ytit,yrange=yr,xrange=xr,/ylog,/nodata
  oplot,engsc,spec,color=5
  xyouts,0.9*xr[1],0.1*yr[1],'ch_synthetic + make_chianti_spec',color=5,/data,align=1,chars=1
  plot,engsi,speci,xtit=xtit,ytit=ytit,yrange=yr,xrange=xr,/ylog,/nodata
  oplot,engsi,speci,color=4
  xyouts,0.9*xr[1],0.1*yr[1],'isothermal',color=4,/data,align=1,chars=1

  device,/close
  set_plot, mydevice
  convert_eps2pdf,figname,/del
  

  stop
end