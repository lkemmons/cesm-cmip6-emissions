; SO2 and sulfate aerosol anthro emissions
; 2.5% of SO2 emitted as SO4.
; energy+industry emissions emitted at 100-300m
; MAM so4_a1, so4_a2 - different sectors in different modes
;  a1_surf: agricult+waste+solvents
;  a1_surf: shipping
;  a1_vert: energy+industrial
;  a2_surf: residential+transport
;  and num
; Reads concatenated file regridded to 1 deg 
; anthro files 1750-2014 (repeat 2014 for 2015)
;
pro make_emis_cmip6_so2so4_anthro

today = bin_date(systime())
todaystr = String(today[0:2],format='(i4,"/",i2.2,"/",i2.2)')
sdate_today = String(today[0:2],format='(i4,i2.2,i2.2)')
creation_date = sdate_today
newdir_date = '20170608'
so2filedate = '20170616'
thisfile = Routine_filepath()

path_a = '/glade/scratch/emmons/Emissions_CMIP6/concat/'
path_out = '/glade/p/cesm/chwg_dev/emmons/CMIP6_emissions_1750_2015_v'+newdir_date+'/'

;make new time array - gregorian calendar, days since 1750
ntim = (2015-1750+1)*12
time = fltarr(ntim)
date = lonarr(ntim)
for yr = 1750,2015 do begin
  for mon = 1,12 do begin
    itim =  (mon-1) + (yr-1750)*12
    date[itim] = yr*10000L + mon*100L + 16
    time[itim] = julday(mon,16,yr) - julday(1,1,1750)
    ;print,itim,mon,yr,date[itim],time[itim]
  endfor
endfor
  
mw_so2 = 64.
mw_so4_bam = 96.
mw_so4_mam = 115.
mw_num = 1.

; read SO2-anthro-sectors
 file_a = path_a+'emissions-cmip6_SO2_anthro-sectors_surface_175001-201412_0.9x1.25_c'+so2filedate+'.nc'
    print,'reading ',file_a
    ncid_a = ncdf_open(file_a)
    ncdf_varget,ncid_a,'lon',lon
    ncdf_varget,ncid_a,'lat',lat
    ncdf_varget,ncid_a,'date',date_a
    ncdf_varget,ncid_a,'emiss_agriculture',agriculture
    ncdf_varget,ncid_a,'emiss_energy',energy
    ncdf_varget,ncid_a,'emiss_industry',industry
    ncdf_varget,ncid_a,'emiss_transport',transport
    ncdf_varget,ncid_a,'emiss_resident',resident
    ncdf_varget,ncid_a,'emiss_solvents',solvents
    ncdf_varget,ncid_a,'emiss_waste',waste
    ncdf_varget,ncid_a,'emiss_shipping',shipping
    ncdf_attget,ncid_a,/global,'source_file',source_file_b
    source_file = String(source_file_b)
    ncdf_close,ncid_a
  ntima = n_elements(date_a)
nlon = n_elements(lon)
nlat = n_elements(lat)

  ;extend anthro, repeating 2014 for 2015
  ; grouping sectors, original SO2
  ag_sol_was = fltarr(nlon,nlat,ntim)
  ship = fltarr(nlon,nlat,ntim)
  res_tran = fltarr(nlon,nlat,ntim)
  ene_ind = fltarr(nlon,nlat,ntim)
  ag_sol_was[*,*,0:ntima-1] = agriculture+solvents+waste
  ag_sol_was[*,*,ntima:ntim-1] = ag_sol_was[*,*,ntima-12:ntima-1]
  ship[*,*,0:ntima-1] = shipping
  ship[*,*,ntima:ntim-1] = ship[*,*,ntima-12:ntima-1]
  res_tran[*,*,0:ntima-1] = resident+transport
  res_tran[*,*,ntima:ntim-1] = res_tran[*,*,ntima-12:ntima-1]
  ene_ind[*,*,0:ntima-1] = energy+industry
  ene_ind[*,*,ntima:ntim-1] = ene_ind[*,*,ntima-12:ntima-1]
  ag_sol_was_so2 = ag_sol_was*0.975
  ag_sol_was_so4 = ag_sol_was*0.025
  ship_so2 = ship*0.975
  ship_so4 = ship*0.025
  res_tran_so2 = res_tran*0.975
  res_tran_so4 = res_tran*0.025
  ene_ind_so2_surf = 0.975*ene_ind
  ene_ind_so4_surf = 0.025*ene_ind

  ; write so4_a1: sulfate accumulation mode for energy+industrial
  ; vertically distributed: 100-300m
  nalt = 8
  dz = 0.05 ;km
  altitude = findgen(nalt) * dz + 0.5*dz 
  altitude_int = findgen(nalt+1) * dz
  ene_ind_so2_vert = fltarr(nlon,nlat,nalt,ntim)
  ene_ind_so4_vert = fltarr(nlon,nlat,nalt,ntim)
  dalt = 2.e4 ;cm = 200m
  for ilon=0,nlon-1 do begin
   for ilat=0,nlat-1 do begin
     emis1 = reform(ene_ind[ilon,ilat,*])
     if (max(emis1) gt 0) then begin
       emis_vol = emis1 /dalt  ;(molecules/cm2/s)/cm
       for ialt = 3,6 do ene_ind_so2_vert[ilon,ilat,ialt,*] = 0.975*emis_vol
       for ialt = 3,6 do ene_ind_so4_vert[ilon,ilat,ialt,*] = 0.025*emis_vol
     endif
   endfor
  endfor
;----------
; calculate number emissions
  ;varname = 'emiss_ag_sol_was'
  diam = 0.134e-6
  rho = 1770.
  mw = 115.
  mass_particle = rho *(!PI/6.) *(diam)^3  ;mass per particle (kg/particle)
  num_ag_sol_was_so4 = ag_sol_was_so4 *mw /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
  ; varname = 'emiss_shipping'
  diam = 0.261e-6
  rho = 1770.
  mw = 115.
  mass_particle = rho *(!PI/6.) *(diam)^3  ;mass per particle (kg/particle)
  num_ship_so4 = ship_so4 *mw /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
  ;varname = 'emiss_res_tran'
  diam = 0.0504e-6
  rho = 1770.
  mw = 115.
  mass_particle = rho *(!PI/6.) *(diam)^3  ;mass per particle (kg/particle)
  num_res_tran_so4 = res_tran_so4 *mw /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
  ;varname = 'emiss_ene_ind'
  diam = 0.261e-6
  rho = 1770.
  mw = 115.
  mass_particle = rho *(!PI/6.) *(diam)^3  ;mass per particle (kg/particle)
  num_ene_ind_so4_surf = ene_ind_so4_surf *mw /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
  num_ene_ind_so4_vert = ene_ind_so4_vert *mw /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
;----------

;---------
;Write all the files
;---------
  ; write SO2 surface - ag+solv+waste, shipping, resid+trans
  spec = 'SO2'
  newfile = path_out+'emissions-cmip6_'+spec+'_anthro-ag-ship-res_surface_1750-2015_0.9x1.25_c'+creation_date+'.nc'
  print,newfile
  ncid = ncdf_create(newfile,/clobber)
  xid = ncdf_dimdef(ncid,'lon',nlon)
  yid = ncdf_dimdef(ncid,'lat',nlat)
  tid = ncdf_dimdef(ncid,'time',/unlimited)
  ; Define variables with attributes
  xvarid = ncdf_vardef(ncid,'lon',[xid],/float)
  ncdf_attput, ncid, xvarid,/char, 'units', 'degrees_east'
  ncdf_attput, ncid, xvarid,/char, 'long_name', 'Longitude'
  yvarid = ncdf_vardef(ncid,'lat',[yid],/float)
  ncdf_attput, ncid, yvarid,/char, 'units', 'degrees_north'
  ncdf_attput, ncid, yvarid,/char, 'long_name', 'Latitude'
  tvarid = ncdf_vardef(ncid,'time',[tid],/float)
  ncdf_attput, ncid, tvarid,/char, 'long_name', 'time'
  ncdf_attput, ncid, tvarid,/char, 'units', 'days since 1750-01-01 00:00:00'
  ncdf_attput, ncid, tvarid,/char, 'calendar', 'Gregorian'
  tvarid = ncdf_vardef(ncid,'date',[tid],/long)
  ncdf_attput, ncid, tvarid,/char, 'units', 'YYYYMMDD'
  ncdf_attput, ncid, tvarid,/char, 'long_name', 'Date'
  varid = ncdf_vardef(ncid,'emiss_ag_sol_was',[xid,yid,tid],/float)
  ncdf_attput,ncid,/char,varid,'units','molecules/cm2/s'
  ncdf_attput,ncid,/char,varid,'long_name','CEDS SO2*0.975 anthro agriculture+solvents+waste emissions'
  ncdf_attput,ncid,/float,varid,'molecular_weight',mw_so2
  varid = ncdf_vardef(ncid,'emiss_ship',[xid,yid,tid],/float)
  ncdf_attput,ncid,/char,varid,'units','molecules/cm2/s'
  ncdf_attput,ncid,/char,varid,'long_name','CEDS SO2*0.975 anthro shipping emissions'
  ncdf_attput,ncid,/float,varid,'molecular_weight',mw_so2
  varid = ncdf_vardef(ncid,'emiss_res_tran',[xid,yid,tid],/float)
  ncdf_attput,ncid,/char,varid,'units','molecules/cm2/s'
  ncdf_attput,ncid,/char,varid,'long_name','CEDS SO2*0.975 anthro residential+transportation emissions'
  ncdf_attput,ncid,/float,varid,'molecular_weight',mw_so2
  ;Define global attributes
  ncdf_attput,ncid,/GLOBAL,/char,'data_title','Anthropogenic surface emissions of '+spec+' for CMIP6'
  ncdf_attput,ncid,/GLOBAL,/float,'molecular_weight',mw_so2
  ncdf_attput,ncid,/GLOBAL,/char,'data_creator','Jean-Francois Lamarque (lamar@ucar.edu) and Louisa Emmons (emmons@ucar.edu)'
  ncdf_attput,ncid,/GLOBAL,/char,'data_summary','Emissions from the Community Emission Data System (CEDS) have been manipulated for use in CESM2 for CMIP6. Original file: '+source_file
  ncdf_attput,ncid,/GLOBAL,/char,'cesm_contact','Louisa Emmons or Simone Tilmes'
  ncdf_attput,ncid,/GLOBAL,/char,'creation_date',creation_date
  ncdf_attput,ncid,/GLOBAL,/char,'update_date',sdate_today
  ncdf_attput,ncid,/GLOBAL,/char,'history','original files have been regridded, units changed, concatenated, some species combined or renamed.'
  ncdf_attput,ncid,/GLOBAL,/char,'data_script',thisfile
  ncdf_attput,ncid,/GLOBAL,/char,'data_source_url','http://www.globalchange.umd.edu/ceds/ceds-cmip6-data/'
  ncdf_attput,ncid,/GLOBAL,/char,'data_reference','Hoesly et al., GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)'

  ncdf_control,ncid,/ENDEF
  ncdf_varput,ncid,'lon',lon
  ncdf_varput,ncid,'lat',lat
  ncdf_varput,ncid,'time',time
  ncdf_varput,ncid,'date',date
  ncdf_varput,ncid,'emiss_ag_sol_was',ag_sol_was_so2
  ncdf_varput,ncid,'emiss_ship',ship_so2
  ncdf_varput,ncid,'emiss_res_tran',res_tran_so2
  ncdf_close,ncid

  ; write SO2 surface - energy+industry
  spec = 'SO2'
  newfile = path_out+'emissions-cmip6_'+spec+'_anthro-ene_surface_1750-2015_0.9x1.25_c'+creation_date+'.nc'
  print,newfile
  ncid = ncdf_create(newfile,/clobber)
  xid = ncdf_dimdef(ncid,'lon',nlon)
  yid = ncdf_dimdef(ncid,'lat',nlat)
  tid = ncdf_dimdef(ncid,'time',/unlimited)
  ; Define variables with attributes
  xvarid = ncdf_vardef(ncid,'lon',[xid],/float)
  ncdf_attput, ncid, xvarid,/char, 'units', 'degrees_east'
  ncdf_attput, ncid, xvarid,/char, 'long_name', 'Longitude'
  yvarid = ncdf_vardef(ncid,'lat',[yid],/float)
  ncdf_attput, ncid, yvarid,/char, 'units', 'degrees_north'
  ncdf_attput, ncid, yvarid,/char, 'long_name', 'Latitude'
  tvarid = ncdf_vardef(ncid,'time',[tid],/float)
  ncdf_attput, ncid, tvarid,/char, 'long_name', 'time'
  ncdf_attput, ncid, tvarid,/char, 'units', 'days since 1750-01-01 00:00:00'
  ncdf_attput, ncid, tvarid,/char, 'calendar', 'Gregorian'
  tvarid = ncdf_vardef(ncid,'date',[tid],/long)
  ncdf_attput, ncid, tvarid,/char, 'units', 'YYYYMMDD'
  ncdf_attput, ncid, tvarid,/char, 'long_name', 'Date'
  varid = ncdf_vardef(ncid,'emiss_ene_ind',[xid,yid,tid],/float)
  ncdf_attput,ncid,/char,varid,'units','molecules/cm2/s'
  ncdf_attput,ncid,/char,varid,'long_name','CEDS SO2*0.975 anthro energy+industrial emissions'
  ncdf_attput,ncid,/float,varid,'molecular_weight',mw_so2
  ;Define global attributes
  ncdf_attput,ncid,/GLOBAL,/char,'data_title','Anthropogenic surface emissions of '+spec+' for CMIP6'
  ncdf_attput,ncid,/GLOBAL,/float,'molecular_weight',mw_so2
  ncdf_attput,ncid,/GLOBAL,/char,'data_creator','Jean-Francois Lamarque (lamar@ucar.edu) and Louisa Emmons (emmons@ucar.edu)'
  ncdf_attput,ncid,/GLOBAL,/char,'data_summary','Emissions from the Community Emission Data System (CEDS) have been manipulated for use in CESM2 for CMIP6. Original file: '+source_file
  ncdf_attput,ncid,/GLOBAL,/char,'cesm_contact','Louisa Emmons or Simone Tilmes'
  ncdf_attput,ncid,/GLOBAL,/char,'creation_date',creation_date
  ncdf_attput,ncid,/GLOBAL,/char,'update_date',sdate_today
  ncdf_attput,ncid,/GLOBAL,/char,'history','original files have been regridded, units changed, concatenated, some species combined or renamed.'
  ncdf_attput,ncid,/GLOBAL,/char,'data_script',thisfile
  ncdf_attput,ncid,/GLOBAL,/char,'data_source_url','http://www.globalchange.umd.edu/ceds/ceds-cmip6-data/'
  ncdf_attput,ncid,/GLOBAL,/char,'data_reference','Hoesly et al., GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)'

  ncdf_control,ncid,/ENDEF
  ncdf_varput,ncid,'lon',lon
  ncdf_varput,ncid,'lat',lat
  ncdf_varput,ncid,'time',time
  ncdf_varput,ncid,'date',date
  ncdf_varput,ncid,'emiss_ene_ind',ene_ind_so2_surf
  ncdf_close,ncid

  ; write SO2 vertical - energy+industry
  spec = 'SO2'
  newfile = path_out+'emissions-cmip6_'+spec+'_anthro-ene_vertical_1750-2015_0.9x1.25_c'+creation_date+'.nc'
  print,newfile
  ncid = ncdf_create(newfile,/clobber)
  xid = ncdf_dimdef(ncid,'lon',nlon)
  yid = ncdf_dimdef(ncid,'lat',nlat)
  zid = ncdf_dimdef(ncid,'altitude',nalt)
  zbid = ncdf_dimdef(ncid,'altitude_int',nalt+1)
  tid = ncdf_dimdef(ncid,'time',/unlimited)
  ; Define variables with attributes
  xvarid = ncdf_vardef(ncid,'lon',[xid],/float)
  ncdf_attput, ncid, xvarid,/char, 'units', 'degrees_east'
  ncdf_attput, ncid, xvarid,/char, 'long_name', 'Longitude'
  yvarid = ncdf_vardef(ncid,'lat',[yid],/float)
  ncdf_attput, ncid, yvarid,/char, 'units', 'degrees_north'
  ncdf_attput, ncid, yvarid,/char, 'long_name', 'Latitude'
  zvarid = ncdf_vardef(ncid,'altitude',[zid],/float)
  ncdf_attput, ncid, zvarid,/char,'units','km'
  ncdf_attput, ncid, zvarid, /char, 'long_name','Altitude'
  zbvarid = ncdf_vardef(ncid,'altitude_int',[zbid],/float)
  ncdf_attput, ncid, zbvarid,/char,'units','km'
  ncdf_attput, ncid, zbvarid, /char, 'long_name','Altitude interfaces'
  tvarid = ncdf_vardef(ncid,'time',[tid],/float)
  ncdf_attput, ncid, tvarid,/char, 'long_name', 'time'
  ncdf_attput, ncid, tvarid,/char, 'units', 'days since 1750-01-01 00:00:00'
  ncdf_attput, ncid, tvarid,/char, 'calendar', 'Gregorian'
  tvarid = ncdf_vardef(ncid,'date',[tid],/long)
  ncdf_attput, ncid, tvarid,/char, 'units', 'YYYYMMDD'
  ncdf_attput, ncid, tvarid,/char, 'long_name', 'Date'
  varid = ncdf_vardef(ncid,'emiss_ene_ind',[xid,yid,zid,tid],/float)
  ncdf_attput,ncid,/char,varid,'units','molecules/cm3/s'
  ncdf_attput,ncid,/char,varid,'long_name','CEDS SO2*0.975 anthro energy+industrial emissions'
  ncdf_attput,ncid,/float,varid,'molecular_weight',mw_so2
  ;Define global attributes
  ncdf_attput,ncid,/GLOBAL,/char,'data_title','Anthropogenic vertical emissions of '+spec+' for CMIP6'
  ncdf_attput,ncid,/GLOBAL,/float,'molecular_weight',mw_so2
  ncdf_attput,ncid,/GLOBAL,/char,'data_creator','Jean-Francois Lamarque (lamar@ucar.edu) and Louisa Emmons (emmons@ucar.edu)'
  ncdf_attput,ncid,/GLOBAL,/char,'data_summary','Emissions from the Community Emission Data System (CEDS) have been manipulated for use in CESM2 for CMIP6. Original file: '+source_file
  ncdf_attput,ncid,/GLOBAL,/char,'cesm_contact','Louisa Emmons or Simone Tilmes'
  ncdf_attput,ncid,/GLOBAL,/char,'creation_date',creation_date
  ncdf_attput,ncid,/GLOBAL,/char,'update_date',sdate_today
  ncdf_attput,ncid,/GLOBAL,/char,'history','original files have been regridded, units changed, concatenated, some species combined or renamed.'
  ncdf_attput,ncid,/GLOBAL,/char,'data_script',thisfile
  ncdf_attput,ncid,/GLOBAL,/char,'data_source_url','http://www.globalchange.umd.edu/ceds/ceds-cmip6-data/'
  ncdf_attput,ncid,/GLOBAL,/char,'data_reference','Hoesly et al., GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)'
  ncdf_control,ncid,/ENDEF
  ncdf_varput,ncid,'lon',lon
  ncdf_varput,ncid,'lat',lat
  ncdf_varput,ncid,'altitude',altitude
  ncdf_varput,ncid,'altitude_int',altitude_int
  ncdf_varput,ncid,'time',time
  ncdf_varput,ncid,'date',date
  ncdf_varput,ncid,'emiss_ene_ind',ene_ind_so2_vert
  ncdf_close,ncid

  ; write new file - so4_a1: sulfate accumulation mode for ag+waste+solvents and shipping
  spec = 'so4_a1'
  newfile = path_out+'emissions-cmip6_'+spec+'_anthro-ag-ship_surface_1750-2015_0.9x1.25_c'+creation_date+'.nc'
  print,newfile
  ncid = ncdf_create(newfile,/clobber)
  xid = ncdf_dimdef(ncid,'lon',nlon)
  yid = ncdf_dimdef(ncid,'lat',nlat)
  tid = ncdf_dimdef(ncid,'time',/unlimited)
  ; Define variables with attributes
  xvarid = ncdf_vardef(ncid,'lon',[xid],/float)
  ncdf_attput, ncid, xvarid,/char, 'units', 'degrees_east'
  ncdf_attput, ncid, xvarid,/char, 'long_name', 'Longitude'
  yvarid = ncdf_vardef(ncid,'lat',[yid],/float)
  ncdf_attput, ncid, yvarid,/char, 'units', 'degrees_north'
  ncdf_attput, ncid, yvarid,/char, 'long_name', 'Latitude'
  tvarid = ncdf_vardef(ncid,'time',[tid],/float)
  ncdf_attput, ncid, tvarid,/char, 'long_name', 'time'
  ncdf_attput, ncid, tvarid,/char, 'units', 'days since 1750-01-01 00:00:00'
  ncdf_attput, ncid, tvarid,/char, 'calendar', 'Gregorian'
  tvarid = ncdf_vardef(ncid,'date',[tid],/long)
  ncdf_attput, ncid, tvarid,/char, 'units', 'YYYYMMDD'
  ncdf_attput, ncid, tvarid,/char, 'long_name', 'Date'
  varid = ncdf_vardef(ncid,'emiss_ag_sol_was',[xid,yid,tid],/float)
  ncdf_attput,ncid,/char,varid,'units','molecules/cm2/s'
  ncdf_attput,ncid,/char,varid,'long_name','CEDS sulfate aerosol (0.025*SO2) anthro agriculture+solvents+waste emissions'
  ncdf_attput,ncid,/float,varid,'molecular_weight',mw_so4_mam
  varid = ncdf_vardef(ncid,'emiss_shipping',[xid,yid,tid],/float)
  ncdf_attput,ncid,/char,varid,'units','molecules/cm2/s'
  ncdf_attput,ncid,/char,varid,'long_name','CEDS sulfate aerosol (0.025*SO2) anthro shipping emissions'
  ncdf_attput,ncid,/float,varid,'molecular_weight',mw_so4_mam
  ;Define global attributes
  ncdf_attput,ncid,/GLOBAL,/char,'data_title','Anthropogenic surface emissions of '+spec+' for CMIP6'
  ncdf_attput,ncid,/GLOBAL,/float,'molecular_weight',mw_so4_mam
  ncdf_attput,ncid,/GLOBAL,/char,'data_creator','Jean-Francois Lamarque (lamar@ucar.edu) and Louisa Emmons (emmons@ucar.edu)'
  ncdf_attput,ncid,/GLOBAL,/char,'data_summary','Emissions from the Community Emission Data System (CEDS) have been manipulated for use in CESM2 for CMIP6.  Sulfate is 0.025*SO2. Original file: '+source_file
  ncdf_attput,ncid,/GLOBAL,/char,'cesm_contact','Louisa Emmons or Simone Tilmes'
  ncdf_attput,ncid,/GLOBAL,/char,'creation_date',creation_date
  ncdf_attput,ncid,/GLOBAL,/char,'update_date',sdate_today
  ncdf_attput,ncid,/GLOBAL,/char,'history','original files have been regridded, units changed, concatenated, some species combined or renamed.'
  ncdf_attput,ncid,/GLOBAL,/char,'data_script',thisfile
  ncdf_attput,ncid,/GLOBAL,/char,'data_source_url','http://www.globalchange.umd.edu/ceds/ceds-cmip6-data/'
  ncdf_attput,ncid,/GLOBAL,/char,'data_reference','Hoesly et al., GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)'
  ncdf_control,ncid,/ENDEF
  ncdf_varput,ncid,'lon',lon
  ncdf_varput,ncid,'lat',lat
  ncdf_varput,ncid,'time',time
  ncdf_varput,ncid,'date',date
  ncdf_varput,ncid,'emiss_ag_sol_was',ag_sol_was_so4
  ncdf_varput,ncid,'emiss_shipping',ship_so4
  ncdf_close,ncid

  ; write new file - num_so4_a1: sulfate accumulation mode for ag+waste+solvents and shipping
  spec = 'num_so4_a1'
  newfile = path_out+'emissions-cmip6_'+spec+'_anthro-ag-ship_surface_1750-2015_0.9x1.25_c'+creation_date+'.nc'
  print,newfile
  ncid = ncdf_create(newfile,/clobber)
  xid = ncdf_dimdef(ncid,'lon',nlon)
  yid = ncdf_dimdef(ncid,'lat',nlat)
  tid = ncdf_dimdef(ncid,'time',/unlimited)
  ; Define variables with attributes
  xvarid = ncdf_vardef(ncid,'lon',[xid],/float)
  ncdf_attput, ncid, xvarid,/char, 'units', 'degrees_east'
  ncdf_attput, ncid, xvarid,/char, 'long_name', 'Longitude'
  yvarid = ncdf_vardef(ncid,'lat',[yid],/float)
  ncdf_attput, ncid, yvarid,/char, 'units', 'degrees_north'
  ncdf_attput, ncid, yvarid,/char, 'long_name', 'Latitude'
  tvarid = ncdf_vardef(ncid,'time',[tid],/float)
  ncdf_attput, ncid, tvarid,/char, 'long_name', 'time'
  ncdf_attput, ncid, tvarid,/char, 'units', 'days since 1750-01-01 00:00:00'
  ncdf_attput, ncid, tvarid,/char, 'calendar', 'Gregorian'
  tvarid = ncdf_vardef(ncid,'date',[tid],/long)
  ncdf_attput, ncid, tvarid,/char, 'units', 'YYYYMMDD'
  ncdf_attput, ncid, tvarid,/char, 'long_name', 'Date'
  varid = ncdf_vardef(ncid,'emiss_ag_sol_was',[xid,yid,tid],/float)
  ncdf_attput,ncid,/char,varid,'units','(particles/cm2/s)(molecules/mole)(g/kg)'
  ncdf_attput,ncid,/char,varid,'long_name','particle number emissions of so4_a1 ag_sol_was'
  ncdf_attput,ncid,/float,varid,'molecular_weight',mw_num
  varid = ncdf_vardef(ncid,'emiss_shipping',[xid,yid,tid],/float)
  ncdf_attput,ncid,/char,varid,'units','(particles/cm2/s)(molecules/mole)(g/kg)'
  ncdf_attput,ncid,/char,varid,'long_name','particle number emissions of so4_a1 shipping'
  ncdf_attput,ncid,/float,varid,'molecular_weight',mw_num
  ;Define global attributes
  ncdf_attput,ncid,/GLOBAL,/char,'data_title','Anthropogenic surface emissions of '+spec+' for CMIP6'
  ncdf_attput,ncid,/GLOBAL,/float,'molecular_weight',mw_num
  ncdf_attput,ncid,/GLOBAL,/char,'data_creator','Jean-Francois Lamarque (lamar@ucar.edu) and Louisa Emmons (emmons@ucar.edu)'
  ncdf_attput,ncid,/GLOBAL,/char,'data_summary','Emissions of number of particles for corresponding aerosol emissions for MAM4 (bc, pom, so4). These MAM number emissions are scaled by 6.02E26 (Avog.*g/kg) to account for scaling applied in CESM/MOZART when read into model. '
  ncdf_attput,ncid,/GLOBAL,/char,'cesm_contact','Louisa Emmons or Simone Tilmes'
  ncdf_attput,ncid,/GLOBAL,/char,'creation_date',creation_date
  ncdf_attput,ncid,/GLOBAL,/char,'update_date',sdate_today
  ncdf_attput,ncid,/GLOBAL,/char,'history','original files have been regridded, units changed, concatenated, some species combined or renamed.'
  ncdf_attput,ncid,/GLOBAL,/char,'data_script',thisfile
  ncdf_attput,ncid,/GLOBAL,/char,'data_source_url','http://www.globalchange.umd.edu/ceds/ceds-cmip6-data/'
  ncdf_attput,ncid,/GLOBAL,/char,'data_reference','Hoesly et al., GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)'
  ncdf_control,ncid,/ENDEF
  ncdf_varput,ncid,'lon',lon
  ncdf_varput,ncid,'lat',lat
  ncdf_varput,ncid,'time',time
  ncdf_varput,ncid,'date',date
  ncdf_varput,ncid,'emiss_ag_sol_was',num_ag_sol_was_so4
  ncdf_varput,ncid,'emiss_shipping',num_ship_so4
  ncdf_close,ncid

  ; write new file - so4_a2: sulfate Aitken mode for residential+transportation
  spec = 'so4_a2'
  newfile = path_out+'emissions-cmip6_'+spec+'_anthro-res_surface_1750-2015_0.9x1.25_c'+creation_date+'.nc'
  print,newfile
  ncid = ncdf_create(newfile,/clobber)
  xid = ncdf_dimdef(ncid,'lon',nlon)
  yid = ncdf_dimdef(ncid,'lat',nlat)
  tid = ncdf_dimdef(ncid,'time',/unlimited)
  ; Define variables with attributes
  xvarid = ncdf_vardef(ncid,'lon',[xid],/float)
  ncdf_attput, ncid, xvarid,/char, 'units', 'degrees_east'
  ncdf_attput, ncid, xvarid,/char, 'long_name', 'Longitude'
  yvarid = ncdf_vardef(ncid,'lat',[yid],/float)
  ncdf_attput, ncid, yvarid,/char, 'units', 'degrees_north'
  ncdf_attput, ncid, yvarid,/char, 'long_name', 'Latitude'
  tvarid = ncdf_vardef(ncid,'time',[tid],/float)
  ncdf_attput, ncid, tvarid,/char, 'long_name', 'time'
  ncdf_attput, ncid, tvarid,/char, 'units', 'days since 1750-01-01 00:00:00'
  ncdf_attput, ncid, tvarid,/char, 'calendar', 'Gregorian'
  tvarid = ncdf_vardef(ncid,'date',[tid],/long)
  ncdf_attput, ncid, tvarid,/char, 'units', 'YYYYMMDD'
  ncdf_attput, ncid, tvarid,/char, 'long_name', 'Date'
  varid = ncdf_vardef(ncid,'emiss_res_tran',[xid,yid,tid],/float)
  ncdf_attput,ncid,/char,varid,'units','molecules/cm2/s'
  ncdf_attput,ncid,/char,varid,'long_name','CEDS sulfate aerosol (0.025*SO2) anthro residential+transportation emissions'
  ncdf_attput,ncid,/float,varid,'molecular_weight',mw_so4_mam
  ;Define global attributes
  ncdf_attput,ncid,/GLOBAL,/char,'data_title','Anthropogenic surface emissions of '+spec+' for CMIP6'
  ncdf_attput,ncid,/GLOBAL,/float,'molecular_weight',mw_so4_mam
  ncdf_attput,ncid,/GLOBAL,/char,'data_creator','Jean-Francois Lamarque (lamar@ucar.edu) and Louisa Emmons (emmons@ucar.edu)'
  ncdf_attput,ncid,/GLOBAL,/char,'data_summary','Emissions from the Community Emission Data System (CEDS) have been manipulated for use in CESM2 for CMIP6.  Sulfate is 0.025*SO2. Original file: '+source_file
  ncdf_attput,ncid,/GLOBAL,/char,'cesm_contact','Louisa Emmons or Simone Tilmes'
  ncdf_attput,ncid,/GLOBAL,/char,'creation_date',creation_date
  ncdf_attput,ncid,/GLOBAL,/char,'update_date',sdate_today
  ncdf_attput,ncid,/GLOBAL,/char,'history','original files have been regridded, units changed, concatenated, some species combined or renamed.'
  ncdf_attput,ncid,/GLOBAL,/char,'data_script',thisfile
  ncdf_attput,ncid,/GLOBAL,/char,'data_source_url','http://www.globalchange.umd.edu/ceds/ceds-cmip6-data/'
  ncdf_attput,ncid,/GLOBAL,/char,'data_reference','Hoesly et al., GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)'
  ncdf_control,ncid,/ENDEF
  ncdf_varput,ncid,'lon',lon
  ncdf_varput,ncid,'lat',lat
  ncdf_varput,ncid,'time',time
  ncdf_varput,ncid,'date',date
  ncdf_varput,ncid,'emiss_res_tran',res_tran_so4
  ncdf_close,ncid

  ; write new file - num_so4_a2: sulfate Aitken mode for residential+transportation
  spec = 'num_so4_a2'
  newfile = path_out+'emissions-cmip6_'+spec+'_anthro-res_surface_1750-2015_0.9x1.25_c'+creation_date+'.nc'
  print,newfile
  ncid = ncdf_create(newfile,/clobber)
  xid = ncdf_dimdef(ncid,'lon',nlon)
  yid = ncdf_dimdef(ncid,'lat',nlat)
  tid = ncdf_dimdef(ncid,'time',/unlimited)
  ; Define variables with attributes
  xvarid = ncdf_vardef(ncid,'lon',[xid],/float)
  ncdf_attput, ncid, xvarid,/char, 'units', 'degrees_east'
  ncdf_attput, ncid, xvarid,/char, 'long_name', 'Longitude'
  yvarid = ncdf_vardef(ncid,'lat',[yid],/float)
  ncdf_attput, ncid, yvarid,/char, 'units', 'degrees_north'
  ncdf_attput, ncid, yvarid,/char, 'long_name', 'Latitude'
  tvarid = ncdf_vardef(ncid,'time',[tid],/float)
  ncdf_attput, ncid, tvarid,/char, 'long_name', 'time'
  ncdf_attput, ncid, tvarid,/char, 'units', 'days since 1750-01-01 00:00:00'
  ncdf_attput, ncid, tvarid,/char, 'calendar', 'Gregorian'
  tvarid = ncdf_vardef(ncid,'date',[tid],/long)
  ncdf_attput, ncid, tvarid,/char, 'units', 'YYYYMMDD'
  ncdf_attput, ncid, tvarid,/char, 'long_name', 'Date'
  varid = ncdf_vardef(ncid,'emiss_res_tran',[xid,yid,tid],/float)
  ncdf_attput,ncid,/char,varid,'units','(particles/cm2/s)(molecules/mole)(g/kg)'
  ncdf_attput,ncid,/char,varid,'long_name','particle number emissions of so4_a1 residential+transportation'
  ncdf_attput,ncid,/float,varid,'molecular_weight',mw_num
  ;Define global attributes
  ncdf_attput,ncid,/GLOBAL,/char,'data_title','Anthropogenic surface emissions of '+spec+' for CMIP6'
  ncdf_attput,ncid,/GLOBAL,/float,'molecular_weight',mw_num
  ncdf_attput,ncid,/GLOBAL,/char,'data_creator','Jean-Francois Lamarque (lamar@ucar.edu) and Louisa Emmons (emmons@ucar.edu)'
  ncdf_attput,ncid,/GLOBAL,/char,'data_summary','Emissions of number of particles for corresponding aerosol emissions for MAM4 (bc, pom, so4). These MAM number emissions are scaled by 6.02E26 (Avog.*g/kg) to account for scaling applied in CESM/MOZART when read into model. '
  ncdf_attput,ncid,/GLOBAL,/char,'cesm_contact','Louisa Emmons or Simone Tilmes'
  ncdf_attput,ncid,/GLOBAL,/char,'creation_date',creation_date
  ncdf_attput,ncid,/GLOBAL,/char,'update_date',sdate_today
  ncdf_attput,ncid,/GLOBAL,/char,'history','original files have been regridded, units changed, concatenated, some species combined or renamed.'
  ncdf_attput,ncid,/GLOBAL,/char,'data_script',thisfile
  ncdf_attput,ncid,/GLOBAL,/char,'data_source_url','http://www.globalchange.umd.edu/ceds/ceds-cmip6-data/'
  ncdf_attput,ncid,/GLOBAL,/char,'data_reference','Hoesly et al., GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)'
  ncdf_control,ncid,/ENDEF
  ncdf_varput,ncid,'lon',lon
  ncdf_varput,ncid,'lat',lat
  ncdf_varput,ncid,'time',time
  ncdf_varput,ncid,'date',date
  ncdf_varput,ncid,'emiss_res_tran',num_res_tran_so4
  ncdf_close,ncid

  ; write so4_a1: sulfate accumulation mode for energy+industrial
  ; at surface
  spec = 'so4_a1'
  newfile = path_out+'emissions-cmip6_'+spec+'_anthro-ene_surface_1750-2015_0.9x1.25_c'+creation_date+'.nc'
  print,newfile
  ncid = ncdf_create(newfile,/clobber)
  xid = ncdf_dimdef(ncid,'lon',nlon)
  yid = ncdf_dimdef(ncid,'lat',nlat)
  tid = ncdf_dimdef(ncid,'time',/unlimited)
  ; Define variables with attributes
  xvarid = ncdf_vardef(ncid,'lon',[xid],/float)
  ncdf_attput, ncid, xvarid,/char, 'units', 'degrees_east'
  ncdf_attput, ncid, xvarid,/char, 'long_name', 'Longitude'
  yvarid = ncdf_vardef(ncid,'lat',[yid],/float)
  ncdf_attput, ncid, yvarid,/char, 'units', 'degrees_north'
  ncdf_attput, ncid, yvarid,/char, 'long_name', 'Latitude'
  tvarid = ncdf_vardef(ncid,'time',[tid],/float)
  ncdf_attput, ncid, tvarid,/char, 'long_name', 'time'
  ncdf_attput, ncid, tvarid,/char, 'units', 'days since 1750-01-01 00:00:00'
  ncdf_attput, ncid, tvarid,/char, 'calendar', 'Gregorian'
  tvarid = ncdf_vardef(ncid,'date',[tid],/long)
  ncdf_attput, ncid, tvarid,/char, 'units', 'YYYYMMDD'
  ncdf_attput, ncid, tvarid,/char, 'long_name', 'Date'
  varid = ncdf_vardef(ncid,'emiss_ene_ind',[xid,yid,tid],/float)
  ncdf_attput,ncid,/char,varid,'units','molecules/cm2/s'
  ncdf_attput,ncid,/char,varid,'long_name','CEDS sulfate aerosol (0.025*SO2) anthro energy+industry emissions'
  ncdf_attput,ncid,/float,varid,'molecular_weight',mw_so4_mam
  ;Define global attributes
  ncdf_attput,ncid,/GLOBAL,/char,'data_title','Anthropogenic surface emissions of '+spec+' for CMIP6'
  ncdf_attput,ncid,/GLOBAL,/float,'molecular_weight',mw_so4_mam
  ncdf_attput,ncid,/GLOBAL,/char,'data_creator','Jean-Francois Lamarque (lamar@ucar.edu) and Louisa Emmons (emmons@ucar.edu)'
  ncdf_attput,ncid,/GLOBAL,/char,'data_summary','Emissions from the Community Emission Data System (CEDS) have been manipulated for use in CESM2 for CMIP6.  Sulfate is 0.025*SO2. Original file: '+source_file
  ncdf_attput,ncid,/GLOBAL,/char,'cesm_contact','Louisa Emmons or Simone Tilmes'
  ncdf_attput,ncid,/GLOBAL,/char,'creation_date',creation_date
  ncdf_attput,ncid,/GLOBAL,/char,'update_date',sdate_today
  ncdf_attput,ncid,/GLOBAL,/char,'history','original files have been regridded, units changed, concatenated, some species combined or renamed.'
  ncdf_attput,ncid,/GLOBAL,/char,'data_script',thisfile
  ncdf_attput,ncid,/GLOBAL,/char,'data_source_url','http://www.globalchange.umd.edu/ceds/ceds-cmip6-data/'
  ncdf_attput,ncid,/GLOBAL,/char,'data_reference','Hoesly et al., GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)'
  ncdf_control,ncid,/ENDEF
  ncdf_varput,ncid,'lon',lon
  ncdf_varput,ncid,'lat',lat
  ncdf_varput,ncid,'time',time
  ncdf_varput,ncid,'date',date
  ncdf_varput,ncid,'emiss_ene_ind',ene_ind_so4_surf
  ncdf_close,ncid

  ; write num_so4_a1: sulfate accumulation mode for energy+industrial
  ; at surface
  spec = 'num_so4_a1'
  newfile = path_out+'emissions-cmip6_'+spec+'_anthro-ene_surface_1750-2015_0.9x1.25_c'+creation_date+'.nc'
  print,newfile
  ncid = ncdf_create(newfile,/clobber)
  xid = ncdf_dimdef(ncid,'lon',nlon)
  yid = ncdf_dimdef(ncid,'lat',nlat)
  tid = ncdf_dimdef(ncid,'time',/unlimited)
  ; Define variables with attributes
  xvarid = ncdf_vardef(ncid,'lon',[xid],/float)
  ncdf_attput, ncid, xvarid,/char, 'units', 'degrees_east'
  ncdf_attput, ncid, xvarid,/char, 'long_name', 'Longitude'
  yvarid = ncdf_vardef(ncid,'lat',[yid],/float)
  ncdf_attput, ncid, yvarid,/char, 'units', 'degrees_north'
  ncdf_attput, ncid, yvarid,/char, 'long_name', 'Latitude'
  tvarid = ncdf_vardef(ncid,'time',[tid],/float)
  ncdf_attput, ncid, tvarid,/char, 'long_name', 'time'
  ncdf_attput, ncid, tvarid,/char, 'units', 'days since 1750-01-01 00:00:00'
  ncdf_attput, ncid, tvarid,/char, 'calendar', 'Gregorian'
  tvarid = ncdf_vardef(ncid,'date',[tid],/long)
  ncdf_attput, ncid, tvarid,/char, 'units', 'YYYYMMDD'
  ncdf_attput, ncid, tvarid,/char, 'long_name', 'Date'
  varid = ncdf_vardef(ncid,'emiss_ene_ind',[xid,yid,tid],/float)
  ncdf_attput,ncid,/char,varid,'units','(particles/cm2/s)(molecules/mole)(g/kg)'
  ncdf_attput,ncid,/char,varid,'long_name','particle number emissions of so4_a1 energy+industry'
  ncdf_attput,ncid,/float,varid,'molecular_weight',mw_num
  ;Define global attributes
  ncdf_attput,ncid,/GLOBAL,/char,'data_title','Anthropogenic surface emissions of '+spec+' for CMIP6'
  ncdf_attput,ncid,/GLOBAL,/float,'molecular_weight',mw_num
  ncdf_attput,ncid,/GLOBAL,/char,'data_creator','Jean-Francois Lamarque (lamar@ucar.edu) and Louisa Emmons (emmons@ucar.edu)'
  ncdf_attput,ncid,/GLOBAL,/char,'data_summary','Emissions of number of particles for corresponding aerosol emissions for MAM4 (bc, pom, so4). These MAM number emissions are scaled by 6.02E26 (Avog.*g/kg) to account for scaling applied in CESM/MOZART when read into model. '
  ncdf_attput,ncid,/GLOBAL,/char,'cesm_contact','Louisa Emmons or Simone Tilmes'
  ncdf_attput,ncid,/GLOBAL,/char,'creation_date',creation_date
  ncdf_attput,ncid,/GLOBAL,/char,'update_date',sdate_today
  ncdf_attput,ncid,/GLOBAL,/char,'history','original files have been regridded, units changed, concatenated, some species combined or renamed.'
  ncdf_attput,ncid,/GLOBAL,/char,'data_script',thisfile
  ncdf_attput,ncid,/GLOBAL,/char,'data_source_url','http://www.globalchange.umd.edu/ceds/ceds-cmip6-data/'
  ncdf_attput,ncid,/GLOBAL,/char,'data_reference','Hoesly et al., GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)'
  ncdf_control,ncid,/ENDEF
  ncdf_varput,ncid,'lon',lon
  ncdf_varput,ncid,'lat',lat
  ncdf_varput,ncid,'time',time
  ncdf_varput,ncid,'date',date
  ncdf_varput,ncid,'emiss_ene_ind',num_ene_ind_so4_surf
  ncdf_close,ncid

  spec = 'so4_a1'
  newfile = path_out+'emissions-cmip6_'+spec+'_anthro-ene_vertical_1750-2015_0.9x1.25_c'+creation_date+'.nc'
  print,newfile
  ncid = ncdf_create(newfile,/clobber)
  xid = ncdf_dimdef(ncid,'lon',nlon)
  yid = ncdf_dimdef(ncid,'lat',nlat)
  zid = ncdf_dimdef(ncid,'altitude',nalt)
  zbid = ncdf_dimdef(ncid,'altitude_int',nalt+1)
  tid = ncdf_dimdef(ncid,'time',/unlimited)
  ; Define variables with attributes
  xvarid = ncdf_vardef(ncid,'lon',[xid],/float)
  ncdf_attput, ncid, xvarid,/char, 'units', 'degrees_east'
  ncdf_attput, ncid, xvarid,/char, 'long_name', 'Longitude'
  yvarid = ncdf_vardef(ncid,'lat',[yid],/float)
  ncdf_attput, ncid, yvarid,/char, 'units', 'degrees_north'
  ncdf_attput, ncid, yvarid,/char, 'long_name', 'Latitude'
  zvarid = ncdf_vardef(ncid,'altitude',[zid],/float)
  ncdf_attput, ncid, zvarid,/char,'units','km'
  ncdf_attput, ncid, zvarid, /char, 'long_name','Altitude'
  zbvarid = ncdf_vardef(ncid,'altitude_int',[zbid],/float)
  ncdf_attput, ncid, zbvarid,/char,'units','km'
  ncdf_attput, ncid, zbvarid, /char, 'long_name','Altitude interfaces'
  tvarid = ncdf_vardef(ncid,'time',[tid],/float)
  ncdf_attput, ncid, tvarid,/char, 'long_name', 'time'
  ncdf_attput, ncid, tvarid,/char, 'units', 'days since 1750-01-01 00:00:00'
  ncdf_attput, ncid, tvarid,/char, 'calendar', 'Gregorian'
  tvarid = ncdf_vardef(ncid,'date',[tid],/long)
  ncdf_attput, ncid, tvarid,/char, 'units', 'YYYYMMDD'
  ncdf_attput, ncid, tvarid,/char, 'long_name', 'Date'
  varid = ncdf_vardef(ncid,'emiss_ene_ind',[xid,yid,zid,tid],/float)
  ncdf_attput,ncid,/char,varid,'units','molecules/cm3/s'
  ncdf_attput,ncid,/char,varid,'long_name','CEDS sulfate aerosol (0.025*SO2) anthro energy+industrial emissions'
  ncdf_attput,ncid,/float,varid,'molecular_weight',mw_so4_mam
  ;Define global attributes
  ncdf_attput,ncid,/GLOBAL,/char,'data_title','Anthropogenic vertical emissions of '+spec+' for CMIP6'
  ncdf_attput,ncid,/GLOBAL,/float,'molecular_weight',mw_so4_mam
  ncdf_attput,ncid,/GLOBAL,/char,'data_creator','Jean-Francois Lamarque (lamar@ucar.edu) and Louisa Emmons (emmons@ucar.edu)'
  ncdf_attput,ncid,/GLOBAL,/char,'data_summary','Emissions from the Community Emission Data System (CEDS) have been manipulated for use in CESM2 for CMIP6.  Sulfate is 0.025*SO2. Original file: '+source_file
  ncdf_attput,ncid,/GLOBAL,/char,'cesm_contact','Louisa Emmons or Simone Tilmes'
  ncdf_attput,ncid,/GLOBAL,/char,'creation_date',creation_date
  ncdf_attput,ncid,/GLOBAL,/char,'update_date',sdate_today
  ncdf_attput,ncid,/GLOBAL,/char,'history','original files have been regridded, units changed, concatenated, some species combined or renamed.'
  ncdf_attput,ncid,/GLOBAL,/char,'data_script',thisfile
  ncdf_attput,ncid,/GLOBAL,/char,'data_source_url','http://www.globalchange.umd.edu/ceds/ceds-cmip6-data/'
  ncdf_attput,ncid,/GLOBAL,/char,'data_reference','Hoesly et al., GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)'
  ncdf_control,ncid,/ENDEF
  ncdf_varput,ncid,'lon',lon
  ncdf_varput,ncid,'lat',lat
  ncdf_varput,ncid,'altitude',altitude
  ncdf_varput,ncid,'altitude_int',altitude_int
  ncdf_varput,ncid,'time',time
  ncdf_varput,ncid,'date',date
  ncdf_varput,ncid,'emiss_ene_ind',ene_ind_so4_vert
  ncdf_close,ncid

  spec = 'num_so4_a1'
  newfile = path_out+'emissions-cmip6_'+spec+'_anthro-ene_vertical_1750-2015_0.9x1.25_c'+creation_date+'.nc'
  print,newfile
  ncid = ncdf_create(newfile,/clobber)
  xid = ncdf_dimdef(ncid,'lon',nlon)
  yid = ncdf_dimdef(ncid,'lat',nlat)
  zid = ncdf_dimdef(ncid,'altitude',nalt)
  zbid = ncdf_dimdef(ncid,'altitude_int',nalt+1)
  tid = ncdf_dimdef(ncid,'time',/unlimited)
  ; Define variables with attributes
  xvarid = ncdf_vardef(ncid,'lon',[xid],/float)
  ncdf_attput, ncid, xvarid,/char, 'units', 'degrees_east'
  ncdf_attput, ncid, xvarid,/char, 'long_name', 'Longitude'
  yvarid = ncdf_vardef(ncid,'lat',[yid],/float)
  ncdf_attput, ncid, yvarid,/char, 'units', 'degrees_north'
  ncdf_attput, ncid, yvarid,/char, 'long_name', 'Latitude'
  zvarid = ncdf_vardef(ncid,'altitude',[zid],/float)
  ncdf_attput, ncid, zvarid,/char,'units','km'
  ncdf_attput, ncid, zvarid, /char, 'long_name','Altitude'
  zbvarid = ncdf_vardef(ncid,'altitude_int',[zbid],/float)
  ncdf_attput, ncid, zbvarid,/char,'units','km'
  ncdf_attput, ncid, zbvarid, /char, 'long_name','Altitude interfaces'
  tvarid = ncdf_vardef(ncid,'time',[tid],/float)
  ncdf_attput, ncid, tvarid,/char, 'long_name', 'time'
  ncdf_attput, ncid, tvarid,/char, 'units', 'days since 1750-01-01 00:00:00'
  ncdf_attput, ncid, tvarid,/char, 'calendar', 'Gregorian'
  tvarid = ncdf_vardef(ncid,'date',[tid],/long)
  ncdf_attput, ncid, tvarid,/char, 'units', 'YYYYMMDD'
  ncdf_attput, ncid, tvarid,/char, 'long_name', 'Date'
  varid = ncdf_vardef(ncid,'emiss_ene_ind',[xid,yid,zid,tid],/float)
  ncdf_attput,ncid,/char,varid,'units','(particles/cm2/s)(molecules/mole)(g/kg)'
  ncdf_attput,ncid,/char,varid,'long_name','particle number emissions of so4_a1 energy+industry'
  ncdf_attput,ncid,/float,varid,'molecular_weight',mw_num
 ;Define global attributes
  ncdf_attput,ncid,/GLOBAL,/char,'data_title','Anthropogenic vertical emissions of '+spec+' for CMIP6'
  ncdf_attput,ncid,/GLOBAL,/float,'molecular_weight',mw_num
  ncdf_attput,ncid,/GLOBAL,/char,'data_creator','Jean-Francois Lamarque (lamar@ucar.edu) and Louisa Emmons (emmons@ucar.edu)'
  ncdf_attput,ncid,/GLOBAL,/char,'data_summary','Emissions of number of particles for corresponding aerosol emissions for MAM4 (bc, pom, so4). These MAM number emissions are scaled by 6.02E26 (Avog.*g/kg) to account for scaling applied in CESM/MOZART when read into model. '
  ncdf_attput,ncid,/GLOBAL,/char,'cesm_contact','Louisa Emmons or Simone Tilmes'
  ncdf_attput,ncid,/GLOBAL,/char,'creation_date',creation_date
  ncdf_attput,ncid,/GLOBAL,/char,'update_date',sdate_today
  ncdf_attput,ncid,/GLOBAL,/char,'history','original files have been regridded, units changed, concatenated, some species combined or renamed.'
  ncdf_attput,ncid,/GLOBAL,/char,'data_script',thisfile
  ncdf_attput,ncid,/GLOBAL,/char,'data_source_url','http://www.globalchange.umd.edu/ceds/ceds-cmip6-data/'
  ncdf_attput,ncid,/GLOBAL,/char,'data_reference','Hoesly et al., GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)'
  ncdf_control,ncid,/ENDEF
  ncdf_varput,ncid,'lon',lon
  ncdf_varput,ncid,'lat',lat
  ncdf_varput,ncid,'altitude',altitude
  ncdf_varput,ncid,'altitude_int',altitude_int
  ncdf_varput,ncid,'time',time
  ncdf_varput,ncid,'date',date
  ncdf_varput,ncid,'emiss_ene_ind',num_ene_ind_so4_vert
  ncdf_close,ncid

end

