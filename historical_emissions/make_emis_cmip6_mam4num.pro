; calculate number emissions for MAM bc,pom,so4
; various parameters for different modes and sectors

pro make_emis_cmip6_mam4num

today = bin_date(systime())
todaystr = String(today[0:2],format='(i4,"/",i2.2,"/",i2.2)')
sdate_today = String(today[0:2],format='(i4,i2.2,i2.2)')
creation_date = '20170608'
;creation_date = sdate_today
infiledate = '20170608'
thisfile = Routine_filepath()

;path_emis = '/glade/p/cesm/chwg_dev/emmons/CMIP6_emissions_1750_2015_v20170322/'
path_emis = '/glade/p/cesm/chwg_dev/emmons/CMIP6_emissions_1750_2015_v20170608/'

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

 avog=6.022e23  ;molecules/mole
 kg_g = 1.e-3   ;kg/g
 mw_num = 1.
 ;amufac = 1.65979e-23      ;(cm2/m2)(kg/g)(mole/molecules)
 ;when emissions are read in CESM, if not in kg/m2/s, they are scaled by mw*amufac
 ; so, emissions files must be scaled by 1./(mw*amufac)
 ; amufac = 1E4*1E-3/avog

 ;specs_num_surf = ['num_so4_a1_anthro','num_so4_a1_ship','num_so4_a2_anthro','num_so4_a1_bb', $
 ;                  'num_bc_a4_anthro','num_bc_a4_bb','num_pom_a4_anthro','num_pom_a4_bb']
 ;specs_num_surf = ['num_pom_a4_anthro','num_pom_a4_bb']
 specs_num_surf = ['num_bc_a4_anthro','num_pom_a4_anthro']
 ;nspsurf = n_elements(specs_num_surf)
 nspsurf = 0
 ;specs_num_vert = ['num_so4_a1_anthro','num_so4_a1_contvolcano','num_so4_a2_contvolcano', $
 ;  'num_bc_a4_aircraft','num_bc_a4_bb','num_so4_a1_bb','num_pom_a4_bb']
 ;specs_num_vert = ['num_pom_a4_bb']
 specs_num_vert = ['num_bc_a4_aircraft']
 nspvert = n_elements(specs_num_vert)
 ;nspvert = 0

 for ispec = 0,nspsurf-1 do begin
 
    spec = specs_num_surf[ispec]
    case spec of
       'num_so4_a1_anthro': begin
          filename = path_emis+'emissions-cmip6_so4_a1_anthro_surface_1750-2015_0.9x1.25_c'+infiledate+'2.nc'
          varname = 'emiss_ag_sol_was'
          diam = 0.134e-6
          rho = 1770.
          mw = 115.
       end
       'num_so4_a1_ship': begin
          filename = path_emis+'emissions-cmip6_so4_a1_anthro_surface_1750-2015_0.9x1.25_c'+infiledate+'.nc'
          varname = 'emiss_shipping'
          diam = 0.261e-6
          rho = 1770.
          mw = 115.
       end       
       'num_so4_a2_anthro': begin
          filename = path_emis+'emissions-cmip6_so4_a2_anthro_surface_1750-2015_0.9x1.25_c'+infiledate+'.nc'
          varname = 'emiss_res_tran'
          diam = 0.0504e-6
          rho = 1770.
          mw = 115.
       end       
       'num_so4_a1_bb': begin
          filename = path_emis+'emissions-cmip6_so4_a1_bb_surface_1750-2015_0.9x1.25_c'+infiledate+'.nc'
          varname = 'emiss_bb'
          diam = 0.134e-6
          rho = 1770.
          mw = 115.
       end       
       'num_bc_a4_anthro': begin
          filename = path_emis+'emissions-cmip6_bc_a4_anthro_surface_1750-2015_0.9x1.25_c'+infiledate+'.nc'
          varname = 'emiss_anthro'
          diam = 0.134e-6
          rho = 1700.
          mw = 12.
       end       
       'num_bc_a4_bb': begin
          filename = path_emis+'emissions-cmip6_bc_a4_bb_surface_1750-2015_0.9x1.25_c'+infiledate+'.nc'
          varname = 'emiss_bb'
          diam = 0.134e-6
          rho = 1700.
          mw = 12.
       end       
       'num_pom_a4_anthro': begin
          filename = path_emis+'emissions-cmip6_pom_a4_anthro_surface_1750-2015_0.9x1.25_c'+infiledate+'.nc'
          varname = 'emiss_anthro'
          diam = 0.134e-6
          rho = 1000.
          ;mw = 16.8
          mw = 12.
       end       
       'num_pom_a4_bb': begin
          filename = path_emis+'emissions-cmip6_pom_a4_bb_surface_1750-2015_0.9x1.25_c'+infiledate+'.nc'
          varname = 'emiss_bb'
          diam = 0.134e-6
          rho = 1000.
          ;mw = 16.8
          mw = 12.
       end       
    endcase
    
    ; read concentration emissions file for species
    print,filename
    ncid = ncdf_open(filename)
    ncdf_varget,ncid,'lon',lon
    ncdf_varget,ncid,'lat',lat
    ncdf_varget,ncid,'date',date
    ncdf_varget,ncid,'time',time
    ncdf_varget,ncid,varname, emis_mol ;molecules/cm2/s
    ncdf_attget,ncid,varname,'long_name',var_hist
    var_hist = string(var_hist)
    print,var_hist
    help,emis_mol
    ncdf_close,ncid
    nlon = n_elements(lon)
    nlat = n_elements(lat)
    ntim = n_elements(time)

    ; calculate number emissions
    mass_particle = rho *(!PI/6.) *(diam)^3  ;mass per particle (kg/particle)
    emis_num = emis_mol *mw  /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)

    ; write new file
    newfile = path_emis+'emissions-cmip6_'+spec+'_surface_1750-2015_0.9x1.25_c'+creation_date+'.nc'
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
  varid = ncdf_vardef(ncid,spec,[xid,yid,tid],/float)
  ncdf_attput,ncid,/char,varid,'units','(particles/cm2/s)(molecules/mole)(g/kg)'
  ncdf_attput,ncid,/char,varid,'long_name','particle number emissions of '+spec
  ncdf_attput,ncid,/char,varid,'history',var_hist
  ncdf_attput,ncid,/float,varid,'molecular_weight',mw_num
  ;Define global attributes
  ncdf_attput,ncid,/GLOBAL,/char,'data_title','Emissions of '+spec+' for CESM-MAM4 CMIP6'
  ncdf_attput,ncid,/GLOBAL,/float,'molecular_weight',mw_num
  ncdf_attput,ncid,/GLOBAL,/char,'data_creator','Jean-Francois Lamarque (lamar@ucar.edu) and Louisa Emmons (emmons@ucar.edu)'
  ncdf_attput,ncid,/GLOBAL,/char,'data_summary','Emissions of number of particles for corresponding aerosol emissions for MAM4 (bc, pom, so4). These MAM number emissions are scaled by 6.02E26 (Avog.*g/kg) to account for scaling applied in CESM/MOZART when read into model. '
  ncdf_attput,ncid,/GLOBAL,/char,'cesm_contact','Louisa Emmons or Simone Tilmes'
  ncdf_attput,ncid,/GLOBAL,/char,'creation_date',creation_date
  ncdf_attput,ncid,/GLOBAL,/char,'update_date',sdate_today
  ncdf_attput,ncid,/GLOBAL,/char,'history','Original aerosol file: '+filename
  ncdf_attput,ncid,/GLOBAL,/char,'data_script',thisfile
  ncdf_attput,ncid,/GLOBAL,/char,'data_source_url','http://www.globalchange.umd.edu/ceds/ceds-cmip6-data/'
  ncdf_attput,ncid,/GLOBAL,/char,'data_reference','Number calculations described in X.Liu et al., GMD, 2012.'
  ncdf_control,ncid,/ENDEF
  ncdf_varput,ncid,'lon',lon
  ncdf_varput,ncid,'lat',lat
  ncdf_varput,ncid,'time',time
  ncdf_varput,ncid,'date',date
  ncdf_varput,ncid,spec,emis_num
  ncdf_close,ncid
endfor

;vertical 
for ispec = 0,nspvert-1 do begin
    spec = specs_num_vert[ispec]
print,spec
    case spec of
       'num_so4_a1_anthro': begin
          filename = path_emis+'emissions-cmip6_so4_a1_anthro_vertical_1750-2015_0.9x1.25_c'+infiledate+'.nc'
          varname = 'emiss_ene_ind'
          diam = 0.261e-6
          rho = 1770.
          mw = 115.
       end
       'num_so4_a1_contvolcano': begin
          filename = path_emis+'emissions-cmip6_so4_a1_contvolcano_vertical_1750-2100_0.9x1.25_c'+infiledate+'.nc'
          varname = 'emiss_volcanoes'
          diam = 0.134e-6
          rho = 1770.
          mw = 115.
       end
       'num_so4_a2_contvolcano': begin
          filename = path_emis+'emissions-cmip6_so4_a2_contvolcano_vertical_1750-2100_0.9x1.25_c'+infiledate+'.nc'
          varname = 'emiss_volcanoes'
          diam = 0.0504e-6
          rho = 1770.
          mw = 115.
       end
       'num_bc_a4_aircraft': begin
          filename = path_emis+'emissions-cmip6_bc_a4_aircraft_vertical_1750-2015_0.9x1.25_c'+infiledate+'.nc'
          varname = 'emiss_aircraft'
          diam = 0.134e-6
          rho = 1700.
          mw = 12.
       end       

       'num_bc_a4_bb': begin
          filename = path_emis+'emissions-cmip6_bc_a4_bb_vertical_1750-2015_0.9x1.25_c'+infiledate+'.nc'
          varname = 'emiss_bb'
          diam = 0.134e-6
          rho = 1700.
          mw = 12.
       end       
       'num_pom_a4_bb': begin
          filename = path_emis+'emissions-cmip6_pom_a4_bb_vertical_1750-2015_0.9x1.25_c'+infiledate+'.nc'
          varname = 'emiss_bb'
          diam = 0.134e-6
          rho = 1000.
          ;mw = 16.8
          mw = 12.
       end       
       'num_so4_a1_bb': begin
          filename = path_emis+'emissions-cmip6_so4_a1_bb_vertical_1750-2015_0.9x1.25_c'+infiledate+'.nc'
          varname = 'emiss_bb'
          diam = 0.134e-6
          rho = 1770.
          mw = 115.
       end       
     endcase
    ; read concentration emissions file for species
    print,filename
    ncid = ncdf_open(filename)
    ncdf_varget,ncid,'lon',lon
    ncdf_varget,ncid,'lat',lat
    ncdf_varget,ncid,'altitude',altitude
    ncdf_varget,ncid,'altitude_int',altitude_int
    ncdf_varget,ncid,'date',date
    ncdf_varget,ncid,'time',time
    ncdf_varget,ncid,varname, emis_mol ;molecules/cm2/s
    ncdf_attget,ncid,varname,'long_name',var_hist
    var_hist = string(var_hist)
    print,var_hist
    help,emis_mol
    ncdf_close,ncid
    nlon = n_elements(lon)
    nlat = n_elements(lat)
    nalt = n_elements(altitude)
    ntim = n_elements(time)

    ; calculate number emissions
    mass_particle = rho *(!PI/6.) *(diam)^3  ;mass per particle (kg/particle)
    emis_num = emis_mol *mw  /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)

    ; write new file
    if (spec eq 'num_so4_a1_contvolcano' or spec eq 'num_so4_a2_contvolcano') then $
      newfile = path_emis+'emissions-cmip6_'+spec+'_vertical_1750-2100_0.9x1.25_c'+creation_date+'.nc' else $
      newfile = path_emis+'emissions-cmip6_'+spec+'_vertical_1750-2015_0.9x1.25_c'+creation_date+'.nc'
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
  varid = ncdf_vardef(ncid,spec,[xid,yid,zid,tid],/float)
  ncdf_attput,ncid,/char,varid,'units','(particles/cm2/s)(molecules/mole)(g/kg)'
  ncdf_attput,ncid,/char,varid,'long_name','particle number emissions of '+spec
  ncdf_attput,ncid,/char,varid,'history',var_hist
  ncdf_attput,ncid,/float,varid,'molecular_weight',mw_num
  ;Define global attributes
  ncdf_attput,ncid,/GLOBAL,/char,'data_title','Emissions of '+spec+' for CESM-MAM4 CMIP6'
  ncdf_attput,ncid,/GLOBAL,/float,'molecular_weight',mw_num
  ncdf_attput,ncid,/GLOBAL,/char,'data_creator','Jean-Francois Lamarque (lamar@ucar.edu) and Louisa Emmons (emmons@ucar.edu)'
  ncdf_attput,ncid,/GLOBAL,/char,'data_summary','Emissions of number of particles for corresponding aerosol emissions for MAM4 (bc, pom, so4). These MAM number emissions are scaled by 6.02E26 (Avog.*g/kg) to account for scaling applied in CESM/MOZART when read into model. '
  ncdf_attput,ncid,/GLOBAL,/char,'cesm_contact','Louisa Emmons or Simone Tilmes'
  ncdf_attput,ncid,/GLOBAL,/char,'creation_date',creation_date
  ncdf_attput,ncid,/GLOBAL,/char,'update_date',sdate_today
  ncdf_attput,ncid,/GLOBAL,/char,'history','Original aerosol file: '+filename
  ncdf_attput,ncid,/GLOBAL,/char,'data_script',thisfile
  ncdf_attput,ncid,/GLOBAL,/char,'data_source_url','http://www.globalchange.umd.edu/ceds/ceds-cmip6-data/'
  ncdf_attput,ncid,/GLOBAL,/char,'data_reference','Number calculations described in X.Liu et al., GMD, 2012.'
  ncdf_control,ncid,/ENDEF
  ncdf_varput,ncid,'lon',lon
  ncdf_varput,ncid,'lat',lat
  ncdf_varput,ncid,'altitude',altitude
  ncdf_varput,ncid,'altitude_int',altitude_int
  ncdf_varput,ncid,'time',time
  ncdf_varput,ncid,'date',date
  ncdf_varput,ncid,spec,emis_num
  ncdf_close,ncid
endfor

end 
