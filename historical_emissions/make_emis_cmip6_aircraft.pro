pro make_emis_cmip6_aircraft

today = bin_date(systime())
todaystr = String(today[0:2],format='(i4,"/",i2.2,"/",i2.2)')
sdate_today = String(today[0:2],format='(i4,i2.2,i2.2)')
creation_date = '20170608'
thisfile = Routine_filepath()

altitude_int = [0, 0.61, 1.22, 1.83, 2.44, 3.05, 3.66, 4.27, 4.88, 5.49, 6.1, $
    6.71, 7.32, 7.93, 8.54, 9.150001, 9.76, 10.37, 10.98, 11.59, 12.2, 12.81, $ 
    13.42, 14.03, 14.64, 15.25]
altitude = [0.305, 0.915, 1.525, 2.135, 2.745, 3.355, 3.965, 4.575, 5.185, $
    5.795, 6.405, 7.015, 7.625, 8.235001, 8.845, 9.455001, 10.065, 10.675, $ 
    11.285, 11.895, 12.505, 13.115, 13.725, 14.335, 14.945]

specs = ['NOx', 'BC', 'SO2'] 
spec_out = ['NO2', 'bc_a4', 'SO2'] 
;path_new = '/glade/p/cesm/chwg_dev/emmons/CMIP6_emissions_1750_2015_v20170322/'
path_new = '/glade/p/cesm/chwg_dev/emmons/CMIP6_emissions_1750_2015_v20170608/'

;v2 emissions have aircraft emissions starting in 1850
; zero them until 1920
; Include 1750 and 1850 and 1919 with 0
  yr1emis=1919
  ntim = (2015-yr1emis+3)*12
  time = fltarr(ntim)
  date = lonarr(ntim)
  for itim=0,11 do begin
    mon = itim+1
    yr = 1750
    date[itim] = yr*10000L + mon*100L + 16
    time[itim] = julday(mon,16,yr)-julday(1,1,1750)
    print,itim,date[itim],time[itim]
  endfor
  for itim=12,23 do begin
    mon = itim-11
    yr = 1850
    date[itim] = yr*10000L + mon*100L + 16
    time[itim] = julday(mon,16,yr)-julday(1,1,1750)
    print,itim,date[itim],time[itim]
  endfor
  for yr = yr1emis,2015 do begin
    for mon = 1,12 do begin
     itim = 24 + (mon-1) + (yr-yr1emis)*12
     date[itim] = yr*10000L + mon*100L + 16
     time[itim] = julday(mon,16,yr) - julday(1,1,1750)
     if (mon eq 1) then print,itim,date[itim],time[itim]
    endfor
  endfor

for ispec = 0,0 do begin

  spec = specs[ispec]

  ;read concatenated CEDS aircraft file
  file1 = '/glade/scratch/emmons/Emissions_CMIP6/concat/emissions-cmip6_'+spec+'_aircraft_190001-201412_0.9x1.25_c20170608.nc'
  ncid = ncdf_open(file1)
  ncdf_varget,ncid,'lon',lon
  ncdf_varget,ncid,'lat',lat
  ncdf_varget,ncid,'level',lev
  ncdf_varget,ncid,'date',date_a
  ncdf_varget,ncid,'emiss_aircraft',emiss_air_cm2  ;molecules/cm2/s
  ncdf_attget,ncid,'emiss_aircraft','molecular_weight',mw
  ncdf_attget,ncid,/global,'source_file',source_file_b
  source_file = String(source_file_b)
  ncdf_close,ncid
  nlon = n_elements(lon)
  nlat = n_elements(lat)
  nlev = n_elements(lev)
  ntima = n_elements(date_a)
help,mw
print,lev,altitude
  ;convert to molecules/cm3/s
  emiss_air = fltarr(nlon,nlat,nlev,ntim)
  for ilev=0,nlev-1 do begin
    dz = (altitude_int[ilev+1] - altitude_int[ilev])*1.e5   ;km -> cm
    ;prior to 1920 just zeros
    i1 = where(date eq 19200116)
    i2 = where(date eq 20141216)
    i1a = where(date_a eq 19200116) 
    emiss_air[*,*,ilev,i1:i2] = emiss_air_cm2[*,*,ilev,i1a:ntima-1] /dz  ;molecules/cm3/s
    ;repeat 2014 for 2015
    emiss_air[*,*,ilev,(ntim-12):(ntim-1)] = emiss_air[*,*,ilev,(ntim-24):(ntim-13)] 
  endfor

  newfile = path_new+'emissions-cmip6_'+spec_out[ispec]+'_aircraft_vertical_1750-2015_0.9x1.25_c'+creation_date+'.nc'
  print,newfile
  ncid = ncdf_create(newfile,/clobber)
  xid = ncdf_dimdef(ncid,'lon',nlon)
  yid = ncdf_dimdef(ncid,'lat',nlat)
  zid = ncdf_dimdef(ncid,'altitude',nlev)
  zbid = ncdf_dimdef(ncid,'altitude_int',nlev+1)
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

  varid = ncdf_vardef(ncid,'emiss_aircraft',[xid,yid,zid,tid],/float)
  ncdf_attput,ncid,/char,'emiss_aircraft','units','molecules/cm3/s'
  ncdf_attput,ncid,/char,'emiss_aircraft','long_name',spec+' CEDS aircraft emissions'
  ncdf_attput,ncid,/float,'emiss_aircraft','molecular_weight',mw

  ;Define global attributes
  ncdf_attput,ncid,/GLOBAL,/char,'data_title','Aircraft emissions of '+spec+' for CMIP6'
  ncdf_attput,ncid,/GLOBAL,/float,'molecular_weight',mw
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
  ncdf_varput,ncid,'emiss_aircraft',emiss_air
  ncdf_close,ncid

help,date,emiss_air

endfor

end

