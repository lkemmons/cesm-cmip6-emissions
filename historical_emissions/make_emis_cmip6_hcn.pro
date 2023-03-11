; HCN and CH3CN emissions
; created from CO biofuel and bb emissions
; Reads files regridded to 1 deg (created with read_write_biofuel.ncl
;   and read_write_bb.nc from original CEDS files)
; Concatenate anthro and bb files for all time slices for each species
; Write separate anthro and bb files
; CEDS provides HCN BB

pro make_emis_cmip6_hcn

today = bin_date(systime())
todaystr = String(today[0:2],format='(i4,"/",i2.2,"/",i2.2)')
sdate_today = String(today[0:2],format='(i4,i2.2,i2.2)')
creation_date = '20170608'
thisfile = Routine_filepath()

path_in = '/glade/scratch/emmons/Emissions_CMIP6/concat/'
path_out = '/glade/p/cesm/chwg_dev/emmons/CMIP6_emissions_1750_2015_v20170608/'

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

;MOZART mechanism species names to create from CO
specs = ['HCN','CH3CN']
types = ['anthro'] ;'bb'

for ispec = 0,n_elements(specs)-1 do begin
  spec = specs[ispec]
  case spec of
   'HCN': begin
     mw = 27.
     sf_co = 0.006
   end
   'CH3CN': begin
      mw = 41.
      sf_co = 0.002
   end 
  endcase 

  for itype=0,0 do begin
  type = types[itype]

  if (type eq 'anthro') then begin

  ; read biofuel file
  file_a = path_in+'emissions-cmip6_CO_biofuel-anthro_surface_175001-201412_0.9x1.25_c20170615.nc'
  print,'reading ',file_a
  ncid_a = ncdf_open(file_a)
  ncdf_varget,ncid_a,'lon',lon
  ncdf_varget,ncid_a,'lat',lat
  ncdf_varget,ncid_a,'date',date_a
  ncdf_varget,ncid_a,'emiss_anthro',emis1
  ncdf_attget,ncid_a,/global,'source_file',source_file_b
  source_file_an = String(source_file_b)
  emiss_anthro = emis1*sf_co 
  ncdf_close,ncid_a
  ntima = n_elements(date_a)
  nlon = n_elements(lon)
  nlat = n_elements(lat)
 
  ;extend anthro, repeating 2014 for 2015
  emiss_anthro15 = fltarr(nlon,nlat,ntim)
  emiss_anthro15[*,*,0:ntima-1] = emiss_anthro
  emiss_anthro15[*,*,ntima:ntim-1] = emiss_anthro[*,*,ntima-12:ntima-1]

; write new anthro file
  newfile = path_out+'emissions-cmip6_'+spec+'_anthro_surface_1750-2015_0.9x1.25_c'+creation_date+'.nc'
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
  varid = ncdf_vardef(ncid,'emiss_anthro',[xid,yid,tid],/float)
  ncdf_attput,ncid,/char,'emiss_anthro','units','molecules/cm2/s'
  ncdf_attput,ncid,/char,'emiss_anthro','long_name','anthro-biofuel emissions'
  ncdf_attput,ncid,/float,'emiss_anthro','molecular_weight',mw
 label = 'CO emissions scaled by '+string(sf_co,format='(f5.3)')
  ncdf_attput,ncid,/char,'emiss_anthro','history',label
 
  ;Define global attributes
  ncdf_attput,ncid,/GLOBAL,/char,'data_title','Anthropogenic emissions of '+spec+' for CMIP6'
  ncdf_attput,ncid,/GLOBAL,/float,'molecular_weight',mw
  ncdf_attput,ncid,/GLOBAL,/char,'data_creator','Jean-Francois Lamarque (lamar@ucar.edu) and Louisa Emmons (emmons@ucar.edu)'
  ncdf_attput,ncid,/GLOBAL,/char,'data_summary','Emissions from the Community Emission Data System (CEDS) have been manipulated for use in CESM2 for CMIP6. Original file: '+source_file_an
  ncdf_attput,ncid,/GLOBAL,/char,'cesm_contact','Louisa Emmons or Simone Tilmes'
  ncdf_attput,ncid,/GLOBAL,/char,'creation_date',creation_date
  ncdf_attput,ncid,/GLOBAL,/char,'update_date',sdate_today
  ncdf_attput,ncid,/GLOBAL,/char,'history','original files have been regridded, units changed, concatenated.'
  ncdf_attput,ncid,/GLOBAL,/char,'data_script',thisfile
  ncdf_attput,ncid,/GLOBAL,/char,'data_source_url','http://www.globalchange.umd.edu/ceds/ceds-cmip6-data/'
  ncdf_attput,ncid,/GLOBAL,/char,'data_reference','Hoesly et al., GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)'
  ncdf_control,ncid,/ENDEF
  ncdf_varput,ncid,'lon',lon
  ncdf_varput,ncid,'lat',lat
  ncdf_varput,ncid,'time',time
  ncdf_varput,ncid,'date',date
  ncdf_varput,ncid,'emiss_anthro',emiss_anthro15
  ncdf_close,ncid
  endif
  if (type eq 'b') then begin
  ; read bb emis
  file_b = path_in+'emissions-cmip6_CO_biomassburning_surface_175001-201512_0.9x1.25_c20170117.nc'
  print,'reading ',file_b
  ncid_b = ncdf_open(file_b)
  ncdf_varget,ncid_b,'lon',lon
  ncdf_varget,ncid_b,'lat',lat
  ncdf_varget,ncid_b,'date',date
  ncdf_varget,ncid_b,'time',time
  ncdf_varget,ncid_b,'time_bnds',time_bnds
  ncdf_varget,ncid_b,'emiss_bb',emis1
  ncdf_attget,ncid_b,/global,'source_file',source_file_b
  source_file_bb = String(source_file_b)
  emiss_bb = emis1*sf_co
  ncdf_close,ncid_b
  ntimb = n_elements(time)
  nlon = n_elements(lon)
  nlat = n_elements(lat)
; write new bb file
 if (spec eq 'HCN') then goto,skipwrite
  newfile = path_out+'emissions-cmip6_'+spec+'_bb_surface_1750-2015_0.9x1.25_c'+creation_date+'.nc'
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
  varid = ncdf_vardef(ncid,'emiss_bb',[xid,yid,tid],/float)
  ncdf_attput,ncid,/char,'emiss_bb','units','molecules/cm2/s'
  ncdf_attput,ncid,/char,'emiss_bb','long_name','biomass burning emissions'
  ncdf_attput,ncid,/float,'emiss_bb','molecular_weight',mw
  label = 'CO emissions scaled by '+string(sf_co,format='(f5.3)')
  help,label
  ncdf_attput,ncid,/char,'emiss_bb','history',label

  ;Define global attributes
  ncdf_attput,ncid,/GLOBAL,/char,'data_title','Biomass burning missions of '+spec+' for CMIP6'
  ncdf_attput,ncid,/GLOBAL,/float,'molecular_weight',mw
  ncdf_attput,ncid,/GLOBAL,/char,'data_creator','Jean-Francois Lamarque (lamar@ucar.edu) and Louisa Emmons (emmons@ucar.edu)'
  ncdf_attput,ncid,/GLOBAL,/char,'data_summary','Emissions from the Community Emission Data System (CEDS) have been manipulated for use in CESM2 for CMIP6. Original file: '+source_file_bb
  ncdf_attput,ncid,/GLOBAL,/char,'cesm_contact','Louisa Emmons or Simone Tilmes'
  ncdf_attput,ncid,/GLOBAL,/char,'creation_date',creation_date
  ncdf_attput,ncid,/GLOBAL,/char,'update_date',sdate_today
  ncdf_attput,ncid,/GLOBAL,/char,'history','original files have been regridded, units changed, concatenated.'
  ncdf_attput,ncid,/GLOBAL,/char,'data_script',thisfile
  ncdf_attput,ncid,/GLOBAL,/char,'data_source_url','http://www.globalchange.umd.edu/ceds/ceds-cmip6-data/'
  ncdf_attput,ncid,/GLOBAL,/char,'data_reference','van Marle et al., GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-32/)'

  ncdf_control,ncid,/ENDEF

  ncdf_varput,ncid,'lon',lon
  ncdf_varput,ncid,'lat',lat
  ncdf_varput,ncid,'time',time
  ncdf_varput,ncid,'date',date
  ncdf_varput,ncid,'emiss_bb',emiss_bb
  ncdf_close,ncid
skipwrite:
endif
endfor
endfor

end

