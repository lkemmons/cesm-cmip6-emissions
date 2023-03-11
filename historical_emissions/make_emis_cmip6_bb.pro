;
; Create biomass burning emissions file for 1750-2015 from CEDS/CMIP6
; Reads files regridded to 1 deg (created with read_write_bb.nc)
; bb files 1750-2015
; Rename to MOZART VOC species

pro make_emis_cmip6_bb

today = bin_date(systime())
todaystr = String(today[0:2],format='(i4,"/",i2.2,"/",i2.2)')
sdate_today = String(today[0:2],format='(i4,i2.2,i2.2)')
creation_date = '20170322'
thisfile = Routine_filepath()

path_concat = '/glade/scratch/emmons/Emissions_CMIP6/concat/'
path_out = '/glade/p/cesm/chwg_dev/emmons/CMIP6_emissions_1750_2015_v20170322/'

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

;MOZART mechanism species names
specs = ['CO','NO','NH3', $   ;'SO2',
  'C2H6','C3H8','BIGALK','C2H4','C3H6','C2H2','BIGENE', $
  'BENZENE','TOLUENE','XYLENES','CH2O','CH3CHO','CH3OH','C2H5OH', $
  'CH3COCH3','MEK','HCOOH','CH3COOH','ISOP','MTERP','DMS','CH3COCHO','GLYALD','HCN', $
  'BC','OC','CB1','CB2','OC1','OC2','bc_a4','pom_a4']
specs = ['pom_a4']
for ispec = 0,n_elements(specs)-1 do begin
  spec = specs[ispec]
  spec_ceds = spec
  sf = 1.0
  mw = 0.
  case spec of
   'NO': begin
      spec_ceds = 'NOx'
      mw = 30.
   end
   'C2H6': spec_ceds = 'NMVOC_C2H6'
   'C3H8': spec_ceds = 'NMVOC_C3H8'
   'BIGALK': spec_ceds = 'NMVOC_Higher_Alkanes'
   'C2H4': spec_ceds = 'NMVOC_C2H4'
   'C3H6': spec_ceds  = 'NMVOC_C3H6'
   'C2H2': spec_ceds  = 'NMVOC_C2H2'
   'BIGENE': spec_ceds  = 'NMVOC_Higher_Alkenes'
   'BENZENE': spec_ceds = 'NMVOC_C6H6'
   'TOLUENE': spec_ceds  = 'NMVOC_C7H8'
   'XYLENES': spec_ceds  = 'NMVOC_C8H10'
   'CH2O': spec_ceds  = 'NMVOC_CH2O'
   'CH3CHO': spec_ceds  = 'NMVOC_C2H4O'
   'CH3OH': spec_ceds  = 'NMVOC_CH3OH'
   'C2H5OH': spec_ceds  = 'NMVOC_C2H5OH'
   'CH3COCH3': spec_ceds  = 'NMVOC_C3H6O'
   'MEK': spec_ceds  = 'NMVOC_MEK'
   'HCOOH': spec_ceds  = 'NMVOC_HCOOH'
   'CH3COOH': spec_ceds  = 'NMVOC_CH3COOH'
   'ISOP': spec_ceds  = 'NMVOC_C5H8'
   'MTERP':  spec_ceds  = 'NMVOC_C10H16'
   'DMS': spec_ceds  = 'NMVOC_C2H6S'
   'CH3COCHO': spec_ceds  = 'NMVOC_CH3COCHO'
   'GLYALD': spec_ceds  = 'NMVOC_HOCH2CHO'
   'HCN': spec_ceds  = 'NMVOC_HCN'
   'CB1': begin
    spec_ceds = 'BC'
    sf = 0.8
    mw = 12.
   end
   'CB2': begin
    spec_ceds = 'BC'
    sf = 0.2
    mw = 12.
   end
   'bc_a4': begin
    spec_ceds = 'BC'
    sf = 1.0
    mw = 12.
   end
   'OC1': begin
    spec_ceds = 'OC'
    sf = 0.5
    mw = 12.
   end
   'OC2': begin
    spec_ceds = 'OC'
    sf = 0.5
    mw = 12.
   end
   'pom_a4': begin
    spec_ceds = 'OC'
    sf = 1.4
    mw = 12.
   end
   else: begin
     spec_ceds = spec
     sf = 1.0
   end
  endcase 

  ; read bb emis
  file_b = path_concat+'emissions-cmip6_'+spec_ceds+'_biomassburning_surface_175001-201512_0.9x1.25_c20170117.nc'
  print,'reading ',file_b
  ncid_b = ncdf_open(file_b)
  ncdf_varget,ncid_b,'lon',lon
  ncdf_varget,ncid_b,'lat',lat
  ncdf_varget,ncid_b,'date',date_b
  ncdf_varget,ncid_b,'time',time_b
  ncdf_varget,ncid_b,'emiss_bb',emiss_bb
  ncdf_attget,ncid_b,'emiss_bb','molecular_weight',mw
  ncdf_attget,ncid_b,'emiss_bb','long_name',var_origname
  ncdf_attget,ncid_b,/global,'source_file',source_file_b
  source_file = String(source_file_b)
  ncdf_close,ncid_b
  nlon = n_elements(lon)
  nlat = n_elements(lat)
  emiss_bb = emiss_bb*sf

; write new file
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
  ncdf_attput, ncid, tvarid,/char, 'long_name', 'Time'
  ncdf_attput, ncid, tvarid,/char, 'units', 'days since 1750-01-01 00:00:00'
  ncdf_attput, ncid, tvarid,/char, 'calendar', 'Gregorian'
  tvarid = ncdf_vardef(ncid,'date',[tid],/long)
  ncdf_attput, ncid, tvarid,/char, 'units', 'YYYYMMDD'
  ncdf_attput, ncid, tvarid,/char, 'long_name', 'Date'
 
  varid = ncdf_vardef(ncid,'emiss_bb',[xid,yid,tid],/float)
  ncdf_attput,ncid,/char,'emiss_bb','units','molecules/cm2/s'
  ncdf_attput,ncid,/char,'emiss_bb','long_name',spec+' CEDS biomass burning emissions'
  ncdf_attput,ncid,/char,'emiss_bb','history','CEDS species: '+spec_ceds+String(sf,format='("*",f3.1)')
  ncdf_attput,ncid,/float,'emiss_bb','molecular_weight',mw

  ;Define global attributes
  ncdf_attput,ncid,/GLOBAL,/char,'data_title','Biomass burning emissions of '+spec+' for CMIP6'
  ncdf_attput,ncid,/GLOBAL,/float,'molecular_weight',mw
  ncdf_attput,ncid,/GLOBAL,/char,'data_creator','Jean-Francois Lamarque (lamar@ucar.edu) and Louisa Emmons (emmons@ucar.edu)'
  ncdf_attput,ncid,/GLOBAL,/char,'data_summary','Emissions from the Community Emission Data System (CEDS) have been manipulated for use in CESM2 for CMIP6. Original file: '+source_file
  ncdf_attput,ncid,/GLOBAL,/char,'cesm_contact','Louisa Emmons or Simone Tilmes'
  ncdf_attput,ncid,/GLOBAL,/char,'creation_date',creation_date
  ncdf_attput,ncid,/GLOBAL,/char,'update_date',sdate_today
  ncdf_attput,ncid,/GLOBAL,/char,'history','original files have been regridded, units changed, concatenated, some species combined or renamed.'
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

endfor

end

