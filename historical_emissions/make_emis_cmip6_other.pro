;
; Make emissions files with ocean and soil emissions
; Reads CCMI soil NOx, ocean emissions, etc.

pro make_emis_cmip6_other

today = bin_date(systime())
todaystr = String(today[0:2],format='(i4,"/",i2.2,"/",i2.2)')
sdate_today = String(today[0:2],format='(i4,i2.2,i2.2)')
creation_date = '20170322'
thisfile = Routine_filepath()

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
nyrs = ntim/12

;species with soil or ocean emissions
specs = ['CO','NO','NH3', 'C2H6','C3H8','C2H4','C3H6','DMS']

for ispec = 0,n_elements(specs)-1 do begin
  spec = specs[ispec]
  print,spec
 
  varname_s = ''
  varname_o = ''

  case spec of 
  'NO': begin
    varname_s = 'emiss_soils'
  end
  'CO': begin
    varname_o = 'emiss_oceans'
   end
  'DMS': begin
    varname_o = 'emiss_oceans'
   end
  'NH3': begin
    varname_s = 'emiss_soils'
    varname_o = 'emiss_oceans'
   end
  'C2H6': begin
    varname_o = 'emiss_oceans'
   end
  'C2H4': begin
    varname_o = 'emiss_oceans'
   end
  'C3H8': begin
    varname_o = 'emiss_oceans'
   end
  'C3H6': begin
    varname_o = 'emiss_oceans'
   end
 endcase

  ; read CCMI for other types (besides biogenic) - ocean, soils
   path1 = '/glade/p/work/emmons/emis/CCMI/'
   file1 = File_search(path1+'ccmi*_1yr_'+spec+'_0.9x1.25_mol.nc')
   file1 = file1[0]
   source_file = file1
   print,file1
   ncid_c = ncdf_open(file1)
   ncdf_varget,ncid_c,'lon',lon
   ncdf_varget,ncid_c,'lat',lat
   ncdf_attget,ncid_c,/global,'molecular_weight',mw
   nlon = n_elements(lon)
   nlat = n_elements(lat)

  newfile = path_out+'emissions-cmip6_'+spec+'_other_surface_1750-2015_0.9x1.25_c'+creation_date+'.nc'  
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
  ncdf_attput, ncid, tvarid,/char, 'units', 'days since 1750-01-01 00:00:00'
  ncdf_attput, ncid, tvarid,/char, 'calendar', 'Gregorian'
  tvarid = ncdf_vardef(ncid,'date',[tid],/long)
  ncdf_attput, ncid, tvarid,/char, 'units', 'YYYYMMDD'
  ncdf_attput, ncid, tvarid,/char, 'long_name', 'Date'

  if (strlen(varname_s) gt 1) then begin
   varid = ncdf_vardef(ncid,varname_s,[xid,yid,tid],/float)
   ncdf_attput,ncid,/char,varname_s,'units','molecules/cm2/s'
   ncdf_attput,ncid,/char,varname_s,'long_name',spec+' soil emissions'
   ncdf_attput,ncid,/float,varname_s,'molecular_weight',mw
  endif
  if (strlen(varname_o) gt 1) then begin
   varid = ncdf_vardef(ncid,varname_o,[xid,yid,tid],/float)
   ncdf_attput,ncid,/char,varname_o,'units','molecules/cm2/s'
   ncdf_attput,ncid,/char,varname_o,'long_name',spec+' ocean emissions'
   ncdf_attput,ncid,/float,varname_o,'molecular_weight',mw
  endif
  
  ;Define global attributes
  ncdf_attput,ncid,/GLOBAL,/char,'data_title','Natural emissions of '+spec+' from CCMI for CMIP6'
  ncdf_attput,ncid,/GLOBAL,/float,'molecular_weight',mw
  ncdf_attput,ncid,/GLOBAL,/char,'data_creator','Jean-Francois Lamarque (lamar@ucar.edu) and Louisa Emmons (emmons@ucar.edu)'
  ncdf_attput,ncid,/GLOBAL,/char,'data_summary','Emissions from CCMI have been manipulated for use in CESM2 for CMIP6. '+source_file
  ncdf_attput,ncid,/GLOBAL,/char,'cesm_contact','Louisa Emmons or Simone Tilmes'
  ncdf_attput,ncid,/GLOBAL,/char,'creation_date',creation_date
  ncdf_attput,ncid,/GLOBAL,/char,'update_date',sdate_today
  ncdf_attput,ncid,/GLOBAL,/char,'history','A single year of emissions (monthly variation) is repeated for all years.'
  ncdf_attput,ncid,/GLOBAL,/char,'data_script',thisfile
  ncdf_attput,ncid,/GLOBAL,/char,'data_source_url','/glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/ccmi_1960-2008/'
  ncdf_attput,ncid,/GLOBAL,/char,'data_reference','Used in CESM1 for CCMI, e.g., Tilmes et al., GMD, 2016.'

  ncdf_control,ncid,/ENDEF

  ncdf_varput,ncid,'lon',lon
  ncdf_varput,ncid,'lat',lat
  ncdf_varput,ncid,'time',time
  ncdf_varput,ncid,'date',date

  if (strlen(varname_s) gt 1) then begin
   ncdf_varget,ncid_c,varname_s,emis1
   emis_yrs = fltarr(nlon,nlat,ntim)
   for iyr = 0,nyrs-1 do emis_yrs[*,*,(iyr*12):(iyr*12+11)] = emis1
   ncdf_varput,ncid,varname_s,emis_yrs
  endif
  if (strlen(varname_o) gt 1) then begin
   if (spec eq 'DMS') then ncdf_varget,ncid_c,'ocean',emis1 else $
   ncdf_varget,ncid_c,varname_o,emis1
   emis_yrs = fltarr(nlon,nlat,ntim)
   for iyr = 0,nyrs-1 do emis_yrs[*,*,(iyr*12):(iyr*12+11)] = emis1
   ncdf_varput,ncid,varname_o,emis_yrs
  endif

  ncdf_close,ncid
  ncdf_close,ncid_c

endfor

end

