;
; Make emissions files with biogenic emissions for 2000 (from CCMI)
;  all years the same
;
pro make_emis_cmip6_biogenic

today = bin_date(systime())
todaystr = String(today[0:2],format='(i4,"/",i2.2,"/",i2.2)')
sdate_today = String(today[0:2],format='(i4,i2.2,i2.2)')

path_out = '/glade/p/cesm/chwg_dev/emmons/CMIP6_emissions_1750_2015/'

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

specs = ['ISOPRENE','TERPENES','CO','CH3OH','CH3COCH3','CH3CHO','CH2O']
varname = 'emiss_bio'

for ispec = 0,n_elements(specs)-1 do begin
  spec = specs[ispec]
  print,spec
 
  ; read CCMI biogenics 2000 file
   path1 = '/glade/p/work/emmons/emis/CCMI/biogenic/'
   file1 = File_search(path1+'ccmi_biogenic_2008_'+spec+'_0.9x1.25_mol.nc')
   file1 = file1[0]
   print,file1
   ncid_c = ncdf_open(file1)
   ncdf_varget,ncid_c,'lon',lon
   ncdf_varget,ncid_c,'lat',lat
   ncdf_attget,ncid_c,varname,'molecular_weight',mw
print,mw
   nlon = n_elements(lon)
   nlat = n_elements(lat)

  newfile = path_out+'emissions-cmip6_'+spec+'_biogenic_surface_1750-2015_0.9x1.25_c'+sdate_today+'.nc'  
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

  varid = ncdf_vardef(ncid,'emiss_biogenic',[xid,yid,tid],/float)
  ncdf_attput,ncid,/char,varid,'units','molecules/cm2/s'
  ncdf_attput,ncid,/char,varid,'long_name',spec+' MEGANv2.12008  biogenic emissions'

  ;Define global attributes
  ncdf_attput,ncid,/GLOBAL,/char,'title','Emissions of '+spec+' from CCMI for CMIP6'
  ncdf_attput,ncid,/GLOBAL,/char,'author','Jean-Francois Lamarque (lamar@ucar.edu) and Louisa Emmons (emmons@ucar.edu)'
  ncdf_attput,ncid,/GLOBAL,/char,'creation_date',sdate_today
  ncdf_attput,ncid,/GLOBAL,/float,'molecular_weight',mw
  ncdf_attput,ncid,/GLOBAL,/char,'source_file',file1

  ncdf_control,ncid,/ENDEF

  ncdf_varput,ncid,'lon',lon
  ncdf_varput,ncid,'lat',lat
  ncdf_varput,ncid,'time',time
  ncdf_varput,ncid,'date',date

  ;repeat 1 yr emissions for all years
  ncdf_varget,ncid_c,varname,emis1
  emis_yrs = fltarr(nlon,nlat,ntim)
  for iyr = 0,nyrs-1 do emis_yrs[*,*,(iyr*12):(iyr*12+11)] = emis1
  ncdf_varput,ncid,'emiss_biogenic',emis_yrs

  ncdf_close,ncid
  ncdf_close,ncid_c

endfor

end

