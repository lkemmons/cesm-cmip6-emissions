; SOAG emissions for MAM4
; SOA precursor of lumped HCs
; mass yields from Liu et al., GMD, 2012
; 5% BIGALK, 5% BIGENE, 15% aromatics, 4% isoprene, 25% monoterpenes

pro make_emis_cmip6_soag

today = bin_date(systime())
todaystr = String(today[0:2],format='(i4,"/",i2.2,"/",i2.2)')
sdate_today = String(today[0:2],format='(i4,i2.2,i2.2)')
creation_date = '20170608'
thisfile = Routine_filepath()

;path_in = '/glade/p/cesm/chwg_dev/emmons/CMIP6_emissions_1750_2015_v20170322/'
path_in = '/glade/p/cesm/chwg_dev/emmons/CMIP6_emissions_1750_2015_v20170608/'
path_biog = '/glade/p/cesm/chwg_dev/emmons/CMIP6_emissions_1750_2015/biogenic_climo/'
path_out = '/glade/p/cesm/chwg_dev/emmons/CMIP6_emissions_1750_2015_v20170322/'
path_out = '/glade/p/cesm/chwg_dev/emmons/CMIP6_emissions_1750_2015_v20170608/'

spec = 'SOAG'
mw_soag = 12.
specs_in = ['BIGALK','BIGENE','TOLUENE','BENZENE','XYLENES','ISOP','MTERP']
yields = [0.05, 0.05, 0.15, 0.15, 0.15, 0.04, 0.25]

;get dims
file0 = path_in+'emissions-cmip6_BIGALK_anthro_surface_1750-2015_0.9x1.25_c'+creation_date+'.nc'
ncid = ncdf_open(file0)
ncdf_varget,ncid,'lon',lon
ncdf_varget,ncid,'lat',lat
ncdf_varget,ncid,'time',time
ncdf_varget,ncid,'date',date
ncdf_close,ncid
nlon = n_elements(lon)
nlat = n_elements(lat)
ntim = n_elements(date)

types = ['anthro','bb','biogenic']

for itype = 0,2 do begin
 type = types[itype]
 print,type
 varname = 'emiss_'+type
 emis_soag = fltarr(nlon,nlat,ntim)  ;sum all species
 hist = 'SOAG='
 for isp = 0,n_elements(specs_in)-1 do begin
  specin = specs_in[isp]
  ; read cmip6 file
  if (type eq 'biogenic') then path_a = path_biog else path_a = path_in
  file_a = File_search(path_a+'emissions-cmip6_'+specin+'_'+type+'_*.nc',count=nf)
  if (nf eq 1) then begin
    print,'reading ',file_a
    ncid_a = ncdf_open(file_a)
    ncdf_varget,ncid_a,varname,emis1
    ncdf_attget,ncid_a,/global,'molecular_weight',mw
    emis_soag = emis_soag + emis1*yields[isp] *mw /mw_soag
    hist=hist + string(yields[isp],format='("+",f4.2,"*")')+specin
    ncdf_close,ncid_a
  endif else if (nf gt 1) then print,'more than 1 anthro file' 
 endfor
print,hist

; write new SOAG file
  newfile = path_out+'emissions-cmip6_'+spec+'_'+type+'_surface_1750-2015_0.9x1.25_c'+creation_date+'.nc'
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
  varid = ncdf_vardef(ncid,varname,[xid,yid,tid],/float)
  ncdf_attput,ncid,/char,varid,'units','molecules/cm2/s'
  ncdf_attput,ncid,/char,varid,'long_name',type+' emissions'
  ncdf_attput,ncid,/float,varid,'molecular_weight',mw_soag
  ncdf_attput,ncid,/char,varid,'history',hist
 
  ;Define global attributes
  ncdf_attput,ncid,/GLOBAL,/char,'data_title','Emissions of '+spec+' for CMIP6'
  ncdf_attput,ncid,/GLOBAL,/float,'molecular_weight',mw_soag
  ncdf_attput,ncid,/GLOBAL,/char,'data_creator','Louisa Emmons (emmons@ucar.edu)'
  ncdf_attput,ncid,/GLOBAL,/char,'data_summary','Lumped HCs precursor of SOA for simple scheme in MAM, based on various fractions of alkanes, aromatics and biogenic HCs.'
  ncdf_attput,ncid,/GLOBAL,/char,'cesm_contact','Louisa Emmons or Simone Tilmes'
  ncdf_attput,ncid,/GLOBAL,/char,'creation_date',creation_date
  ncdf_attput,ncid,/GLOBAL,/char,'update_date',sdate_today
  ncdf_attput,ncid,/GLOBAL,/char,'history',' '
  ncdf_attput,ncid,/GLOBAL,/char,'data_script',thisfile
  ncdf_attput,ncid,/GLOBAL,/char,'data_source_url',' '
  ncdf_attput,ncid,/GLOBAL,/char,'data_reference',' '

  ncdf_control,ncid,/ENDEF

  ncdf_varput,ncid,'lon',lon
  ncdf_varput,ncid,'lat',lat
  ncdf_varput,ncid,'time',time
  ncdf_varput,ncid,'date',date
  ncdf_varput,ncid,varname,emis_soag
  ncdf_close,ncid

; also write files with SOAG*1.5
  newfile = path_out+'emissions-cmip6_'+spec+'x1.5_'+type+'_surface_1750-2015_0.9x1.25_c'+creation_date+'.nc'
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
  varid = ncdf_vardef(ncid,varname,[xid,yid,tid],/float)
  ncdf_attput,ncid,/char,varid,'units','molecules/cm2/s'
  ncdf_attput,ncid,/char,varid,'long_name',type+' emissions'
  ncdf_attput,ncid,/float,varid,'molecular_weight',mw_soag
  ncdf_attput,ncid,/char,varid,'history','1.5*'+hist

 ;Define global attributes
  ncdf_attput,ncid,/GLOBAL,/char,'data_title','Emissions of '+spec+' for CMIP6, scaled by 1.5'
  ncdf_attput,ncid,/GLOBAL,/float,'molecular_weight',mw_soag
  ncdf_attput,ncid,/GLOBAL,/char,'data_creator','Louisa Emmons (emmons@ucar.edu)'
  ncdf_attput,ncid,/GLOBAL,/char,'data_summary','Lumped HCs precursor of SOA for simple scheme in MAM, based on various fractions of alkanes, aromatics and biogenic HCs. Scaled by 1.5.'
  ncdf_attput,ncid,/GLOBAL,/char,'cesm_contact','Louisa Emmons or Simone Tilmes'
  ncdf_attput,ncid,/GLOBAL,/char,'creation_date',creation_date
  ncdf_attput,ncid,/GLOBAL,/char,'update_date',sdate_today
  ncdf_attput,ncid,/GLOBAL,/char,'history',' '
  ncdf_attput,ncid,/GLOBAL,/char,'data_script',thisfile
  ncdf_attput,ncid,/GLOBAL,/char,'data_source_url',' '
  ncdf_attput,ncid,/GLOBAL,/char,'data_reference',' '
  ncdf_control,ncid,/ENDEF

  ncdf_varput,ncid,'lon',lon
  ncdf_varput,ncid,'lat',lat
  ncdf_varput,ncid,'time',time
  ncdf_varput,ncid,'date',date
  emis15 = emis_soag*(1.5)
  ncdf_varput,ncid,varname,emis15
  ncdf_close,ncid

 endfor
end
