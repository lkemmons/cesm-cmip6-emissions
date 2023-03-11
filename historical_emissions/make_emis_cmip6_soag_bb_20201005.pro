; SOAG emissions for MAM4
; SOA precursor of lumped HCs
; mass yields from Liu et al., GMD, 2012
; 5% BIGALK, 5% BIGENE, 15% aromatics, 4% isoprene, 25% monoterpenes

pro make_emis_cmip6_soag_bb_20201005

today = bin_date(systime())
todaystr = String(today[0:2],format='(i4,"/",i2.2,"/",i2.2)')
sdate_today = String(today[0:2],format='(i4,i2.2,i2.2)')
creation_date = sdate_today
thisfile = Routine_filepath()

;path_in = '/glade/p/cesm/chwg_dev/emmons/CMIP6_emissions_1750_2015_v20170322/'
;date_in = '20170322'
;path_in = '/glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/emissions_ssp370/'
;date_in = '20190222'
;path_in = '/glade/p/cesm/chwg_dev/emmons/CMIP6_emissions_1750_2015_FINAL/'
;scenario = 'historical'
;date_range = '1750-2015'

;path_in = '/glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/emissions_ssp370/'
;scenario = 'ScenarioMIP_IAMC-AIM-ssp370-1-1'
;path_in = '/glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/emissions_ssp126/'
;scenario = 'ScenarioMIP_IAMC-IMAGE-ssp126-1-1'
;path_in = '/glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/emissions_ssp245/'
;scenario = 'ScenarioMIP_IAMC-MESSAGE-GLOBIOM-ssp245-1-1'
;path_in = '/glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/emissions_ssp370-lowNTCF/'
;scenario = 'AerChemMIP_IAMC-AIM-ssp370-lowNTCF-1-1'
;path_in = '/glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/emissions_ssp534_over/'
;scenario = 'ScenarioMIP_IAMC-REMIND-MAGPIE-ssp534-over-1-1'
;path_in = '/glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/emissions_ssp585/'
;scenario = 'ScenarioMIP_IAMC-REMIND-MAGPIE-ssp585-1-1'
;date_range = '175001-210101'
;resol = '1.9x2.5'
;path_in = '/glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/CMIP6_emissions_1750_2015_2deg/'
;scenario = ''
resol = '0.9x1.25'
scenario = ''
;date_range = '1750-2015'
;path_in = '/glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/CMIP6_emissions_1750_2015/'
path_in = '/glade/p/cesmdata/cseg/inputdata/atm/cam/chem/emis/CMIP6_emissions_2000climo/'
date_range = '2000climo'

path_out = '/glade/p/cesm/chwg_dev/emmons/SOAG_corrected/'

spec = 'SOAGx1.5'
mw_soag = 12.
specs_in = ['BIGALK','BIGENE','TOLUENE','BENZENE','XYLENES','ISOP','MTERP']
yields = [0.05, 0.05, 0.15, 0.15, 0.15, 0.04, 0.25]

;get dims
file0 = File_search(path_in+'emissions-cmip*_'+specs_in[0]+'_bb_surface*'+date_range+'_*.nc')
print,file0
ncid = ncdf_open(file0[0])
ncdf_varget,ncid,'lon',lon
ncdf_varget,ncid,'lat',lat
ncdf_varget,ncid,'time',time
ncdf_varget,ncid,'date',date
ncdf_close,ncid
nlon = n_elements(lon)
nlat = n_elements(lat)
ntim = n_elements(date)

type = 'bb'
print,type
 varname = 'emiss_'+type
 emis_soag = fltarr(nlon,nlat,ntim)  ;sum all species
 hist = 'SOAG=1.5*('
 hist_files = 'files used:'
 for isp = 0,n_elements(specs_in)-1 do begin
  specin = specs_in[isp]
  ; read cmip6 file
  file_a = File_search(path_in+'emissions-cmip*_'+specin+'_bb_surface*'+date_range+'_*.nc',count=nf)
  ifile=0
  if (nf gt 1) then begin
    for i=0,nf-1 do print,i,' ',file_a[i]
    read,ifile,prompt='pick file'
  endif
  file1 = file_a[ifile]
  hist_files = hist_files+' '+file1
  print,'reading ',file1
  ncid_a = ncdf_open(file1)
  ncdf_varget,ncid_a,varname,emis1
  ncdf_attget,ncid_a,/global,'molecular_weight',mw
  emis_soag = emis_soag + emis1*yields[isp] *mw /mw_soag
  hist=hist + string(yields[isp],mw,mw_soag,format='("+",f4.2,"*",i0,"/",i0,"*")')+specin
  ncdf_close,ncid_a
 endfor
 hist = hist+')'
help,hist
help,hist_files

; write new SOAG file
  if (strlen(scenario) gt 1) then $
  newfile = path_out+'emissions-cmip6-'+scenario+'_'+spec+'_bb_surface_'+date_range+'_0.9x1.25_c'+creation_date+'.nc' $
  else newfile = path_out+'emissions-cmip6_'+spec+'_bb_surface_'+date_range+'_'+resol+'_c'+creation_date+'.nc'

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
  ncdf_attput,ncid,/GLOBAL,/char,'history',hist_files
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

end
