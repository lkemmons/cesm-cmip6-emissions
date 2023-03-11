;Create SVOC and IVOC emissions for VBS-SOA
;
; "As suggested by Jathar et al. [2014, Table 1], these precursor
;species were emitted as 0.6 x POA emissions for the IVOC fraction
;considered as lost by evaporation, and as 0.2 x NMVOC emissions for
;the unspeciated IVOC fraction of organic carbon mass. The
;corresponding SOA yields (Table 2) are derived from the GECKO-A model
;(Generator of Explicit Chemistry and Kinetics of Organics in the
;Atmosphere, [Aumont et al., 2005]) for low and high NOx conditions
;considering a mixture of n-alkane species shown in Table 3." 
;
; IVOC and SVOC emissions
;IVOC mass emissions =
;0.2 * sum of mass emissions of following HCs
;C3H6, C3H8, C2H6, C2H4, BIGENE, BIGALK
;CH3COCH3, MEK, CH3CHO, CH2O
;BENZENE, TOLUENE, XYLENES
;SVOC mass emissions = 
;0.6 * sum of mass of emissions of hydrophilic and hydrophobic POA
;You need to be careful to convert OC if emissions are in OC to OA (we used 1.4 factor).
;IVOC_emissions_total = IVOC mass emissions + (SVOC mass emissions)*1.25
;And these are only anthropogenic sources, see Jathar et al. [2014, Table 1]

pro make_emis_cmip6_svocivoc

today = bin_date(systime())
todaystr = String(today[0:2],format='(i4,"/",i2.2,"/",i2.2)')
sdate_today = String(today[0:2],format='(i4,i2.2,i2.2)')
creation_date = '20170608'
thisfile = Routine_filepath()

path_in = '/glade/p/cesm/chwg_dev/emmons/CMIP6_emissions_1750_2015_v20170608/'
path_out = '/glade/p/cesm/chwg_dev/emmons/CMIP6_emissions_1750_2015_v20170608/'

specs_ivoc = ['C3H6', 'C3H8', 'C2H6', 'C2H4', 'BIGENE', 'BIGALK', $
 'CH3COCH3', 'MEK', 'CH3CHO', 'CH2O', 'BENZENE', 'TOLUENE', 'XYLENES']

;get dims
;file0 = path_in+'emissions-cmip6_BIGALK_anthro_surface_1750-2015_0.9x1.25_c20170322.nc'
file0 = path_in+'emissions-cmip6_BIGALK_anthro_surface_1750-2015_0.9x1.25_c20170608.nc'
ncid = ncdf_open(file0)
ncdf_varget,ncid,'lon',lon
ncdf_varget,ncid,'lat',lat
ncdf_varget,ncid,'time',time
ncdf_varget,ncid,'date',date
ncdf_close,ncid
nlon = n_elements(lon)
nlat = n_elements(lat)
ntim = n_elements(date)

types = ['anthro'] ;,'bb']
ntypes = 1

mw_ivoc = 184.
mw_svoc = 310.

for itype = 0,ntypes-1 do begin
 type = types[itype]

 ;IVOC = 0.2*(HCs)
 emis_ivoc = fltarr(nlon,nlat,ntim)
 for ispec = 0,n_elements(specs_ivoc)-1 do begin
   spec_hc = specs_ivoc[ispec]
   ;file_hc = path_in+'emissions-cmip6_'+spec_hc+'_'+type+'_surface_1750-2015_0.9x1.25_c20170322.nc'
   file_hc = path_in+'emissions-cmip6_'+spec_hc+'_'+type+'_surface_1750-2015_0.9x1.25_c20170608.nc'
   ncid = ncdf_open(file_hc)
   varname = 'emiss_'+type
   ncdf_varget,ncid,varname,emis_hc1
   ncdf_attget,ncid,/global,'molecular_weight',mw_hc
   ncdf_close,ncid
   emis_ivoc = emis_ivoc + 0.2*emis_hc1*mw_hc/mw_ivoc
 endfor
 hist = 'IVOC=0.2*('+strjoin(specs_ivoc,"+")+')'
 spec = 'IVOC'
; write new file
  ;newfile = path_out+'emissions-cmip6_'+spec+'_'+type+'_surface_1750-2015_0.9x1.25_c20170322.nc'
  newfile = path_out+'emissions-cmip6_'+spec+'_'+type+'_surface_1750-2015_0.9x1.25_c20170608.nc'
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
  ncdf_attput,ncid,/float,varid,'molecular_weight',mw_ivoc
  ncdf_attput,ncid,/char,varid,'history',hist
  ;Define global attributes
  ncdf_attput,ncid,/GLOBAL,/char,'data_title','Emissions of '+spec+' for CMIP6'
  ncdf_attput,ncid,/GLOBAL,/float,'molecular_weight',mw_ivoc
  ncdf_attput,ncid,/GLOBAL,/char,'data_creator','Jean-Francois Lamarque (lamar@ucar.edu) and Louisa Emmons (emmons@ucar.edu)'
  ncdf_attput,ncid,/GLOBAL,/char,'data_summary','Lumped HCs precursor of SOA for VBS scheme, based on fractions of various HCs.'
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
  ncdf_varput,ncid,varname,emis_ivoc
  ncdf_close,ncid


 ;SVOC = 0.6*pom_a4
   ;file_pom = path_in+'emissions-cmip6_pom_a4_'+type+'_surface_1750-2015_0.9x1.25_c20170322.nc'
   file_pom = path_in+'emissions-cmip6_pom_a4_'+type+'_surface_1750-2015_0.9x1.25_c20170608.nc'
   ncid = ncdf_open(file_pom)
   varname = 'emiss_'+type
   ncdf_varget,ncid,varname,emis_pom1
   ncdf_attget,ncid,/global,'molecular_weight',mw_pom
   ncdf_close,ncid
   emis_svoc = emis_pom1 * 0.6 *mw_pom/mw_svoc

 spec = 'SVOC'
; write new file
  ;newfile = path_out+'emissions-cmip6_'+spec+'_'+type+'_surface_1750-2015_0.9x1.25_c20170322.nc'
  newfile = path_out+'emissions-cmip6_'+spec+'_'+type+'_surface_1750-2015_0.9x1.25_c20170608.nc'
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
  ncdf_attput,ncid,/float,varid,'molecular_weight',mw_svoc
  ncdf_attput,ncid,/char,varid,'history','0.6*pom_a4'
  ;Define global attributes
  ncdf_attput,ncid,/GLOBAL,/char,'data_title','Emissions of '+spec+' for CMIP6'
  ncdf_attput,ncid,/GLOBAL,/float,'molecular_weight',mw_svoc
  ncdf_attput,ncid,/GLOBAL,/char,'data_creator','Jean-Francois Lamarque (lamar@ucar.edu) and Louisa Emmons (emmons@ucar.edu)'
  ncdf_attput,ncid,/GLOBAL,/char,'data_summary','Lumped HCs precursor of SOA for VBS scheme, based on fractions of various HCs.'
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
  ncdf_varput,ncid,varname,emis_svoc
  ncdf_close,ncid

 endfor
end
