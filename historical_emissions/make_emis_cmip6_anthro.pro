;
; Create anthro emissions file for 1750-2015 from CEDS/CMIP6 
; Reads files regridded to 1 deg (created with read_write_anthro.ncl)
;   with all anthro sectors summed (industry+residential+ ...)
; anthro files 1750-2014 (repeat 2014 for 2015)
; Make MOZART VOC species (combine or scale anthro VOCs)

pro make_emis_cmip6_anthro

today = bin_date(systime())
todaystr = String(today[0:2],format='(i4,"/",i2.2,"/",i2.2)')
sdate_today = String(today[0:2],format='(i4,i2.2,i2.2)')
creation_date = '20170608'
;creation_date = sdate_today
thisfile = Routine_filepath()

path_concat = '/glade/scratch/emmons/Emissions_CMIP6/concat/'
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

;MOZART mechanism species names
;specs = ['CO','NO','NH3', $  ;'SO2',
;  'BC','OC','CB1','CB2','OC1','OC2','bc_a4','pom_a4', $
specs = ['BC','OC','CB1','CB2','OC1','OC2']
;VOCs 
;specs = [ 'C2H6','C3H8','BIGALK','C2H4','C3H6','C2H2','BIGENE', $
;  'BENZENE','TOLUENE','XYLENES','CH2O','CH3CHO','CH3OH','C2H5OH', $
;  'CH3COCH3','MEK','HCOOH','CH3COOH']
;specs = ['CO','NO','NH3','bc_a4','pom_a4','SO2','SO4']

for ispec = 0,n_elements(specs)-1 do begin
  spec = specs[ispec]
  spec_ceds = [spec]
  sf = 1.0
  mw = 0.
  case spec of
   'NO': begin
      spec_ceds = ['NOx']
      mw = 30.
   end
   'C2H6': begin
     spec_ceds = ['VOC02-ethane']
     mw = 30.
    end
   'C3H8': begin
     spec_ceds = ['VOC03-propane']
     mw = 44.
   end
   'BIGALK': begin
      spec_ceds = ['VOC04-butanes','VOC05-pentanes','VOC06-hexanes-pl','VOC18-esters','VOC19-ethers']
      mw = 72.
   end
   'C2H4': begin
      spec_ceds = ['VOC07-ethene']
      mw = 28.
   end
   'C3H6': begin
      spec_ceds  = ['VOC08-propene']
      mw = 42.
   end
   'C2H2': begin
      spec_ceds  = ['VOC09-ethyne']
      mw = 26.
   end
   'BIGENE': begin
      spec_ceds  = ['VOC12-other-alke']
      mw = 56.
   end
   'BENZENE': begin
      spec_ceds = ['VOC13-benzene']
      mw = 78.
   end
   'TOLUENE': begin
      spec_ceds  = ['VOC14-toluene']
      mw = 92.
   end
   'XYLENES': begin
      spec_ceds  = ['VOC15-xylene','VOC16-trimethylb','VOC17-other-arom']
      mw = 106.
   end
   'CH2O': begin
      spec_ceds  = ['VOC21-methanal']
      mw = 30.
   end
   'CH3CHO': begin
      spec_ceds  = ['VOC22-other-alka']
      mw = 44.
   end
   'CH3OH': begin
      spec_ceds  = ['VOC01-alcohols']
      sf = 0.15
      mw = 32.
   end 
   'C2H5OH': begin
      spec_ceds  = ['VOC01-alcohols']
      sf = 0.85
      mw = 46.
   end 
   'CH3COCH3': begin
      spec_ceds  = ['VOC23-ketones']
      sf = 0.2
      mw = 58.
   end 
   'MEK': begin
      spec_ceds  = ['VOC23-ketones']
      sf = 0.8
      mw = 72.
   end 
   'HCOOH': begin
      spec_ceds  = ['VOC24-acids']
      sf = 0.5
      mw = 46.
   end 
   'CH3COOH': begin
      spec_ceds  = ['VOC24-acids']
      sf = 0.5
      mw = 60.
   end 
   'CB1': begin
    spec_ceds = ['BC']
    sf = 0.8
    mw = 12.
   end
   'CB2': begin
    spec_ceds = ['BC']
    sf = 0.2
    mw = 12.
   end
   'bc_a4': begin
    spec_ceds = ['BC']
    sf = 1.0
    mw = 12.
   end
   'OC1': begin
    spec_ceds = ['OC']
    sf = 0.5
    mw = 12.
   end
   'OC2': begin
    spec_ceds = ['OC']
    sf = 0.5
    mw = 12.
   end
   'pom_a4': begin
    spec_ceds = ['OC']
    sf = 1.4
    mw = 12.
   end
   'SO2': begin
    spec_ceds = ['SO2']
    sf = 0.975
    mw = 64.
   end
   'SO4': begin
    spec_ceds = ['SO2']
    sf = 0.025
    mw = 96.
   end
   else: begin
     spec_ceds = [spec]
     sf = 1.0
   end
  endcase 

  str1 = strjoin(spec_ceds,'+')
  if (sf lt 1) then str2 = string(sf,format='(f4.2,"*")') else str2=''
  hist_spec = 'CEDS species: '+str2+str1
help,hist_spec


  ; sum CEDS species to get MOZ species, if needed
  source_file = 'Original files: '
  for ivar = 0,n_elements(spec_ceds)-1 do begin
    file_a = path_concat + 'emissions-cmip6_'+spec_ceds[ivar]+'_anthro_surface_175001-201412_0.9x1.25_c20170608.nc'
    print,'reading ',file_a
    ncid_a = ncdf_open(file_a)
    ncdf_varget,ncid_a,'lon',lon
    ncdf_varget,ncid_a,'lat',lat
    ncdf_varget,ncid_a,'date',date_a
    ncdf_varget,ncid_a,'emiss_anthro',emis1
    if (ivar eq 0) then emiss_anthro = emis1*sf $
                   else emiss_anthro = emiss_anthro + emis1*sf
    if (mw lt 1) then ncdf_attget,ncid_a,'emiss_anthro','molecular_weight',mw
    ncdf_attget,ncid_a,/global,'source_file',source_file_b
    source_file = source_file + string(source_file_b)+', '
    ncdf_close,ncid_a
  endfor
  ntima = n_elements(date_a)
  nlon = n_elements(lon)
  nlat = n_elements(lat)
help,mw

  ;extend anthro, repeating 2014 for 2015
  emiss_anthro15 = fltarr(nlon,nlat,ntim)
  emiss_anthro15[*,*,0:ntima-1] = emiss_anthro
  emiss_anthro15[*,*,ntim-12:ntim-1] = emiss_anthro[*,*,ntima-12:ntima-1]

; write new file
  newfile = path_out+'emissions-cmip6_'+spec+'_anthro_surface_1750-2015_0.9x1.25_c'+creation_date+'.nc'
  ;newfile = path_out+'emissions-cmip6_'+spec+'_anthro-all_surface_1750-2015_0.9x1.25_c'+creation_date+'.nc'
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

  varid = ncdf_vardef(ncid,'emiss_anthro',[xid,yid,tid],/float)
  ncdf_attput,ncid,/char,'emiss_anthro','units','molecules/cm2/s'
  ncdf_attput,ncid,/char,'emiss_anthro','long_name',spec+' CEDS anthropogenic emissions'
  ncdf_attput,ncid,/float,'emiss_anthro','molecular_weight',mw
  ncdf_attput,ncid,/char,'emiss_anthro','history',hist_spec

  ;Define global attributes
  ncdf_attput,ncid,/GLOBAL,/char,'data_title','Anthropogenic emissions of '+spec+' for CMIP6'
  ncdf_attput,ncid,/GLOBAL,/float,'molecular_weight',mw
  ncdf_attput,ncid,/GLOBAL,/char,'data_creator','Jean-Francois Lamarque (lamar@ucar.edu) and Louisa Emmons (emmons@ucar.edu)'
  ncdf_attput,ncid,/GLOBAL,/char,'data_summary','Emissions from the Community Emission Data System (CEDS) have been manipulated for use in CESM2 for CMIP6. '+source_file
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
  ncdf_varput,ncid,'emiss_anthro',emiss_anthro15
  ncdf_close,ncid

endfor

end

