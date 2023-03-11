; make vertically distributed emissions
;  - evenly distributed to altitude as function of PFT

pro make_emis_cmip6_vertbbztop

today = bin_date(systime())
todaystr = String(today[0:2],format='(i4,"/",i2.2,"/",i2.2)')
sdate_today = String(today[0:2],format='(i4,i2.2,i2.2)')
creation_date = '20170322'
thisfile = Routine_filepath()

file_pfts = '/glade/p/work/emmons/fire_inputs/fgeos5.cesm1_5_beta02_fire_pftfrac.nc'
ncid = ncdf_open(file_pfts)
ncdf_varget,ncid,'lon',lon_p
ncdf_varget,ncid,'lat',lat_p
ncdf_varget,ncid,'PFT_tot',pft_tot2
pft_tot = reform(pft_tot2[*,*,0])
nlon_p = n_elements(lon_p)
nlat_p = n_elements(lat_p)
pftfracs = fltarr(nlon_p,nlat_p,18)
for ipft = 1,17 do begin
 varname = 'PFT_'+string(ipft,format='(i2.2)')
 ncdf_varget,ncid,varname,pft1
 pftfracs[*,*,ipft] = reform(pft1[*,*,0])
endfor
ncdf_close,ncid
help,pftfracs
for ilon=0,nlon_p-1 do begin
 for ilat=0,nlat_p-1 do begin
   if (abs(Total(pftfracs[ilon,ilat,*])- pft_tot[ilon,ilat]) gt 0.001) then begin
     print,lon_p[ilon],lat_p[ilat],Total(pftfracs[ilon,ilat,*]),pft_tot[ilon,ilat]
   endif
 endfor
endfor

file_injht = '/glade/p/work/emmons/fire_inputs/fire_injheight_CLM_PFTs_c20151130.nc'
ncid = ncdf_open(file_injht)
ncdf_varget,ncid,'FireVertDistribution',fire_injht
ncdf_varget,ncid,'PFT_Num', pft_num
ncdf_varget,ncid,'PFT_Name', pftname
ncdf_varget,ncid,'Region', region
ncdf_varget,ncid,'Reg_Name', regname
ncdf_varget,ncid,'Month', month
ncdf_varget,ncid,'Altitude', alt_injht
ncdf_close,ncid
;print,alt_injht
nmon = n_elements(month)
nalt = n_elements(alt_injht)-1
nreg = n_elements(region)  ;7
npft = n_elements(pft_num) ;79


;ztop = fltarr(18)
;ztop = [0, 1.5, 1.5, 1.5, 0.7, 1.0, 0.7, 1.0, 1.5, 1.0, 1.0, 1.0, 0.7, 0.7, 0.7, 0.7, 0.7]
;not_vegetated    500 m                       
;needleleaf_evergreen_temperate_tree     4000 m
;needleleaf_evergreen_boreal_tree    4000 m
;needleleaf_deciduous_boreal_tree    3000 m    
;broadleaf_evergreen_tropical_tree     2500 m  
;broadleaf_evergreen_temperate_tree   3000 m   
;broadleaf_deciduous_tropical_tree     2500 m  
;broadleaf_deciduous_temperate_tree  3000 m    
;broadleaf_deciduous_boreal_tree      3000 m   
;broadleaf_evergreen_shrub   2000 m            
;broadleaf_deciduous_temperate_shrub  2000 m  
;broadleaf_deciduous_boreal_shrub    2000 m    
;c3_arctic_grass   1000 m                      
;c3_non-arctic_grass  1000 m              
;c4_grass   1000 m                             
;c3_crop      1000 m
;(and all new crops: 1000m)
ztop = [0.5, 4.0, 4.0, 3.0, 2.5, 3.0, 2.5, 3.0, 3.0, 2.0, 2.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0]
for ipft = 0,16 do begin
  print,ipft,': ',ztop[ipft],' : ',String(pftname[*,ipft])
endfor


dalt = 0.5E5 ;cm (0.5 km)
nalt = 10
altitude_int = findgen(nalt+1)*(dalt/1.e5) ;km
altitude_mid = findgen(nalt)*(dalt/1.e5) + 0.5*(dalt/1.e5)
print,altitude_mid
;ztop_alt_ind = Fix(ztop /(dalt/1.e5))
;print,ztop_alt_ind
;print,ztop

for ipft = 1,16 do begin
  prof1 = fltarr(nalt)
  ind = where(altitude_mid lt ztop[ipft],nz)
  prof1[ind] = 1./nz
  print,ipft,':',(prof1)
endfor
stop


;specs = ['bc_a4','pom_a4','so4_a1','SO2','SOAG','SOAGx1.5', $
;specs = [  'NH3','DMS','CO','NO',
;specs = ['C2H2','C2H4','C2H6','C3H6','C3H8','BIGALK','BIGENE', $
;  'BENZENE','TOLUENE','XYLENES','CH3OH','C2H5OH','CH2O','CH3CHO','CH3COCH3','MEK', $
;  'HCOOH','CH3COOH','HCN','CH3CN','CH3COCHO','GLYALD',
specs = ['ISOP','MTERP','IVOC','SVOC']

;path_orig = '/glade/p/cesm/chwg_dev/emmons/CMIP6_emissions_1750_2015/'
path_orig = '/glade/p/cesm/chwg_dev/emmons/CMIP6_emissions_1750_2015_v20170322/'
path_new = path_orig

for ispec = 0,n_elements(specs)-1 do begin
 spec = specs[ispec]
 print,spec
 file_emis = File_search(path_orig+'emissions-cmip6_'+spec+'_bb_surface_1750-2015_0.9x1.25_c20170322.nc',count=nf)
 print,file_emis
 source_file = file_emis
 if (nf ne 1) then stop
 ncid = ncdf_open(file_emis[0])
 ncdf_varget,ncid,'lon',lon
 ncdf_varget,ncid,'lat',lat
 ncdf_varget,ncid,'date',date
 ncdf_varget,ncid,'time',time
 ncdf_varget,ncid,'emiss_bb',emis_bb_surf
 ncdf_attget,ncid,/global,'molecular_weight',mw
 ncdf_close,ncid
 nlon = n_elements(lon)
 nlat = n_elements(lat)
 ntim = n_elements(date)

 emis_bb_vert = fltarr(nlon,nlat,nalt,ntim)

 for ilon=0,nlon-1 do begin
  for ilat=0,nlat-1 do begin

   for itim=0,ntim-1 do begin
    emis_surf = emis_bb_surf[ilon,ilat,itim]
    if (emis_surf gt 0) then begin
     emis_vert = fltarr(nalt)
     imon= Fix( Strmid(String(date[itim],format='(i8)'),4,2) )-1
     ;print,total(pftfracs[ilon,ilat,*]),reform(pftfracs[ilon,ilat,*])
     for ipft=1,16 do begin
      ;pftfr1 = pftfracs[ilon,ilat,ipft]/pft_tot[ilon,ilat]
      pftfr1 = pftfracs[ilon,ilat,ipft]/Total(pftfracs[ilon,ilat,*])
      ;print,ipft,pftfr1
      if (pftfr1 gt 0.0) then begin
       injprof1 = fltarr(nalt)
       ind = where(altitude_mid lt ztop[ipft],nz)
       injprof1[ind] = 1./nz
       emis_vert[*] = emis_vert[*] + emis_surf * pftfr1 * injprof1 /dalt
      endif 
     endfor 
     emis_bb_vert[ilon,ilat,*,itim] = emis_vert
     ;stop
    endif 
   endfor

  endfor
 endfor

 file_vert = path_new+'emissions-cmip6_'+spec+'_bb_vertical_1750-2015_0.9x1.25_c'+creation_date+'.nc'
 print,'creating ',file_vert

ncid = ncdf_create(file_vert,/clobber)
    lonid = ncdf_dimdef(ncid, 'lon', nlon)
    latid = ncdf_dimdef(ncid, 'lat', nlat)
    timeid = ncdf_dimdef(ncid, 'time', /unlimited) 
    levidint = ncdf_dimdef(ncid, 'altitude_int', nalt+1)
    levidmid = ncdf_dimdef(ncid, 'altitude', nalt)
 ; Define dimension variables with attributes
    xvarid = ncdf_vardef(ncid,'lon',[lonid],/float)
    ncdf_attput, ncid, xvarid, 'units',/char, 'degrees_east'
    ncdf_attput, ncid, xvarid, 'long_name',/char, 'Longitude'
    yvarid = ncdf_vardef(ncid,'lat',[latid],/float)
    ncdf_attput, ncid, yvarid, 'units',/char, 'degrees_north'
    ncdf_attput, ncid, yvarid, 'long_name',/char, 'Latitude'
    tvarid = ncdf_vardef(ncid,'time',[timeid],/float)
    ncdf_attput, ncid, tvarid, 'units',/char, 'days since 1750-01-01 00:00:00'
    ncdf_attput, ncid, tvarid, 'long_name',/char, 'Time'
    ncdf_attput, ncid, tvarid,/char, 'calendar', 'Gregorian'
    tvarid = ncdf_vardef(ncid,'date',[timeid],/long)
    ncdf_attput, ncid, tvarid, 'units',/char, 'YYYYMMDD'
    ncdf_attput, ncid, tvarid, 'long_name',/char, 'Date'
    xvarid = ncdf_vardef(ncid,'altitude_int',[levidint],/float)
    ncdf_attput, ncid, xvarid, 'units',/char, 'km'
    ncdf_attput, ncid, xvarid, 'long_name',/char, 'Altitude Interfaces'
    xvarid = ncdf_vardef(ncid,'altitude',[levidmid],/float)
    ncdf_attput, ncid, xvarid, 'units',/char, 'km'
    ncdf_attput, ncid, xvarid, 'long_name',/char, 'Altitude Midlevel'
    varid = ncdf_vardef(ncid, 'emiss_bb', [lonid,latid,levidmid,timeid], /FLOAT)
    ncdf_attput,ncid,varid, 'long_name', /char, 'CEDS '+spec+' biomass burning emissions distributed vertically'
    ncdf_attput,ncid,varid, 'units', /char, 'molecules/cm3/s'
    ncdf_attput,ncid,varid,/float,'molecular_weight',mw

  ;Define global attributes
  ncdf_attput,ncid,/GLOBAL,/char,'data_title','Biomass burning vertically distributed emissions of '+spec+' for CMIP6'
  ncdf_attput,ncid,/GLOBAL,/float,'molecular_weight',mw
  ncdf_attput,ncid,/GLOBAL,/char,'data_creator','Jean-Francois Lamarque (lamar@ucar.edu) and Louisa Emmons (emmons@ucar.edu)'
  ncdf_attput,ncid,/GLOBAL,/char,'data_summary','Emissions from the Community Emission Data System (CEDS) have been manipulated for use in CESM2 for CMIP6. Vertical distribution uniform to a max height dependent on PFT fractions (using file: '+file_injht+'. Original emissions file: '+source_file
  ncdf_attput,ncid,/GLOBAL,/char,'cesm_contact','Louisa Emmons or Simone Tilmes'
  ncdf_attput,ncid,/GLOBAL,/char,'creation_date',creation_date
  ncdf_attput,ncid,/GLOBAL,/char,'update_date',sdate_today
  ncdf_attput,ncid,/GLOBAL,/char,'history','original files have been regridded, units changed, concatenated, some species combined or renamed.'
  ncdf_attput,ncid,/GLOBAL,/char,'data_script',thisfile
  ncdf_attput,ncid,/GLOBAL,/char,'data_source_url','http://www.globalchange.umd.edu/ceds/ceds-cmip6-data/'
  ncdf_attput,ncid,/GLOBAL,/char,'data_reference','van Marle et al., GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-32/)'

    ncdf_control,ncid, /endef
    ncdf_varput,ncid,'lon',lon
    ncdf_varput,ncid,'lat',lat
    ncdf_varput,ncid,'time',time
    ncdf_varput,ncid,'date',date
    ncdf_varput,ncid,'altitude_int', altitude_int
    ncdf_varput,ncid,'altitude', altitude_mid
    ncdf_varput,ncid,'emiss_bb',emis_bb_vert
ncdf_close,ncid

endfor

end
