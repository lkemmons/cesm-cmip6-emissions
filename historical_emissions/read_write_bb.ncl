;
; read original CEDS/CMIP6 biomass burning emissions and regrid to 1deg
;
begin
;
; constants
;
  avog = 6.022e23
  earth_rad  = 6.371e8 ; in cm
  earth_area = 4. * acos(-1.) * earth_rad * earth_rad
;
; input directory
;
  dir_inp = "/glade/p/cesm/wawg/inputdata/CMIP6/emis/"
;
; output directory
;
  dir_out = "/glade/scratch/emmons/Emissions_CMIP6/regrid_mol/"
;
; template for lat-lon CESM grid
;
  resolution = "0.9x1.25"
  in = addfile("/glade/p/cesmdata/cseg/inputdata/atm/cam/inic/fv/template_" + resolution + "_L32.nc","r")
  lat_out = in->lat
  lon_out = in->lon
  delete(in)
;
; create weights from lat (to be used for global average)
;
  gw = NormCosWgtGlobe(lat_out)
;
; make sure output dir exists
;
  system("mkdir -p " + dir_out)
;
; loop over species
;
;  do nspc=0,30
;  do nspc=4,30
  do nspc=31,31
;
    if ( nspc.eq.31 ) then
      spc_name = "CO2"
      mw = 44.
    end if
    if ( nspc.eq.0 ) then
      spc_name = "BC"
      mw = 12.
    end if
    if ( nspc.eq.1 ) then
      spc_name = "CO"
      mw = 28.
    end if
    if ( nspc.eq.2 ) then
      spc_name = "NH3"
      mw = 17.
    end if
    if ( nspc.eq.3 ) then
      spc_name = "NOx"
;      mw = 46.  ;6/11/2018 - bb NOx emitted as NO
      mw = 30.
    end if
    if ( nspc.eq.4 ) then
      spc_name = "OC"
      mw = 12.
    end if
    if ( nspc.eq.5 ) then
      spc_name = "SO2"
      mw = 64.
    end if
    if ( nspc.eq.6 ) then
      spc_name = "NMVOC_C10H16"
      mw = 136.
    end if
    if ( nspc.eq.7 ) then
      spc_name = "NMVOC_C2H2"
      mw = 26.
    end if
    if ( nspc.eq.8 ) then
      spc_name = "NMVOC_C2H4"
      mw = 28.
    end if
    if ( nspc.eq.9 ) then
      spc_name = "NMVOC_C2H4O"
      mw = 44.
    end if
    if ( nspc.eq.10 ) then
      spc_name = "NMVOC_C2H5OH"
      mw = 46.
    end if
    if ( nspc.eq.11 ) then
      spc_name = "NMVOC_C2H6"
      mw = 30.
    end if
    if ( nspc.eq.12 ) then
      spc_name = "NMVOC_C2H6S"
      mw = 62.
    end if
    if ( nspc.eq.13 ) then
      spc_name = "NMVOC_C3H6"
      mw = 42.
    end if
    if ( nspc.eq.14 ) then
      spc_name = "NMVOC_C3H6O"  ;acetone, CH3COCH3
      mw = 58.
    end if
    if ( nspc.eq.15 ) then
      spc_name = "NMVOC_C3H8"
      mw = 44.
    end if
    if ( nspc.eq.16 ) then
      spc_name = "NMVOC_C5H8"  ;isoprene
      mw = 68.
    end if
    if ( nspc.eq.17 ) then
      spc_name = "NMVOC_C6H6"  ;Benzene
      mw = 78.
    end if
    if ( nspc.eq.18 ) then
      spc_name = "NMVOC_C7H8"  ;Toluene
      mw = 92.
    end if
    if ( nspc.eq.19 ) then
      spc_name = "NMVOC_C8H10"  ;Xylenes
      mw = 106.
    end if
    if ( nspc.eq.20 ) then
      spc_name = "NMVOC_CH2O"
      mw = 30.
    end if
    if ( nspc.eq.21 ) then
      spc_name = "NMVOC_CH3COCHO"
      mw = 72.
    end if
    if ( nspc.eq.22 ) then
      spc_name = "NMVOC_CH3COOH"
      mw = 60.
    end if
    if ( nspc.eq.23 ) then
      spc_name = "NMVOC_CH3OH"
      mw = 32.
    end if
    if ( nspc.eq.24 ) then
      spc_name = "NMVOC_HCN"
      mw = 27.
    end if
    if ( nspc.eq.25 ) then
      spc_name = "NMVOC_HCOOH"
      mw = 46.
    end if
    if ( nspc.eq.26 ) then
      spc_name = "NMVOC_Higher_Alkanes"
;      mw = 72.  ;6/11/2018 - discovered orig files use mw=12
      mw = 12.
    end if
    if ( nspc.eq.27 ) then
      spc_name = "NMVOC_Higher_Alkenes"
      mw = 56.  ;6/11/2018 - discovered orig files use mw=12
      mw = 12.
    end if
    if ( nspc.eq.28 ) then
      spc_name = "NMVOC_HOCH2CHO"  ;hydroxyacetaldehyde
      mw = 60.
    end if
    if ( nspc.eq.29 ) then
      spc_name = "NMVOC_MEK"
      mw = 72.
    end if
    if ( nspc.eq.30 ) then
      spc_name = "NMVOC_Toluene_lump"  ;AROMATICS
      mw = 92.
    end if
;
; get bb files (2 files per species)
;
    files = systemfunc("ls " + dir_inp + spc_name + "-em-biomassburning*")
    nfils = dimsizes(files)
;
    do nfil=0,nfils-1
;
      file_inp = files(nfil)
      in = addfile(file_inp,"r")
;
      print("Processing " + file_inp)
;
; find emissions variable name
;
      vnames = getfilevarnames(in)
      nnames = dimsizes(vnames)
      do n=0,nnames-1
        var = in->$vnames(n)$
        dimv = dimsizes(var)
        rank = dimsizes(dimv)
        if (rank.eq.3) then 
           var_name = vnames(n)
        end if
        delete(var)
        delete(dimv)
        delete(rank)
      end do
      print(var_name)
;
; read coordinates
;
      lat_inp   = in->latitude
      lon_inp   = in->longitude
      time      = in->time
      time_bnds = in->time_bnds
;
; read bb emissions
;
      emiss_inp = in->$var_name$
;
; reverse lat
;
      lat_inp = lat_inp(::-1)
      emiss_inp = emiss_inp(:,::-1,:)
;
; convert from kg/m2/s to molec/cm2/s
;
;                            kg to g         m2 to cm2
;                               |                |
      emiss_inp = emiss_inp * 1.e3 * avog/mw * 1.e-4
      emiss_inp@_FillValue = 1.e20
;      printVarSummary(emiss_inp)
;
      emiss_inp = where(ismissing(emiss_inp),0.,emiss_inp)
;
; regrid to CESM lat-lon grid
;
;      print(dim_max(dim_max(dim_max(emiss_inp))))
      emiss_out = area_conserve_remap_Wrap(lon_inp,lat_inp,emiss_inp,lon_out,lat_out,False)
;      print(dim_max(dim_max(dim_max(emiss_out))))
      delete(emiss_out@units)
      emiss_out@units = "molecules/cm2/s"
      emiss_out@molecular_weight = mw
      emiss_out@molecular_weight_units = "g/mole"
;
; compute global average (in molec/cm2/s)
; annual average in Tg/year
;
      emiss_out_gl = wgt_areaave(emiss_out,gw,1.0,1)
      dims = dimsizes(emiss_out_gl)
      nyr = dims(0)/12
      emiss_out_gl_ann = emiss_out_gl(0:nyr-1)
      do n=0,nyr-1
        n1 = n * 12
        n2 = n1 + 11
        emiss_out_gl_ann(n) = dim_avg_Wrap(emiss_out_gl(n1:n2))
      end do
      emiss_out_gl_ann = emiss_out_gl_ann * earth_area * 86400. * 365 * mw * 1.e-12/avog
;
      print(emiss_out_gl_ann(0))
;
; extract time range string from input file name
;
      length = strlen(file_inp)
      char   = stringtochar(file_inp)
      time_range = char(length-16:length-4)
      delete(char)
;
; define creation date
;
      creation_date = systemfunc ("date +%Y%m%d")
;
; output to netcdf file
;
      file_out = dir_out + "emissions-cmip6_" + spc_name + "_biomassburning_surface_" + time_range + "_mol_" + resolution + "_c" + tostring(creation_date) + ".nc"
      system("/bin/rm -f " + file_out) ; remove any pre-existing file
      ncdf = addfile(file_out,"c")
;
; find file name (without path)
;
      strs = str_split(files(nfil),"/")
      indx = dimsizes(strs)
      file_name = strs(indx-1)
      delete(strs)
;
; define global attributes
;
      file_att               = True
      file_att@source_file   = file_name
      file_att@creation_date = creation_date
      file_att@author        = "Jean-Francois Lamarque (lamar@ucar.edu) and Louisa Emmons (emmons@ucar.edu)"
      fileattdef(ncdf,file_att)
      filedimdef(ncdf,"time",-1,True)
;
; define date array
;
      utc_date = cd_calendar(time, 0)
      date = tointeger(time)
      date@units = "YYYYMMDD"
      date@long_name = "date"
      date!0 = "time"
      date&time = time
      date = tointeger(utc_date(:,0)) * 10000 + tointeger(utc_date(:,1)) * 100 + tointeger(utc_date(:,2))
;
;      print("before write")
      ncdf->emiss_bb = emiss_out
;      print("after write")
      ncdf->date         = date
      ncdf->time         = time
      ncdf->time_bnds    = time_bnds
;
; clean up
;
      delete(emiss_inp)
      delete(emiss_out)
      delete(emiss_out_gl)
      delete(emiss_out_gl_ann)
      delete(in)
      delete(file_inp)
      delete(time)
      delete(time_bnds)
      delete(utc_date)
      delete(date)
      delete(ncdf)
;
; end loop over files
;
    end do
;
; clean up
;
    delete(mw)
    delete(spc_name)
;
; end loop over species
;
  end do   
;
end
