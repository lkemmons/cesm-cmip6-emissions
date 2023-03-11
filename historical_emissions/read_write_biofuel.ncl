;
; read original CEDS/CMIP6 emissions files; regrid to 1 deg; sum all anthro sectors to one variable
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
  ;dir_inp = "/glade/p/work/mmills/inputdata/emis/CMIP6/"
  dir_inp = "/glade/p/cesm/wawg/inputdata/CMIP6/emis/2017-05-18/"
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
      spc_name = "CO"
;      mw = 30.  ;;error discovered 2018-05-04
      mw = 28.
;
; get anthro files
;
    files = systemfunc("ls " + dir_inp + spc_name + "-em-SOLID-BIOFUEL-anthro*")
    nfils = dimsizes(files)
;
; there are 7 files per species (1750-99,1800-49,1850,1851-99,1900-49,1950-99,2000-2014)
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
        if (rank.eq.4) then 
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
      lat_inp   = in->lat
      lon_inp   = in->lon
      time      = in->time
      time_bnds = in->time_bnds
;
; read sectoral emissions
;
      wk = in->$var_name$
;
; sum over all sectors
;
      emiss_inp = dim_sum_n_Wrap(wk,1)
;
; convert from kg/m2/s to molec/cm2/s
;
;                            kg to g         m2 to cm2
;                               |                |
      emiss_inp = emiss_inp * 1.e3 * avog/mw * 1.e-4
;
; regrid to CESM lat-lon grid
;
      emiss_out = area_conserve_remap_Wrap(lon_inp,lat_inp,emiss_inp,lon_out,lat_out,False)
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
;     print(emiss_out_gl_ann(0))
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
      file_out = dir_out + "emissions-cmip6_" + spc_name + "_biofuel-anthro_surface_" + time_range + "_" + resolution + "_c" + tostring(creation_date) + ".nc"
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
      date = doubletointeger(time)
      date@units = "YYYYMMDD"
      date@long_name = "date"
      date!0 = "time"
      date&time = time
      date = tointeger(utc_date(:,0)) * 10000 + tointeger(utc_date(:,1)) * 100 + tointeger(utc_date(:,2))
;
      ncdf->emiss_anthro = emiss_out
      ncdf->date         = date
      ncdf->time         = time
      ncdf->time_bnds    = time_bnds
;
; clean up
;
      delete(wk)
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
    delete(var_name)
    delete(mw)
    delete(spc_name)
;
end
