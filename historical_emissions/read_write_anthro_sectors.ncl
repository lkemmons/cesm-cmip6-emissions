;
; read original CEDS/CMIP6 emissions files; regrid to 1 deg; keep all anthro sector variables
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
;  dir_inp = "/glade/p/work/mmills/inputdata/emis/CMIP6/"
;  dir_inp = "/glade/p/cesm/wawg/inputdata/CMIP6/emis/"
  dir_inp = "/glade/p/cesm/wawg/inputdata/CMIP6/emis/hist/"
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
; for SO2,NOx
;    spc_name = "SO2"
;    mw = 64.
;    var_name = "SO2_em_anthro"
    spc_name = "NOx"
    mw = 30.
    var_name = "NOx_em_anthro"
;
; Sector names:
;   "0: Agriculture; 1: Energy; 2: Industrial; 3: Transportation; 4: Residential, Commercial, Other;
;      5: Solvents production and application; 6: Waste; 7: International Shipping"
;
; get anthro files
;
    files = systemfunc("ls " + dir_inp + spc_name + "-em-anthro*")
    nfils = dimsizes(files)
;
; there are 7 files per species (1750-99,1800-49,1850,1851-99,1900-49,1950-99,2000-2014)
;    do nfil=0,nfils-1
    do nfil=6,6
;
      file_inp = files(nfil)
      in = addfile(file_inp,"r")
;
      print("Processing " + file_inp)
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
      emis_sect = in->$var_name$
;
; convert from kg/m2/s to molec/cm2/s
;
;                            kg to g         m2 to cm2
;                               |                |
      emis_sect = emis_sect * 1.e3 * avog/mw * 1.e-4
;
   printVarSummary(emis_sect)
;
      sect0 = emis_sect(:,0,:,:)
      sect1 = emis_sect(:,1,:,:)
      sect2 = emis_sect(:,2,:,:)
      sect3 = emis_sect(:,3,:,:)
      sect4 = emis_sect(:,4,:,:)
      sect5 = emis_sect(:,5,:,:)
      sect6 = emis_sect(:,6,:,:)
      sect7 = emis_sect(:,7,:,:)
;
;
; regrid to CESM lat-lon grid
;
      emiss_s0 = area_conserve_remap_Wrap(lon_inp,lat_inp,sect0,lon_out,lat_out,False)
      emiss_s0@units = "molecules/cm2/s"
      emiss_s0@molecular_weight = mw
      emiss_s0@molecular_weight_units = "g/mole"
;
      emiss_s1 = area_conserve_remap_Wrap(lon_inp,lat_inp,sect1,lon_out,lat_out,False)
      emiss_s1@units = "molecules/cm2/s"
      emiss_s1@molecular_weight = mw
      emiss_s1@molecular_weight_units = "g/mole"
;
      emiss_s2 = area_conserve_remap_Wrap(lon_inp,lat_inp,sect2,lon_out,lat_out,False)
      emiss_s2@units = "molecules/cm2/s"
      emiss_s2@molecular_weight = mw
      emiss_s2@molecular_weight_units = "g/mole"
;
      emiss_s3 = area_conserve_remap_Wrap(lon_inp,lat_inp,sect3,lon_out,lat_out,False)
      emiss_s3@units = "molecules/cm2/s"
      emiss_s3@molecular_weight = mw
      emiss_s3@molecular_weight_units = "g/mole"
;
      emiss_s4 = area_conserve_remap_Wrap(lon_inp,lat_inp,sect4,lon_out,lat_out,False)
      emiss_s4@units = "molecules/cm2/s"
      emiss_s4@molecular_weight = mw
      emiss_s4@molecular_weight_units = "g/mole"
;
      emiss_s5 = area_conserve_remap_Wrap(lon_inp,lat_inp,sect5,lon_out,lat_out,False)
      emiss_s5@units = "molecules/cm2/s"
      emiss_s5@molecular_weight = mw
      emiss_s5@molecular_weight_units = "g/mole"
;
      emiss_s6 = area_conserve_remap_Wrap(lon_inp,lat_inp,sect6,lon_out,lat_out,False)
      emiss_s6@units = "molecules/cm2/s"
      emiss_s6@molecular_weight = mw
      emiss_s6@molecular_weight_units = "g/mole"
;
      emiss_s7 = area_conserve_remap_Wrap(lon_inp,lat_inp,sect7,lon_out,lat_out,False)
      emiss_s7@units = "molecules/cm2/s"
      emiss_s7@molecular_weight = mw
      emiss_s7@molecular_weight_units = "g/mole"
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
      file_out = dir_out + "emissions-cmip6_" + spc_name + "_anthro-sectors_surface_" + time_range + "_" + resolution + "_c" + tostring(creation_date) + ".nc"
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
      ncdf->date         = date
      ncdf->time         = time
      ncdf->time_bnds    = time_bnds
      ncdf->emiss_agriculture = emiss_s0
      ncdf->emiss_energy = emiss_s1
      ncdf->emiss_industry = emiss_s2
      ncdf->emiss_transport = emiss_s3
      ncdf->emiss_resident = emiss_s4
      ncdf->emiss_solvents = emiss_s5
      ncdf->emiss_waste = emiss_s6
      ncdf->emiss_shipping = emiss_s7
;
; clean up
;
      delete(emis_sect)
      delete(sect0)
      delete(sect1)
      delete(sect2)
      delete(sect3)
      delete(sect4)
      delete(sect5)
      delete(sect6)
      delete(sect7)
      delete(emiss_s0)
      delete(emiss_s1)
      delete(emiss_s2)
      delete(emiss_s3)
      delete(emiss_s4)
      delete(emiss_s5)
      delete(emiss_s6)
      delete(emiss_s7)
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
