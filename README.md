# cesm-cmip6-emissions
Routines used to create the CMIP6 emissions used in CESM2 (Community Earth System Model).

The forcing data sets for CMIP6 are available on the ESGF, as described here: https://pcmdi.llnl.gov/projects/input4mips/
Historical emissions are from the Community Emissions Data System (CEDS):
http://www.globalchange.umd.edu/ceds/
http://www.globalchange.umd.edu/ceds/ceds-cmip6-data/ 
CMIP6 Forcing Datasets Summary: http://goo.gl/r8up31  

## Regrid, change units, and rename original CEDS/CMIP6 files

NCL routines were used to regrid the original files to 0.9x1.25 degrees, and convert to molecules/cm2/s:
read_write_anthro.ncl; 
read_write_anthro_sectors.ncl (for SO2); 
read_write_bb.ncl; 
read_write_biofuel.ncl (for CO, to make HCN, CH3CN emissions); 
read_write_aircraft.ncl

## Anthropogenic emissions

Created files for each MOZART species of anthro emissions 1750-2015: make_emis_cmip6_anthro.pro. 
Anthro emissions are only provided through 2014, so 2014 is used for 2015.
VOC species are renamed, combined or split, to MOZART species: 
C2H6 = VOC02-ethane; 
C3H8 = VOC03-propane; 
BIGALK = VOC04-butanes + VOC05-pentanes + VOC06-hexanes-pl + VOC18-esters + VOC19-ethers; 
C2H4 = VOC07-ethene; 
C3H6 = VOC08-propene; 
C2H2 = VOC09-ethyne; 
BIGENE = VOC12-other-alke; 
BENZENE = VOC13-benzene; 
TOLUENE = VOC14-toluene; 
XYLENES = VOC15-xylene + VOC16-trimethylb + VOC17-other-arom; 
CH2O = VOC21-methanal; 
CH3CHO = VOC22-other-alka; 
CH3OH = 0.15*VOC01-alcohols; 
C2H5OH = 0.85*VOC01-alcohols; 
CH3COCH3 = 0.2*VOC23-ketones; 
MEK = 0.8*VOC23-ketones; 
HCOOH = 0.5*VOC24-acids; 
CH3COOH = 0.5*VOC24-acids; 

Special case for SO2 and sulfate:  2.5% of SO2 emitted as sulfate, industry+power plants are vertically distributed, other sectors at surface, and different sectors have different parameters for determining number: make_emis_cmip6_so2so4_anthro.pro

## Biomass burning emissions

Created files for each species of biomass burning emissions 1750-2015: make_emis_cmip6_bb.pro.
Special case for SO2 and sulfate: make_emis_cmip6_so2so4_bb.pro

## Volcanic SO2 & sulfate emissions

Continuously outgassing volcanic SO2 emissions - based on GEIA inventory; vertically distributed over various altitude depths at mountaintop elevations. Same for all months, years, so only have monthly values for a few years (1750,1850, 2000, 2100).  Total emissions are 25 Tg SO2/yr. (Jean-Francois Lamarque created final files on pressure levels).

## Soil and Ocean emissions

Soil and ocean emissions have monthly variation, but same for each year; copied from CCMI: make_emis_cmip6_other.pro
Ocean emissions for DMS, CO, NH3, C2H6, C2H4, C3H8, C3H6. 
Soil emissions for NO, NH3. 

## Biogenic emissions

CESM can simulate biogenic emissions online using MEGANv2.1 algorithms in CLM.  Climatological emissions were created for use in calculating SOAG emissions (for the simple MAM SOA scheme).  The monthly varying climatology was copied from CCMI: make_emis_cmip6_biogenic.pro

## SOAG emissions

make_emis_cmip6_soag.pro creates anthro, bb and biogenic emissions with fractions of BIGALK, BIGENE, TOLUENE, BENZENE, XYLENES, ISOP, MTERP.
From Table S2 [Liu, GMD, 2012], SOAG = 0.05*BIGALK + 0.05*BIGENE + 0.15*TOLUENE + 0.04*ISOP + 0.25*MTERP - these are mass yields, so scale emissions by mw_spec/mw_soag(=12):  SOAG = 0.05*BIGALK*(72/12) + 0.05*BIGENE*(56/12) + …

## IVOC, SVOC emissions

IVOC and SVOC are SOA precursors in the full MOZART-TS1 or -TSMLT mechanisms. 
make_emis_cmip6_svocivoc.pro creates anthro and bb emissions.

## HCN, CH3CN emissions

make_emis_cmip6_hcn.pro creates HCN and CH3CN emissions by scaling CO biofuel and biomass burning emissions.  
This uses the separate biofuel files (CO-em-SOLID-BIOFUEL-anthro_*) from CEDS. HCN is available directly from CEDS bb files.

## Aircraft emissions

make_emis_cmip6_aircraft.pro creates aircraft emissions files for NO2, bc_a4, SO2 (other species available from CEDS, but are very small).  
Emissions are zero before 1920; files include monthly values (0) for 1750, 1850, 1919. 

## MAM4 emissions

make_emis_cmip6_mam4num.pro

