;This procedure reads in ARM sonde netCDF files and creates a file (TAPE5) containing the atmospheric inputs for MONORTM:
;  - pressure, temperature, altitude and RH are taken from the sonde file.
;  - CO2 is set to 380 ppmv
;  - all other molecules are set to the values of the standard atmosphere chosen by the user through the "iatm" parameters
;   1 = tropical
;   2 = midlatitude summer
;   3 = midlatitude winter
;   4 = subarctic summer
;   5 = subarctic winter
;   6 = U.S. standard

; default is U.S. standard

; Inputs:

; SONDE_FILE: ARM netcdf filename
; iatm:       index of selected atmosphere for molecules other than H2O and CO2

PRO create_monortm_input_from_sonde, SONDE_FILE,iatm=iatm

out_file = 'TAPE5'
openw, unit_2, out_file,/get_lun
 
; extract time, pressure (hPa), temp(C), rh and alt (m) from netCDF file
ncid = NCDF_OPEN(sonde_file)

varid = NCDF_VARID(ncid,'pres' )
NCDF_VARGET, ncid, varid, pres_cdf 

varid = NCDF_VARID(ncid,'tdry' )
NCDF_VARGET, ncid, varid, tdry_cdf 
tdry_cdf =  tdry_cdf + 273.16

varid = NCDF_VARID(ncid,'rh' )
NCDF_VARGET, ncid, varid, rh_cdf 

varid = NCDF_VARID(ncid,'alt' )
NCDF_VARGET, ncid, varid, alt_cdf 
alt_cdf = alt_cdf/1000.

NCDF_CLOSE, ncid

; set constants for TAPE5 output
jcharp='A'
jchart='A'
hmod = 'User defined Profile'


; set output formats
format_a = '(f10.5,g10.8,f10.4,5x,a1,a1,1x,a1,1x,28a1)'
format_3_5 = '(A10,1X,A9,A10,5x,a1,a1,1x,a1,1x,A7)'
format_b = '(8e10.3)'
format_c = '(8e15.8)'
format_d = '(i5,a24)'
jlong=' '

;;altitude we will average in groups of 5 values together to produce one unique value
;;Start averaging above 15 km
int_num = N_ELEMENTS(alt_cdf)

w_top = where(alt_cdf GE 15.000, ct_top) 
i_top = w_top(0)
i_ct = 0l

FOR p=i_top, int_num -6, 5 DO BEGIN 
    new_index = i_top + i_ct
    alt_cdf(new_index) = (alt_cdf(p)+alt_cdf(p+1)+alt_cdf(p+2)+alt_cdf(p+3)+alt_cdf(p+4))/5. 
    pres_cdf(new_index) = (pres_cdf(p)+pres_cdf(p+1)+pres_cdf(p+2)+pres_cdf(p+3)+pres_cdf(p+4))/5.
    rh_cdf(new_index) = (rh_cdf(p)+rh_cdf(p+1)+rh_cdf(p+2)+rh_cdf(p+3)+rh_cdf(p+4))/5.
    tdry_cdf(new_index) = (tdry_cdf(p)+tdry_cdf(p+1)+tdry_cdf(p+2)+tdry_cdf(p+3)+tdry_cdf(p+4))/5.
    i_ct = i_ct + 1
ENDFOR

;;remove any left over indices
alt_cdf = alt_cdf(0:new_index)
pres_cdf = pres_cdf(0:new_index)
rh_cdf = rh_cdf(0:new_index)
tdry_cdf = tdry_cdf(0:new_index)


;;final check and removal of any levels with the same pressure values
num2 = N_ELEMENTS(pres_cdf)
w_good = where((pres_cdf(0:num2-2) - pres_cdf(1:num2-1)) ge 1.e-5) ;pick out all the good indices
w_bad = where((pres_cdf(0:num2-2) - pres_cdf(1:num2-1)) lt 1.e-5)  ;pick out the bad indices
;;Print out warning of removed points
IF w_bad(0) GT -1 THEN BEGIN
  print, 'Removed level with the same pressure values: '
  print,  pres_cdf(w_bad)
  print, pres_cdf(w_bad+1)
ENDIF
;;keeping the good levels
alt_cdf = alt_cdf(w_good)
pres_cdf = pres_cdf(w_good)
rh_cdf = rh_cdf(w_good)
tdry_cdf = tdry_cdf(w_good)

; set up and fill molecule array
nmol=7
nlev = n_elements(alt_cdf)
molec = fltarr(7,nlev)
molec[0,*] = rh_cdf
molec[1,*] = replicate(380.,nlev)
jcharm='HA     '
for im=0,4 do strput,jcharm,string(iatm,format='(i1)'),im+2


;;Record 3.5 LBLRTM user defined profile
ZM_str = STRTRIM(STRING(alt_cdf,format='(f10.3)'),2)
PM_str = STRTRIM(STRING(pres_cdf,format='(f15.8)'),2)
TM_str = STRTRIM(STRING(tdry_cdf,format='(f10.3)'),2)


;;print record 3.4
printf, unit_2, nlev, hmod, format=format_d

; print record 3.5
FOR n=0, nlev-1 DO BEGIN
    printf, unit_2, zm_str(n), pm_str(n), tm_str(n), $
      jcharp, jchart, jlong, jcharm, format=format_3_5
    IF jlong EQ 'L' THEN BEGIN 
        printf, unit_2, molec(*,n), format=format_c
    ENDIF ELSE BEGIN
        printf, unit_2, molec(*,n), format=format_b
    ENDELSE
ENDFOR 

free_lun, unit_2
close,  unit_2

stop

end
