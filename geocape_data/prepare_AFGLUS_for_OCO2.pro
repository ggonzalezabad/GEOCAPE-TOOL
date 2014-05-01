; This is to prepare an atmospheric input profile for OCO-2
; sensitivity studies

; first read the AFGL US standard profiles
close, 1  & openr, 1, '/data/scibo/xliu/ATMOS/AFGL/afglus.dat'
;# AFGL atmospheric constituent profile. U.S. standard atmosphere 1976. ( AFGL-TR-86-0110)       
;#     z(km)      p(mb)        T(K)    air(cm-3)    o3(cm-3)  o2(cm-3)    h2o(cm-3)    co2(cm-3)     no2(cm-3)
afglus = dblarr(9, 50)
lines=strarr(2)
readf, 1, lines
readf, 1, afglus

da = where(afglus(0, *) lt 65.0)
afglus = afglus(*, da)

for i = 0, 8 do begin
    afglus(i, *) = reverse(reform(afglus(i, *)))
endfor

; this profile does not have methane, read another TES profile and use
; TES profile format, but with fewer layers
;TES_File_ID = "L2: True Full State Vector Data"
;Data_Size = 86 x 15
;surfaceTemperature = 275.99756
;PSUR = 1014.41
;boresightNadirAngle = 0.
;latitude = -48.90103
;longitude = 186.48939
;End_of_Header  ****  End_of_Header  ****  End_of_Header  ****  End_of_Header
;Level Pressure     Altitude    TATM     H2O       CO2    O3        N2O        CO        CH4       O2     NO         NO2        HNO3       OCS
;-     hPa          m           K        VMR       VMR    VMR       VMR        VMR       VMR       VMR    VMR        VMR        VMR        VMR
;  1   1014.4100000 0.00000000  2.7942E2 6.2727E-3 3.8E-4 3.0990E-8 3.1423E-7  8.8651E-8 1.6357E-6 1.E-20 4.9541E-12 1.1741E-11 2.2375E-11 1.E-20


close, 1  & openr, 1, '/data/scibo/xliu/GEOCAPE/RETGRP/VISnew/geocape_data/TES_Data/AtmProfiles_Data_Run551C_Seq0023_Scan002.asc'
lines=strarr(10)
readf, 1, lines
tes=dblarr(15, 86)
readf, 1, tes
close, 1


; produce the new input data
nl = n_elements(afglus(0, *))
afglustes = dblarr(15, nl)
afglustes(0, *) = findgen(nl)+1
afglustes(1, *) = afglus(1, *)
afglustes(2, *) = afglus(0, *)
afglustes(3, *) = afglus(2, *)
hasidxs = [0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0]
yidxs   = [0, 0, 0, 0, 6, 7, 4, 0, 0, 0, 5, 0, 8, 0, 0]

for i = 4, 14 do begin
    if hasidxs(i) gt 0 then begin
        j = yidxs(i)
        afglustes(i, *) = afglus(j, *) / afglus(3, *)
    endif else begin
        afglustes(i, *) = spline(tes(2, *), tes(i, *), afglustes(2, *))
    endelse
endfor

; use AFGL CH4 profile which is not in the used AFGL US standard
; profile
afglustes(9, 0:8) = 1.70
afglustes(9, 9:37) = [1.69, 1.69, 1.68, 1.66, 1.65, 1.63, 1.61, 1.58, 1.55, 1.52, 1.46, $
                      1.42, 1.36, 1.27, 1.19, 1.12, 1.06, 0.987, 0.914, 0.830, 0.746,   $
                      0.662, 0.564, 0.461, 0.363, 0.277, 0.210, 0.165, 0.150]
afglustes(5, *) = afglustes(5, *) * 380. / 330. ; scale to 380 ppm
afglustes(9, *) = afglustes(9, *) * 1.0E-6

close, 2 & openw, 2, 'AtmProfiles_Data_AFGLUS.asc'
printf, 2, lines, format='(A)'
for i = 0, nl - 1 do begin
    printf, 2, fix(afglustes(0, i)), afglustes(1:14, i), format='(I3, F14.7,F12.8,12E12.4)'
endfor
close, 2

stop

end


