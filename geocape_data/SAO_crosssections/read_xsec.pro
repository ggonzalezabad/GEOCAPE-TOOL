; write cross sections for OClO, IO, GLYX, H2O, O2, Solar

datadir = '/home/xliu/OzoneFit/tbl/'
infiles = datadir + $
           ['Unused/OMSAO_OClO_213K.dat', $
            'Unused/OMSAO_IO_298K_UBremen.dat', $
            'Unused/OMSAO_CHOCHO_296K_Volkamer.dat', $
            'Unused/OMSAO_H2O_0.9atm-280k_390-800nm_0.04nmFWHM_0.01nm_i0corr.dat', $
            'Unused/OMSAO_O2_0.5atm-253k_600-800nm_0.04nmFWHM_0.01nm_i0corr.dat', $
            'OMSAO_NKPNO_SolarSpec_JQSRT2011.dat']
            
            
nheaders = [21, 79, 53, 1, 1, 8]
outfiles = ['oclo_213k_sciafm.dat', $
            'io_298k_bremen.dat', $
            'glyx_296k_volkamer.dat', $
            'h2o_0.9atm_280k_hitran.dat', $
            'o2_0.5atm_253k_hitran.dat', $
            'chance_solarspec_jqsrt2011.dat']

nfile = n_elements(infiles)
nl = 0
swav=0.0
ewav=0.0
scale=0.0

for i = 0, nfile - 1 do begin
    close, 1 & openr, 1, infiles(i)
    nh = nheaders(i)
    lines=strarr(nh)
    readf, 1, lines
    readf, 1, nl, swav, ewav, scale
    data=dblarr(2, nl)
    readf, 1, data
    close, 1
    
    close, 2 & openw, 2, outfiles(i)
    for j = 0, nl - 1 do begin
        printf, 2, data(0, j), data(1, j)*scale, format='(F12.6, E16.6)' 
    endfor
    close, 2

endfor

stop
end
