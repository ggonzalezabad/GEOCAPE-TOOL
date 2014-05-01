pgf90 -O3 -o cross-sections-from-hitran.x cross-sections-from-hitran.f90
pgf90 -O3 -o cross-sections-from-hitran_kelly.x cross-sections-from-hitran_kelly.f90
#
# test for IR calculations of solar absorption
#
cross-sections-from-hitran.x <<EOF
../geocape_data/HITRAN/hitran08.part
cross-section-test.out
1, 4000., 0.02, 50001., 0.500, 300., 500, 0.50, 1
EOF

#cross-sections-from-hitran.x <<EOF
#../geocape_data/HITRAN/hitran08.part
#cross-section-test.out
#1, 4000., 0.02, 50001., 0.500, 300., 500, 0.50, 1
#EOF

#cross-sections-from-hitran_kelly.x <<EOF
#../geocape_data/HITRAN/hitran08.part
#cross-section-test_kelly.out
#1, 4000., 0.02, 50001., 0.500, 300., 500, 0.50, 1
#EOF
