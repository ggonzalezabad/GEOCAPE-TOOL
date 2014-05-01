program gauss_uneven

implicit real*8 (a - h, o - z)
parameter (maxpts = 50000)
dimension pos (maxpts), spec (maxpts), specmod (maxpts), &
  temp (maxpts)
character*40 input, output
specmod = 0.d0
temp = 0.d0

write (*, '(5x, a)') 'enter filespec of input file.'
read (*, '(a)') input
write (*, '(5x, a)') 'enter filespec of output file.'
read (*, '(a)') output
open (unit = 1, file = input, status = 'old')
open (unit = 2, file = output, status = 'unknown')

write (*, '(5x, a)') 'enter hw1/e.'
read (*, *) hw1e
emult = - 1.d0 / (hw1e**2)

i = 1
20 read (1, *, end = 30) pos (i), spec (i)
   i = i + 1
   go to 20
30 npoints = i - 1
write (*, *) ' npoints = ', npoints

do i = 1, npoints
  sum = 0.d0
  do j = 1, npoints
    slit = dexp (emult * (pos (j) - pos (i))**2)
    temp (j) = spec (i) * slit
    sum = sum + slit
  end do
  specmod = specmod + temp / sum
end do

do i = 1, npoints
  write (2, '(f10.5, 1p2e13.5)') pos (i), specmod (i), spec (i)
end do

close (unit = 1)
close (unit = 2)
stop
end
