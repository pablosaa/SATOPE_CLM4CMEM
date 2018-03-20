PROGRAM TESTBED

use clm4cmem, only: SATELLITE
implicit none

character(len=300) :: clmin, inpar
type(SATELLITE) :: SAT

clmin = "/run/media/pablo/3022-C347/clmoas.clm2.h0.2008-03-15-00900.nc"
inpar = "input"
call sat_clm4cmem(trim(clmin),trim(inpar),SAT)

print*,SAT%name
print*,SAT%TBSAT_HV(1:10,1,1,1)
print*,size(SAT%TBSAT_HV)

END PROGRAM TESTBED
