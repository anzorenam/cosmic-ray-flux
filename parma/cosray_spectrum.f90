!  Calculate cosmic-ray fluxes in the atmosphere based on PARMA model
program cosray_spectrum
parameter(npart=33) ! number of applicable particle
parameter(nebin=36) ! number of energy mesh (divided by log)
parameter(nabin=90) ! number of angle mesh (divided by linear)

implicit real*8 (a-h, o-z)
dimension IangPart(0:npart)
dimension ehigh(0:nebin),emid(nebin) ! higher and middle point of energy bin
dimension ahigh(0:nabin),amid(nabin) ! higher and middle point of angular bin
dimension flux(0:nabin,0:nebin)      ! probability table (0.0 for 0, 1.0 for nabin)
data IangPart/1,2,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,4,5,5,6/ ! Angular particle ID

ip=33 ! ID particula (0:neutron, 1-28:H-Ni, 29-30:muon+-, 31:e-, 32:e+, 33:photon)
iyear=2017 ! Fecha YYYY/MM/DD
imonth=11
iday=14
glat=18.984 ! Latitude (deg), -90 =< glat =< 90
glong=-97.3 ! Longitude (deg), -180 =< glong =< 180
Alti=4.578  ! Altitude (km)

g=0.2 ! Local geometry parameter, 0=< g =< 1: water weight fraction, 10:no-earth, 100:blackhole, -10< g < 0: pilot, g < -10: cabin
s=getHP(iyear,imonth,iday,idummy) ! Solar activity (W index) of the day
r=getr(glat,glong)                ! Vertical cut-off rigidity (GV)
d=getd(alti,glat)                 ! Atmospheric depth (g/cm2), set glat = 100 for use US Standard Atmosphere 1976.

emin=1.0e1 ! Valores minimos y maximos de energia
emax=1.0e4
amin=0.0 ! Valores minimos y maximos del coseno del angulo
amax=1.0

! Arreglos para energía y ángulos
do k=1,4
  do j=0,8
    ie=j+(k-1)*9
    ehigh(ie)=(emin**k)*(j+1)
    if(ie/=0) emid(ie)=(ehigh(ie)+ehigh(ie-1))*0.5
  enddo
enddo

astep=(amax-amin)/nabin
do ia=0,nabin
  ahigh(ia)=amin+astep*ia
  if(ia/=0) amid(ia)=(ahigh(ia)+ahigh(ia-1))*0.5
enddo

flux(:,:)=0.0d0
do ie=1,nebin-1
  do ia=1,nabin
    flux(ia,ie)=getSpec(ip,s,r,d,emid(ie),g)*getSpecAngFinal(iangpart(ip),s,r,d,emid(ie),g,amid(ia))
    write(*, '(Xes10.4)',advance='no') flux(ia,ie)
  enddo
  write(*,'(A)') ' '
enddo

end program cosray_spectrum
