
function rayon1, siz, inc, pa, centre, asc,rmax
COMPILE_OPT IDL2
ON_ERROR, 2

; siz=a vector giving the size of the 2D sky image [nx,ny]
; inc=inclination in degres
; pa=position angle of the major axis in degrees
; centre = a vector [xc,yc] giving the sky position of the centre of the
disk
; asc = angular scale in arcsec/pixel of a sky-plane pixel
; rmax = maximum galactocentric radius to reach (in arcsec)


x0=centre[0]
y0=centre[1]
pa0=pa*!pi/180.0
inc0=inc*!pi/180.0

rg=fltarr(siz[0],siz[1]) & xg=rg & yg=rg & psig=rg
  j=replicate(1,siz[0]) # indgen(siz[1])
  i=indgen(siz[0]) # replicate(1,siz[1])

xg=-(i-x0)*asc*sin(pa0)+(j-y0)*cos(pa0)*asc
yg=-((i-x0)*asc*cos(pa0)+(j-y0)*asc*sin(pa0))/cos(inc0)


rg=sqrt((xg^2+yg^2))
tanpsi=yg/xg
psig=atan(tanpsi)
psig[where(xg lt 0.0)]+=!pi
psig[where(xg gt 0.0 and yg lt 0.0)]+=2.0*!pi
uu=where(xg eq 0.0 and yg gt 0.0)
if uu[0] ne -1 then psig[uu]=!pi/2.0
uu2=where(xg eq 0.0 and yg lt 0.0)
if uu2[0] ne -1 then psig[uu2]=3*!pi/2.0


uu=where(rg gt rmax)
if uu[0] ne -1 then begin
rg[uu]=!values.f_nan
yg[uu]=!values.f_nan
xg[uu]=!values.f_nan
psig[uu]=!values.f_nan
endif

vec=fltarr(siz[0],siz[1],4)
vec[*,*,0]=rg ; in arcsec
vec[*,*,1]=xg ; in arcsec
vec[*,*,2]=yg ; in arcsec
vec[*,*,3]=psig ; in radians
  return, vec

end


;;;;;;;; EXAMPLE USAGE ;;;;;;;;;;

rmax=100.0 ; arcsec
asc= 1.0 ; pixel size of the sky image in arcsec
naxis1=256
naxis2=256
xc0=128 ; coordinates of
yc0=128 ; the disk centre on the sky
inclination=60d ; in degrees
posang= 30d ; in degrees, counterclockwise from Y sky (North) axis

toto=rayon1([naxis1,naxis2], inclination, posang, [xc0,yc0], asc,rmax)
radiusmap=reform(toto[*,*,0])
tetamap=reform(toto[*,*,3]) ; in radians
xgmap=reform(toto[*,*,1])
ygmap=reform(toto[*,*,2])


; create fits file output for e.g. the radius
writefits,'radiusgal.fits',radiusmap
