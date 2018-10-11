function galactocentrictosky,rg,theta,incli,pa,centre,asc
COMPILE_OPT IDL2
ON_ERROR, 2


; incli=inclinaison en degres
; pa=positin angle en degres
; centre = vecteur [xc,yc] centre de image
; asc = angular scale en arcsec/pixel
; theta = angle en radians
; rg=rayon galactocentric en arcsec


x0=centre[0]
y0=centre[1]
pa0=pa*!pi/180.0
inc0=incli*!pi/180.0

xsky=-rg/asc*(cos(theta)*sin(pa/180.0*!pi)+sin(theta)*cos(pa*!pi/180.0)*cos(incli*!pi/180.0))+x0
ysky=rg/asc*(cos(theta)*cos(pa/180.0*!pi)-sin(theta)*sin(pa*!pi/180.0)*cos(incli*!pi/180.0))+y0

vec=fltarr(n_elements(rg),2)
vec[*,0]=xsky
vec[*,1]=ysky
return,vec

end


;;;;;;;;;;;;Example usage

rgal=10.0 ; galactocentric radius in arcsec
  psig=!pi ; azimuthal angle in rad
inclination=60d ; in degrees
posang= 30d ; in degrees, counterclockwise from Y sky (North) axis
asc= 1.0 ; pixel size of the sky image in arcsec
xc0=128 ; coordinates of
yc0=128 ; the disk centre on the sky

  result=galactocentrictosky(rg,psig,inclination,posang,[xc0,yc0],asc)
  print,result