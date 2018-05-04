pro pconv,datfile,relcfile,timfile,hc,notemp=notemp
;
;plots the convergence log (from unit 9) for stellar atmospheres
; (output from TLUSTY)
;
;datfile, relcfile - names for units 7 and 9 outputs, respectively:
;if RELCFILE is missing, and only DATFILE is present, then the program
; assumes that the names are: DATFILE.7 and DATFILE.9 (i.e. the parameter
; DATFILE contains the core of the name, without any extensions);
;if both are present, they contains the names in full (useful for owerwriting 
; the standard notation convention);
;if both are absent, the program assumes the names being: fort.7 and fort.9;
;
;use hc for hardcopy (0=default, 6=portrait,5=landscape)
;use ymin and ymax to override default y-axis limits

if n_params(0) lt 4 then hc=0
if n_params(0) eq 1 then begin
   datf=datfile+'.7'
   relcfile=datfile+'.9'
   timfile=datfile+'.69'
   datfile=datf
endif
if n_params(0) eq 0 then begin
   datfile='fort.7'
   relcfile='fort.9'
   timfile='fort.69'
endif


;GET MASS DEPTH ARRAY, NITER, ND

get_lun,lun
a=''
openr,lun,datfile

data=intarr(2)
readf,lun,data
nd=data(0)
mass=fltarr(nd)
readf,lun,mass
free_lun,lun

get_lun,lun
openr,lun,relcfile

CT=FLTARR(ND,200) & CP=CT & CNE=CT & CMAX=CT   ;DEFINE OUTPUT ARRAYS
tmax=fltarr(200) & mmax=tmax & itrn=tmax

for i=1,3 do readf,lun,a   ;skip 3 lines
data=fltarr(9)
k=-1
while not eof(lun) do begin
k=k+1
 for j=1,nd do begin
  readf,lun,data
  ct(j-1,k)=data(2)
  cp(j-1,k)=data(4)
  cne(j-1,k)=data(3)
  cmax(j-1,k)=data(6)
 end
endwhile

niter=k+1   ;actual number of iterations

if niter eq 1 then ct=fltarr(nd,1)+ct(*,0) else  ct=ct(*,0:k)
if niter eq 1 then cp=fltarr(nd,1)+cp(*,0) else  cp=cp(*,0:k)
if niter eq 1 then cne=fltarr(nd,1)+cne(*,0) else  cne=cne(*,0:k)
if niter eq 1 then cmax=fltarr(nd,1)+cmax(*,0) else  cmax=cmax(*,0:k)


;REVERSE TO AGREE WITH MASS DEPTH COORDINATE

for k=0,niter-1 do begin
ct(0,k)=reverse(ct(*,k))
cp(0,k)=reverse(cp(*,k))
cne(0,k)=reverse(cne(*,k))
cmax(0,k)=reverse(cmax(*,k))
end
if n_elements(notemp) eq 0 then lct=alog10(abs(ct))
lcmax=alog10(abs(cmax))
for k=0,niter-1 do begin
  if n_elements(notemp) eq 0 then tmax(k)=max(lct(*,k))
  mmax(k)=max(lcmax(*,k))
  itrn(k)=k+1.
end

free_lun,lun

 ITER=STRTRIM(NITER,2)

x=alog10(mass)
!xtitle='Log Depth (Mass)

if hc ne 0 then set_plot,'ps'
if hc eq 6 then device,/landscape

;!p.multi=[0,3,2]

;1st plot
if n_elements(notemp) eq 0 then begin
set_viewport,.15,.38,.575,.875
!mtitle='Temperature
 !ytitle='Relative Change
 plot,x,ct(*,0),yrange=[min(ct),max(ct)]
 if niter gt 1 then for k=1,niter-1 do oplot,x,ct(*,k),line=k
!noeras=1
endif

;2nd plot
if n_elements(notemp) eq 0 then begin
set_viewport,.44,.67,.575,.875
 !ytitle=''
  yt='Log !9!!!3 Relative Change !9!!!3'
 plot,x,lct(*,0),yrange=[min(lct),max(lct)]
 if niter gt 1 then for k=1,niter-1 do oplot,x,lct(*,k),line=k
endif

;3rd plot
 set_viewport,.15,.38,.15,.45
!mtitle='Maximum in State Vector
 !ytitle='Relative Change
 plot,x,cmax(*,0),yrange=[min(cmax),max(cmax)]
 if niter gt 1 then for k=1,niter-1 do oplot,x,cmax(*,k),line=k
!noeras=1

;4th plot
set_viewport,.44,.67,.15,.45
!mtitle='Maximum in State Vector
 !ytitle=''
 yt='Log !9!!!3 Relative Change !9!!!3'
 plot,x,lcmax(*,0),yrange=[min(lcmax),max(lcmax)]
 if niter gt 1 then for k=1,niter-1 do oplot,x,lcmax(*,k),line=k

;5th plot
 if n_elements(notemp) eq 0 then begin
set_viewport,.73,.96,.575,.875
!mtitle='Temperature
!xtitle='iteration'
 !ytitle=''
 yt='Log !9!!!3 Relative Change !9!!!3'
 plot,itrn(0:niter-1),tmax(0:niter-1),psym=-1
endif

;6th plot
set_viewport,.73,.96,.15,.45
!mtitle='Maximum in State Vector
!xtitle='iteration'
 !ytitle=''
 yt='Log !9!!!3 Relative Change !9!!!3'
 plot,itrn(0:niter-1),mmax(0:niter-1),psym=-1

;
tit=relcfile
spawn,'rm -f tmp'
spawn,'date >tmp'
close,2
openr,2,'tmp'
a=''
readf,2,a
;if hc eq 0 then xyouts,0.8,0.98,a,/normal else $
;xyouts,0.8,0.98,a,/normal,size=0.5
close,2
spawn,'rm -f tmp'
;xyouts,0.5,0.98,tit,alignment=0.5,size=2,/normal
;
;!p.multi=0
!p.title=''
!noeras=0
!xtitle=''
set_viewport
;
openr,2,timfile
while not eof(2) do readf,2,t1,t2,t3
tit=tit+'  time='+strcompress(string(format='(i)',(t3)))+' sec'
xyouts,0.5,0.92,tit,alignment=0.5,size=2,/normal
close,2
;

 if hc ne 0 then begin
  device,/close
  set_plot,'x'
  spawn,'mac_ps idl.ps'
; spawn,'rm -f idl.ps'
 end



return
end


