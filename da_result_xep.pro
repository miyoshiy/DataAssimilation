pro DA_result_xep

np=1000L

a=read_ascii("avg_flux_20201230.dat")
;stop

flux_da=fltarr(305,20)
flux_mds=fltarr(305,20)

k=0L
i=0L
j=0L

for i=0,19 do begin
 for j=0,304 do begin
   k=i*305+j
   flux_da[j,i]=a.field1[2,k]
   flux_mds[j,i]=a.field1[3,k]
 endfor
endfor

mape=fltarr(305)

for i=0,304 do begin
  mape(i)=0.
  for j=6,19 do begin
    mape(i)=mape(i)+abs((flux_da[i,j]-flux_mds[i,j])/flux_mds[i,j]*100)
  endfor
  mape(i)=mape(i)/14.
endfor

xL=fltarr(20)
doy=fltarr(305)

for i=0,19 do begin
 xL[i]=2.1+i*0.2
endfor

for j=0,304 do begin
  doy[j]=time_double('2018-03-01/00:00:00')+(j-1)*3600.*24.
endfor

flux_mds[where(flux_mds le -10)] = !values.f_nan


store_data,'da_result',data={x:doy,y:flux_da,v:xL}
options,'da_result',spec=1
options,'da_result',ytitle='DA_1641 keV',ysubtitle='L',ztitle='/cm2 sec str MeV'

store_data,'mds_result',data={x:doy,y:flux_mds,v:xL}
options,'mds_result',spec=1
options,'mds_result',ytitle='Arase/XEP_1641 keV',ysubtitle='L',ztitle='/cm2 sec str MeV'


store_data,'mape',data={x:doy,y:mape}
options,'mape',ytitle='MAPE [%]'


b=read_ascii("parameters_20201230.dat")
;stop
D0=fltarr(305,np)
Bwave=fltarr(305,np)

for i=0L,304L do begin
  for j=0L,np-1L do begin
    k=i*np+j
    D0(i,j)=b.field1[3,k]
    Bwave(i,j)=b.field1[4,k]
  endfor
endfor

c=read_ascii("2018Kp.dat")
xkp=fltarr(8760)
for i=0,8759 do begin
xkp[i]=time_double('2018-01-01/00:00')+(i-1)*60.*60.
endfor
store_data,'Kp',data={x:xkp,y:transpose(c.field1[3,*]/10.)}
store_data,'D0_BA',data={x:xkp,y:10^(0.506*transpose(c.field1[3,*]/10.)-9.325)}
options,'D0_BA','colors',fsc_color('red')


window,1
wset,1
!p.multi=[0,1,2]
x=findgen(305)+61.
plot,x,D0(*,0),/ylog,ytitle='D0'

for j=1,np-1 do begin
  oplot,x,D0(*,j),col=3
endfor

D0_med=fltarr(305)
D0_plus=fltarr(305)
D0_minus=fltarr(305)
for i=0,304 do begin
  D0_med(i)=mean(alog10(D0(i,*)))
  D0_plus(i)=D0_med(i)+stddev(alog10(D0(i,*)))
  D0_minus(i)=D0_med(i)-stddev(alog10(D0(i,*)))
endfor

plot,x,Bwave(*,0),/ylog,ytitle='Bw(pT)'

for j=1,np-1 do begin
  oplot,x,Bwave(*,j)
endfor

Bwave_med=fltarr(305)
Bwave_plus=fltarr(305)
Bwave_minus=fltarr(305)
for i=0,304 do begin
  Bwave_med(i)=mean(Bwave(i,*))
  Bwave_plus(i)=Bwave_med(i)+stddev(Bwave(i,*))
  Bwave_minus(i)=Bwave_med(i)-stddev(Bwave(i,*))
endfor

med_time=time_double('2018-03-01')+(x-60.)*24.*60.*60.
store_data,'D0_med',data={x:med_time,y:10^D0_med}
options,'D0_med','colors',fsc_color('blue')
store_data,'D0_plus',data={x:med_time,y:10^D0_plus}
store_data,'D0_minus',data={x:med_time,y:10^D0_minus}

store_data,'Bwave_med',data={x:med_time,y:Bwave_med}
options,'Bwave_med','colors',fsc_color('blue')
store_data,'Bwave_plus',data={x:med_time,y:Bwave_plus}
store_data,'Bwave_minus',data={x:med_time,y:Bwave_minus}


store_data,'D0_med_std',data=['D0_med','D0_plus','D0_minus']
options,'D0_med_std',ytitle='D0'

store_data,'D0_med_std_BA',data=['D0_med','D0_plus','D0_minus','D0_BA']
options,'D0_med_std_BA',ytitle='D0_BA'


store_data,'Bwave_med_std',data=['Bwave_med','Bwave_plus','Bwave_minus']
options,'Bwave_med_std',ytitle='Bwave [pT]'

options,'OMNI_HRO_1min_SYN_H',ytitle='Sym-H [nT]'

stop
end
 