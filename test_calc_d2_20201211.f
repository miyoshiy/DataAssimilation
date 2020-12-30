
c	2006/06/01 revised
c
	include 'dll_unit.f'
	include 'source_unit.f'
	include 'rc_unit.f'
	include 'wave_unit.f'
	include 'ae8.f'
	include 'trmfun.f'
	include "./gasdev_mel.f"
	include "./mt19937.f"
        include "./sub_pwsplot.f"
c
	parameter (nn_p=1000)
	real L,dt,dL,B0
	real E65,j65,s(100,100)
	real mec2,LL,LoE,LLx
	real a,b,c,d,f(100),e(100)
	real jj,j_max,j_min,xjj(3,3,2000),jj_sum(nn_p,100)
	real ukp(40000),udst(40000),xdst(40000),bz(40000),v(40000)
	real s_B(40000),s_M(40000)
c.......DA
	real gasdev_mel,D0(nn_p),s_p(nn_p,100,100),s_i(100,100)
	real x_ipow(nn_p)
	real Bwave0(nn_p),xhi(nn_p)
	real s_p_buf(nn_p,100,100)
        real D0_buf(nn_p),Bwave0_buf(nn_p)
        real x_ipow_buf(nn_p),xhi_buf(nn_p)
	real*8 prob(nn_p),sum_prob,diff(nn_p),diff_max
        real*8 uu,grnd
	real*8 diff_1
	integer n_number(nn_p)
c.......MDS
	real Ele(2,400,100)
	integer idum,num_sum(nn_p,100)
	character cE*6,ccmin*4,cdate*9

	do ii=1,2
	do jj=1,400
	do kk=1,100
	 Ele(ii,jj,kk)=-999.
	enddo
	enddo
	enddo
	
c.......reading DATA
	open(1,file="./XEP_1641_2018.dat",status="old")
c	read(1,*)

	do i=1,20000
	 read(1,*,end=203) iyear,imonth,idd,xLmi,E1
         idoy=montodoy(iyear-2000,imonth)+idd
	 iLrange=(anint(xLmi*10.)-10)/2+1
	 Ele(1,idoy,iLrange)=E1*1.e3
c      write(6,*) "xep ",iyear,idoy,xLmi,iLrange,Ele(1,idoy,iLrange)
	enddo
203	close(1)
c503	format(i4,1x,i2,1x,i2,f13.5,1x,)

c.......
	write(6,*) "reading end"


	open(2,file="calc_dump_20201230.dat",status="unknown")

	iflag_source=0  ! 0 -> off source, 1 -> on source
	iflag_loss=1    ! 0 -> no loss, 1 -> loss
	iflag_sample=0  ! 0 -> no sample, 1 -> write sample
	iflag_kp=0      ! 0 -> real Kp,   >=1 -> dummy_kp
	sigma=0.8
	if(iflag_source.eq.0) then
	 write(6,*) 'no inner source'
	elseif(iflag_source.eq.1) then
	 write(6,*) 'inner source'
	endif
	if(iflag_loss.eq.0) then
	 write(6,*) 'without loss'
	elseif(iflag_loss.eq.1) then
	 write(6,*) 'with loss'
	endif
	if(iflag_sample.eq.0) then
	 write(6,*) 'no sample mode'
	elseif(iflag_sample.eq.1) then
	 write(6,*) 'sampling out'
	endif
	if(iflag_kp.eq.0) then
	 write(6,*) 'real Kp'
	elseif(iflag_kp.ne.1) then
	 write(6,*) 'dummy Kp = ', iflag_kp
	endif
	write(6,*) 'sigma= (1:implicit,0.5:c-n)',sigma
c**
c	period of simulation
	xmax=365.0

	
c********************************
c	constant
c********************************
	dt=1.0/24.0  ! time step (day)
	dL=0.1       ! scale lenght
	m=55
	B0=0.352     ! Gauss
	mec2=0.511   ! MeV
	xkp=6.0
	Bwave=35.0
c	start_L=2.6
c	end_L=6.4
c	start_L=1.2
c	end_L=5.5
	start_L=2.1
	end_L=6.1
	step_L=0.2

	s_day=1.0
	s_soday=798.
	e_soday=804.

C # changed
	ndst=0
	open(1,file='./18sw.dat',status='old')
        read(1,*)
	do i=1,10000
	  read(1,501,end=200) bz(i),v(i),ukp(i),udst(i)
	  ukp(i)=ukp(i)/10.0
	  xdst(i)=i/24.0+1.0-s_day
	  if(bz(i).gt.999.) bz(i)=bz(i-1)
	  if(v(i).gt.5000.)  v(i)=v(i-1)
	  ndst=ndst+1
	enddo
200	close(1)
501	format(11x,f6.1,15x,f4.0,1x,f2.0,1x,f5.0)

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c	spetrum file open
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C # changed
	open(1,file='./18_spect.dat',status='old')

	n_spe=0
	do j=30720,39480 ! 720+24*365
	  read(1,502,end=201) iiiday,xut,s_B(j),s_M(j)
	  n_spe=n_spe+1
c          write(6,*) iiiday,xut,s_B(j),s_M(j),j
	enddo
201	close(1)
502	format(i8,1x,f8.0,1x,f12.5,f13.5)
c        write(6,*) n_spe,iiiday
        n_spe=n_spe-1
c        stop

c*****************************************:
c	energy loop
c******************************************
c	do im=1,2
c	if(im.eq.1) then
c	 E0=0.4       ! MeV
c	elseif(im.eq.2) then
c	 E0=2.5
c	endif

	E0=1.6        ! MeV


c***************************************************
c	main routine
c***************************************************
	nn=0
c	do L=1.0,6.5,0.1
c	 nn=nn+1
c	 L=3.0
c	 nn=30
c	 n_ELE=11
c***********************
c	initialize
c***********************
	 day=0.0
	 time=0.0
	 E65=0.0
	 j65=0.0
c	 xmu=E0*(E0+2.0*mec2)/(2.0*mec2*(B0/L**3))  ! MeV/G

	 do L=start_L,end_L,step_L
	  nn=int(L*10)
	  n_ELE=int((nint(L*10)-10)/2)+1
   	  xmu=E0*(E0+2.0*mec2)/(2.0*mec2*(B0/L**3))  ! MeV/G
	  do n_p=1,nn_p
	   jj_sum(n_p,nn)=0.0
	   num_sum(n_p,nn)=0
	   D0(n_p)=5.06E-10
	   prob(n_p)=0.
	   Bwave0(n_p)=35.0
	   x_ipow(n_p)=10.
	  enddo

	  do j=1,m+1
	    s(nn,j)=0.0
	  enddo

c.......AE8 as initial
	 do j=1,m+1
	   LLx=1.0+(j-1)*dL
	   xxE=mec2*(sqrt((2.0*xmu*B0)/(mec2*(LLx**3))+1.0)-1.0)
	   call ae8(LLx,1.0,xxE,xint,xdiff)
	   s(nn,j)=xdiff/(xxE*(xxE+2.0*mec2))
	 enddo
c.......

	  E65=mec2*(sqrt((2.0*xmu*B0)/(mec2*(6.5**3))+1.0)-1.0)
	  j65=(10.**8.75)*((E65*1.0e3)**(-2.56))  ! /keV sec str cm^2
	  j65=j65*1.e3                     ! /MeV sec str cm^2

	  s(nn,m+1)=j65/(E65*(E65+2.0*mec2))
	 enddo ! L loop
c**********************
c	calc
c**********************
	 do kk=30000,n_spe+30720  ! for AE8

	 time=time+dt
	 do L=start_L,end_L,step_L
           nn=int(L*10)
	   n_ELE=int((nint(L*10)-10)/2)+1
   	   xmu=E0*(E0+2.0*mec2)/(2.0*mec2*(B0/L**3))  ! MeV/G

	   e(2)=0.0
	   f(2)=0.0
c*********************
c	boundary
c*********************
	   s(nn,1)=0.0
	   E65=mec2*(sqrt((2.0*xmu*B0)/(mec2*(6.5**3))+1.0)-1.0)
	   if(kk.ge.30720) then
	    j65=(10.**s_B(kk))*((E65*1.0e3)**(s_M(kk)))
	    j65=j65*1.0e3   ! /keV -> /MeV
	   endif
	   s(nn,m+1)=j65/(E65*(E65+2.0*mec2))
c	   write(6,*) kk,s(nn,m+1),j65,s_B(kk),s_M(kk)
c           write(6,*) "goes",j65,E65,L,kk,xmu
c	   s(nn,m+1)=0.1*s(nn,m+1)

c	   s(nn,m+1)=j65/(E65*(E65+2.0*mec2))
c	   if(s(nn,m+1).eq.0.) s(nn,m+1)=1.e-2
c	   if(s(nn,m+1).gt.0.) then
c	    yy0=log10(s(nn,m+1))+gasdev_mel(1)*log10(2.0)
c           s(nn,m+1)=10.**yy0
c	   endif
	    
c******************************
c	storm phase Kp
c******************************
	   xkpmax=0.0
	   if(kk.ge.30720) then
	    xkp=ukp(int(s_day-1.0)*24+kk-30719)
	    bbz=bz(int(s_day-1.0)*24+kk-30719)
	    vv=v(int(s_day-1.0)*24+kk-30719)
	    do iij=kk-24,kk
	      xkpmax=max(xkp,xkpmax)
	    enddo
	   endif
	   if(xkp.le.1.0) xkp=1.0
	   if(iflag_kp.ne.0) then
	    xkp=real(iflag_kp)
	   endif

	   if(kk.le.32136) then
	    do j=2,m
	     q=0.0
	     LL=1.0+(j-1)*dL
	     LoE=mec2*(sqrt(2.0*xmu*B0/(LL**3)/mec2+1.0)-1.0)

	     xne=16400.*exp(-0.875*LL)              ! Abel 1999 /cm3
	     xpp=5.6-0.46*xkpmax
	     call xlifetime(LoE,LL,Bwave,xxlife)
	     tau1=xxlife
	     tau2=3.0e8*((LoE*1.0e3)**(1.5))/xne/86400.  ! Coulomb
	     tau3=
     1       LL**4/9.66*sqrt((4.0*LL-3.0)/4./LL)
     2       *(LoE+1)/sqrt(LoE*(LoE+2))
c            ! strong diffusion (Summers et al., 2004)
	     tau3=tau3*1.e4
	     tau3=tau3/60./60./24. ! sec -> day
	     tau4=5.0/xkp
	     

	     if(xkp.ge.6.0) xkp=6.0
	     dll_1=dll_E(xkp,LL+dL/2.0,xmu)+dll_B(xkp,LL+dl/2.0)
	     dll_2=dll_E(xkp,LL-dL/2.0,xmu)+dll_B(xkp,LL-dl/2.0)


	     a=-dll_1*(LL*LL)/(dL*dL)*dt*sigma
	      if(LL.le.xpp) then
	        if(tau1.lt.1000.0.and.kk.gt.10000) then
	         b=(dll_1+dll_2)*(LL*LL)/(dL*dL)*dt*sigma+1.0
     1            +dt/tau1+dt/tau2            !loss
	        else
	         b=(dll_1+dll_2)*(LL*LL)/(dL*dL)*dt*sigma+1.0
     1            +dt/tau2            !loss
	        endif
	      else
	       b=(dll_1+dll_2)*(LL*LL)/(dL*dL)*dt*sigma+1.0
     1          +dt/tau3+dt/tau4
	      endif
	     c=-dll_2*(LL*LL)/(dL*dL)*dt*sigma
	     d=s(nn,j)+So*dt+
     1       (1.-sigma)*(s(nn,j+1)-2.*s(nn,j)+s(nn,j-1))*
     2       dt*(dll_1+dll_2)/2.*(LL*LL)/(dL*dL)
	     e(j+1)=-a/(c*e(j)+b)
	     f(j+1)=(d-c*f(j))/(c*e(j)+b)
	   enddo

	   do j=m,2,-1
	     s(nn,j)=e(j+1)*s(nn,j+1)+f(j+1)
	   enddo
	  endif

	  if(kk.eq.32136) then
	   do j=m,2,-1
	     s_i(nn,j)=s(nn,j)
	   enddo
	  endif

	  if(kk.gt.32136) then ! 32136= doy=89
	   do n_p=1,nn_p
	    idum=n_p+kk-32000
	    if(kk.eq.32137) then
	       D0(n_p)=(0.506*xkp-9.325)
	       D0(n_p)=10**D0(n_p)
	       Bwave0(n_p)=log10(35.0)
	       Bwave0(n_p)=10.**Bwave0(n_p)	
	       x_ipow(n_p)=10.
	       xhi_tmp=xhi(n_p)
	    endif

	    if(kk.gt.32136) then
	     do j=m,2,-1
	      s(nn,j)=s_p(n_p,nn,j)
	     enddo
	    elseif(kk.eq.32136) then
	     do j=m,2,-1
	      s(nn,j)=s_i(nn,j)
	     enddo
	    endif

	    do j=2,m
	     q=0.0
	     LL=1.0+(j-1)*dL
	     LoE=mec2*(sqrt(2.0*xmu*B0/(LL**3)/mec2+1.0)-1.0)

	     xne=16400.*exp(-0.875*LL)              ! Abel 1999 /cm3
	     xpp=5.6-0.46*xkpmax
	     call xlifetime(LoE,LL,Bwave0(n_p),xxlife)
	     tau1=xxlife
	     tau2=3.0e8*((LoE*1.0e3)**(1.5))/xne/86400.  ! Coulomb
	     tau3=
     1       LL**4/9.66*sqrt((4.0*LL-3.0)/4./LL)
     2       *(LoE+1)/sqrt(LoE*(LoE+2))
c            ! strong diffusion (Summers et al., 2004)
	     tau3=tau3*1.e4
	     tau3=tau3/60./60./24. ! sec -> day
	     tau4=xhi(n_p)/xkp

c	     dll_1=D0(n_p)*(LL+dl/2.0)**10
c	     dll_2=D0(n_p)*(LL-dl/2.0)**10
	     dll_1=D0(n_p)*(LL+dl/2.0)**x_ipow(n_p)
	     dll_2=D0(n_p)*(LL-dl/2.0)**x_ipow(n_p)
c             write(6,*) "dLL1_",dll_1,dll_2,tau1,tau2,tau3

	     a=-dll_1*(LL*LL)/(dL*dL)*dt*sigma
	      if(LL.le.xpp) then
	        if(tau1.lt.1000.0.and.kk.gt.10000) then
	         b=(dll_1+dll_2)*(LL*LL)/(dL*dL)*dt*sigma+1.0
     1            +dt/tau1+dt/tau2            !loss
	        else
	         b=(dll_1+dll_2)*(LL*LL)/(dL*dL)*dt*sigma+1.0
     1            +dt/tau2            !loss
	        endif
	      else
	       b=(dll_1+dll_2)*(LL*LL)/(dL*dL)*dt*sigma+1.0
     1          +dt/tau3+dt/tau4
	      endif
	     c=-dll_2*(LL*LL)/(dL*dL)*dt*sigma
	     d=s(nn,j)+So*dt+
     1       (1.-sigma)*(s(nn,j+1)-2.*s(nn,j)+s(nn,j-1))*
     2       dt*(dll_1+dll_2)/2.*(LL*LL)/(dL*dL)
	     e(j+1)=-a/(c*e(j)+b)
	     f(j+1)=(d-c*f(j))/(c*e(j)+b)
	   enddo

	   do j=m,2,-1
	     s(nn,j)=e(j+1)*s(nn,j+1)+f(j+1)
             if(s(nn,j).lt.0.0) then
                 s(nn,j)=0.0
             endif
	     s_p(n_p,nn,j)=s(nn,j)
	   enddo
	  enddo  ! n_p loop
	  endif

	  enddo  ! L loop


	  if(kk.gt.32136) then
	   xx=(time-720./24.)+1.0
	     if(int(xx).eq.int(xx_old)) then
	       do n_p=1,nn_p
	        do L=start_L,end_L,step_L
	          nn=int(L*10)
	          n_ELE=int((nint(L*10)-10)/2)+1
	          jj=s_p(n_p,nn,int((L-1.0)/0.1+1))*E0*(E0+2.0*mec2)
c                  write(6,*) xx,jj,L,n_p,E0
	          if(jj.gt.0.) then
	           jj=log10(jj)
	           jj_sum(n_p,nn)=jj_sum(n_p,nn)+jj
	           num_sum(n_p,nn)=num_sum(n_p,nn)+1
	          endif
	        enddo
	       enddo
	     endif
	     if(int(xx).ne.int(xx_old)) then
	      do n_p=1,nn_p
	       diff(n_p)=0.
	       do L=start_L,end_L,step_L
	        nn=int(L*10)
	        n_ELE=int((nint(L*10)-10)/2)+1
	        if(num_sum(n_p,nn).ne.0) then
	         if(int(L*10).eq.nn) then
	          if(Ele(1,int(xx_old),n_ELE).gt.0.) then
	           if(L.le.5.0) rr=0.5
	           if(L.gt.5.0) rr=3.0
	           diff_1=
     1             log10(Ele(1,int(xx_old),n_ELE))
     2             -jj_sum(n_p,nn)/num_sum(n_p,nn)
	           diff(n_p)=diff(n_p)-(diff_1**2/rr/rr)
	          else
	           goto 826
	          endif
	         endif
	        else
c	         write(6,*) L,int(xx_old),-999.,nn
	        endif
	        jj=s_p(n_p,nn,int((L-1.0)/0.1+1))*E0*(E0+2.0*mec2)
	        if(jj.gt.0.) then
	         jj=log10(jj)
	        endif
	        jj_sum(n_p,nn)=jj
	        num_sum(n_p,nn)=1
	       enddo
	      enddo

	      diff_max=-9999.
	      do n_p=1,nn_p
	        diff_max=max(diff_max,diff(n_p))
	      enddo
	      do n_p=1,nn_p
	        prob(n_p)=exp(diff(n_p)-diff_max)
	      enddo

c...........  skip resampling
c	      goto 826
	      sum_prob=0.
	      do n_p=1,nn_p
	       sum_prob=sum_prob+prob(n_p)
	      enddo
	      if(sum_prob.eq.0.) then
	      write(6,*) "00000"
              goto 824
	      endif
	      prob(1)=prob(1)/sum_prob
	      do n_p=2,nn_p
	       prob(n_p)=prob(n_p-1)+prob(n_p)/sum_prob
	      enddo
	      do n_p=1,nn_p
	       D0_buf(n_p)=D0(n_p)
	       Bwave0_buf(n_p)=Bwave0(n_p)
	       x_ipow_buf(n_p)=x_ipow(n_p)
	       xhi_buf(n_p)=xhi(n_p)
	       do L=start_L,end_L,step_L
	        nn=int(L*10)
	        do j=2,m
	         s_p_buf(n_p,nn,j)=s_p(n_p,nn,j)
	        enddo
	       enddo 
	      enddo

c..........initialize
	      do n_p=1,nn_p
	        n_number(n_p)=n_p
	      enddo

	      do n_p=1,nn_p
	       idum=n_p+kk-32000
c825	       uu=grnd()
825	       uu=gasdev_mel(idum)
	       if(uu.lt.0.or.uu.gt.1.0) goto 825
	       do jjj=1,nn_p
	        if(prob(jjj).gt.uu) then
	         D0_buf(n_p)=D0(jjj)	
	         Bwave0_buf(n_p)=Bwave0(jjj)
	         x_ipow_buf(n_p)=x_ipow(jjj)
	         xhi_buf(n_p)=xhi(jjj)
	         n_number(n_p)=jjj
	         do L=start_L,end_L,step_L
	          nn=int(L*10)
	          do j=2,m
	           s_p_buf(n_p,nn,j)=s_p(jjj,nn,j)
	          enddo
	         enddo
	         goto 823	
	        endif
	       enddo
823	       continue
	      enddo
	      do n_p=1,nn_p
	       D0(n_p)=D0_buf(n_p)
	       Bwave0(n_p)=Bwave0_buf(n_p)
	       x_ipow(n_p)=x_ipow_buf(n_p)
	       xhi(n_p)=xhi_buf(n_p)
	       do L=start_L,end_L,step_L
	        nn=int(L*10)
	        do j=2,m
	         s_p(n_p,nn,j)=s_p_buf(n_p,nn,j)
	        enddo
	       enddo
	      enddo
826	      do n_p=1,nn_p
	       do L=start_L,end_L,step_L
	        nn=int(L*10)
	        n_ELE=int((nint(L*10)-10)/2)+1
	        jj=s_p(n_p,nn,int((L-1.0)/0.1+1))*E0*(E0+2.0*mec2)
	        if(Ele(1,int(xx_old),n_ELE).gt.0.0) then
	        write(2,*) 
c     1         int(xx_old),jj_sum(n_p)/num_sum(n_p),
     1          int(xx_old),jj,
     2          log10(Ele(1,int(xx_old),n_ELE)),D0(n_p),
     3          Bwave0(n_p),x_ipow(n_p),prob(n_p),n_number(n_p),
     4          n_p,L,n_ELE
	       else
	        write(2,*) 
     1           int(xx_old),jj,
     2           Ele(1,int(xx_old),n_ELE),D0(n_p),
     3           Bwave0(n_p),x_ipow(n_p),prob(n_p),n_number(n_p),
     4           n_p,L,n_ELE
	       endif
	      enddo ! L
	      idum=n_p+kk-32000
c             D0(n_p)=log10(5.06E-10)
c	      D0(n_p)=log10(D0(n_p))+gasdev_mel(idum)
c	      D0(n_p)=(0.506*xkp-9.325)+gasdev_mel(idum)*3.0
611           urandom=gasdev_mel(idum)
              if(urandom.gt.1.0.or.urandom.lt.-1.0) then 
                goto 611
              endif
	      D0(n_p)=(0.506*xkp-9.325)+3.0*urandom
	      D0(n_p)=10**D0(n_p)
c	      Bwave0(n_p)=log10(Bwave0(n_p))+gasdev_mel(idum)
c	      Bwave0(n_p)=log10(35.0)
	      Bwave0(n_p)=log10(35.0)+2.*urandom
	      Bwave0(n_p)=10.**Bwave0(n_p)
c              write(6,*) 'D0&Bwave ',
c     1 n_p,D0(n_p),Bwave0(n_p),urandom
c	      x_ipow(n_p)=10.+2.*gasdev_mel(idum)
	      x_ipow(n_p)=10.
	      xhi_tmp=xhi(n_p)+gasdev_mel(idum)
827	      continue
	      if(xhi_tmp.le.0.) then
	        xhi_tmp=xhi(n_p)+gasdev_mel(idum)
	        goto 827
	      else
	         xhi(n_p)=xhi_tmp
	      endif
	     enddo ! n_p
	    endif
824	continue
	   xx_old=xx
	 endif
	enddo    ! kk loop


	close(2)
c       enddo     ! L loop

c       enddo     ! m loop


	stop
	end

	     
        subroutine year(cdate,iyear)
c
        character cdate*9
        iyear=ichara(cdate(9:9),1)
        return
        end

        subroutine doy(cdate,idoy)
c
        character cdate*9

        if(cdate(4:6).eq."Jan") idoy=0
        if(cdate(4:6).eq."Feb") idoy=31
        if(cdate(4:6).eq."Mar") idoy=59
        if(cdate(4:6).eq."Apr") idoy=90
        if(cdate(4:6).eq."May") idoy=120
        if(cdate(4:6).eq."Jun") idoy=151
        if(cdate(4:6).eq."Jul") idoy=181
        if(cdate(4:6).eq."Aug") idoy=212
        if(cdate(4:6).eq."Sep") idoy=243
        if(cdate(4:6).eq."Oct") idoy=273
        if(cdate(4:6).eq."Nov") idoy=303
        if(cdate(4:6).eq."Dec") idoy=334
        idoy_m=ichara(cdate(1:2),2)
        if(idoy_m.le.0) idoy_m=ichara(cdate(2:2),1)
        idoy=idoy+idoy_m
        return
        end

	     

	  
	  


	
