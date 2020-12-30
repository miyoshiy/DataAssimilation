

	parameter (nn_p=1000)
	integer doy(400),n_number(400,nn_p)
	real flux_mds(400,nn_p,100),flux_calc(400,nn_p,100)
        real bwave(400,nn_p,100)
        real D0(400,nn_p,100)
	real*8 prob(400,nn_p,100)
	real s(400,nn_p,100),s1(400,nn_p),s2(400,nn_p)
        real s_buf(400,nn_p,100),s1_buf(400,nn_p)
        real s2_buf(400,nn_p)
	real L,xL(400,100)
	real ave(400,100)
	real xx(10000),kp(10000),dst(10000),DLLm(10000)
	integer num_ave(400,100)

c	open(1,file="calc_dump_20190519.dat",status="old")
c	open(1,file="calc_dump_20201228.dat",status="old")
	open(1,file="calc_dump_20201230.dat",status="old")


c	start_L=2.6
c	end_L=4.1
c	start_L=2.6
c	end_L=6.0
	start_L=2.1
        end_L=6.1
	step_L=0.2


c	nn_L=0
c	do L=start_L,end_L,step_L
c	nn_L=nn_L+1
c	do i=1,nn_p
c	 read(1,*)
c	enddo
c	enddo

c	do i=1,9000
c	 read(1,*)
c	enddo

c	nn_L=18
	nn_L=20
        nskip=6000*10

        do i=1,nskip
	  read(1,*,end=201,err=822) 
        enddo

	n_dat=0
	do i=1,1000
	 do j=1,nn_p
	  do iL=1,nn_L
	  read(1,*,end=201,err=822) 
     1    doy(i),flux_calc(i,j,iL),flux_mds(i,j,iL),
     2    D0(i,j,iL),bwave(i,j,iL),x_in,
     3    prob(i,j,iL),n_number(i,j),n_p,xL(i,iL),n_ELE
	  enddo
	 enddo
	 n_dat=n_dat+1
	enddo
201	close(1)
822	continue

c	do iL=1,nn_L
c	write(6,*) iL,xL(1,iL)
c	enddo

c	goto 826
	iLength=10
	do iL=1,nn_L
	do i=1,n_dat-iLength
	  do j=1,nn_p
	    s(i,j,iL)=flux_calc(i,j,iL)
	    s1(i,j)=D0(i,j,iL)
	    s2(i,j)=bwave(i,j,iL)
	  enddo
	  do k=1,iLength
	   do j=1,nn_p
	    nn=n_number(i+k,j)
	    s_buf(i,j,iL)=s(i,nn,iL)
	    s1_buf(i,j)=s1(i,nn)
	    s2_buf(i,j)=s2(i,nn)
	   enddo
	   do j=1,nn_p
	    s(i,j,iL)=s_buf(i,j,iL)
	    s1(i,j)=s1_buf(i,j)
	    s2(i,j)=s2_buf(i,j)
	   enddo
	  enddo
	  do j=1,nn_p
	   flux_calc(i,j,iL)=s(i,j,iL)
	   D0(i,j,iL)=s1(i,j)
	   bwave(i,j,iL)=s2(i,j)
	  enddo
	enddo
	enddo
826	continue

        open(2,file="parameters_20201230.dat",status="unknown")
          do i=1,n_dat
           do j=1,nn_p
            write(2,*) iL,i,j,D0(i,j,1),bwave(i,j,1)
           enddo
          enddo
        close(2)



	do iL=1,nn_L
	 do i=1,n_dat
	  ave(i,iL)=0.
	  num_ave(i,iL)=0
	  do j=1,nn_p
	   if(flux_calc(i,j,iL).gt.0.0) then
	    ave(i,iL)=ave(i,iL)+log10(flux_calc(i,j,iL))
	    num_ave(i,iL)=num_ave(i,iL)+1
	   endif
	  enddo
	 enddo
	enddo

	do iL=1,nn_L
	 do i=1,n_dat
	  if(num_ave(i,iL).ne.0.0) then
	   ave(i,iL)=ave(i,iL)/num_ave(i,iL)
	  else
	   ave(i,iL)=-99.9
	  endif
	 enddo
	enddo

        open(2,file="avg_flux_20201230.dat",status="unknown")
        do iL=1,nn_L
         do i=1,n_dat
          write(2,*) iL,i,ave(i,iL),flux_mds(i,1,iL)
         enddo
        enddo
        close(2)



	 
        stop
        end
