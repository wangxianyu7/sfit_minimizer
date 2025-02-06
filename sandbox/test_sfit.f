	implicit real*8 (a-h,o-z)
        parameter (ndiv=500)
 	parameter (nobs=2) ! Change
	parameter (nparms=9)
	parameter (ntop=3*nobs + nparms)
	parameter (ntab=10000)
	dimension b0tab(ntab),b1tab(ntab),db0tab(ntab),db1tab(ntab)
        dimension epslog(ndiv),gvec(ndiv),dgvec(ndiv)
	character*80 line
	character*1 linechar(80)
	character*8 word
	logical found
	character*20 linemoa
	character*1 char
	character*12 char12
	character*4 char4
	real*8 a(ntop),b(ntop,ntop),c(ntop,ntop),d(ntop)
     *    ,dum(ntop,ntop),f(ntop),s(ntop),da(ntop)
	dimension date(40000,nobs),flux(40000,nobs),err(40000,nobs)
	dimension see(40000,nobs)
	dimension chi2(nobs),ktot(nobs),knt(nobs)
        dimension chi2max(nobs),kmax(nobs)
	dimension kfix(ntop),kgamma(nobs)
	dimension qlong(nobs),qlat(nobs)
	logical bad(40000,nobs),badd(40000,nobs)
	character*21 filename
	character*3  rootname
	character*1  bandname
	character*200 direc
	integer dirlen
	logical getlat
	logical dave
c 12714.52
c        data a/5763.326250,  0.002631,239.351848,  0.000100,  0.440000,
c     *    0.720000,  0.260000,  0.102964, -0.058164,  0.210228,
c     *    4.122066,  0.0000,0.210228,
c     *    4.122066,  0.0000,0.210228,
c     *    4.122066,  0.0000/
c 12714.7897      6.8510      7.4321
c	data a/1697.5, 0.3, 25.,  0.000100,  0.440000,
c	data a/8645.1, 0.005, 25.,  0.00600,  0.440000,
	data a/1697.2, 0.006, 25.,  0.006500,  0.440000,
     *    0.720000,  0.260000,  0., 0.,  1.3,
     *    0.0,  0.000000,  1.0, 0.0, 0.0/
c       Results
c       data a/1697.00, 0.00499546365, 25.2149772,
c     +    0.00669937561, 0.44, 0.72, 0.26, 0., 0.,
c     +    1.19886310, 0.101130762, 0.0, 2.34866444,
c     +    0.650272159, 0.0/
	
c	data a/1697., 0.25, 25.1,  0.000100,  0.440000,
c     *    0.720000,  0.260000,  0.35, 0.62,  1.2,
c     *    0.1,  0.000000,  2.35, 0.67, 0.0/
! Change
	
	data pi/3.14159265358979323846d0/
	direc = '/Users/jyee/PycharmProjects/sfit_minimizer/data/MMTest/'
	dirlen=Index(direc, 'Test/')
	dirlen=dirlen+4
	
	getlat = .false.
        if(getlat)call getlatlon(qlat,qlong,nobs)
	alpha =  (18 + 00./60 + 00./3600)*15
	delta = -(30 + 00./60 + 00./3600)
c	t0par = 8645.
	t0par = a(1)
	kntfac = 0
	nobmoa = -1
	offmoa = 0
c       set nobbad=-1 for regular use, and =nob to find bad data
	nobbad = -1
	open(1,file='/Users/jyee/Tables/Microlensing/b0b1.dat',
     *      status='old')
	do 176 i=1,ntab
	   read(1,*)z,b0tab(i),b1tab(i),db0tab(i),db1tab(i)
 176	continue
	close(1)
	do 177 i=1,ntop
	   kfix(i) = 0
 177	continue
	dave = .false.
c	write(6,*)'enter u0'
c	read(5,*)a(2)
c	write(6,*)'enter rho'
c	read(5,*)a(4)
	! Change
c	kfix(2) = 1 ! u0
c	kfix(3) = 1 ! tE
c	a(3) = 277
c	kfix(4) = 1 ! rho
	kfix(5) = 1
	kfix(6) = 1
	kfix(7) = 1
	kfix(8) = 1 ! piEN
	kfix(9) = 1 ! piEE
c	kfix(10) = 1 ! fix source flux for first dataset
c	kfix(11) = 1 ! fix blending for first dataset
c	kfix(13) = 1 ! fix source flux for second dataset
	do 178 i=1,nobs
	   kfix(9 + 3*i)=1
 178	continue
	   if(nobmoa.gt.0)kfix(9 + 3*nobmoa)=0
	i = 2
c	kfix(9 + 3*i-1)=1
c	a(9 + 3*i-1)=0
	i = 3
c	kfix(9 + 3*i-1)=1
c	a(9 + 3*i-1)=0
        epstop = 50
        epsbot = 1./200.
        call mkgvec(epslog,gvec,dgvec,ndiv,epstop,epsbot)
	do 5 nob = 1,nobs
	   close(1)
	   if(nob.ne.nobmoa)then
	      if(nob.eq.1)filename='FSPL_Obs_1_I.pho'
	      if(nob.eq.2)filename='FSPL_Obs_2_V.pho'
	   else
	      filename = 'moa025.phot'
	   endif
	   open(1,file=direc(1:dirlen)//filename,status='old')
	   nform = 1                      ! OGLE-uFUN
	   if(nob.eq.nobmoa)nform = 2     ! MOA
	   call getrootband(qlat1,qlong1,bandname,rootname,filename)
	   if(nob.eq.1)bandname = 'I'
	   if(nob.eq.2)bandname = 'V'
	   if(.not.getlat)then
	      qlat(nob) = qlat1
	      qlong(nob)= qlong1
	   endif
c	   read(filename,181)bandname,rootname
c 181	   format(8x,a1,1x,a3)
	   if(rootname.eq.'raw')nform = 3 ! PLANET
	   if(rootname.eq.'pys')nform = 4 ! PLANET
	   kgamma(nob)                    = 1 ! I band
	   if(bandname.eq.'V')kgamma(nob) = 2 ! V band
	   if(bandname.eq.'H')kgamma(nob) = 3 ! H band
	   if(bandname.eq.'R')kgamma(nob) = 4 ! R band
	   if(bandname.eq.'U')kgamma(nob) = 4 ! unfiltered = R 
	   i = 0
	   kk = 0
 3	   continue
 	   read(1,1,end=4)line
	   kk = kk+1
c	   write(6,1)line
 1	   format(a80)
	   read(line,2)char
 2	   format(80a1)
	   if(char.eq.'#')go to 3
	   read(line,2)linechar
	   length=5
	   word=' nan '
	   read(line,*)dat
	   datp = dat
	   if(datp.gt.12000)datp = datp - 2450000d0
c	   if(datp.lt.8500)go to 3
c	   if(datp.gt.8700)go to 3
c	   if(dat.gt.10000.and.dat.lt.2458000)go to 3
	   call vetword2(found,word,length,linechar,1)
	   if(found)go to 3
	   if(nob.eq.nobmoa)then
	      read(line,*)dat,qi,e,char12,char4
c	      if(dat.lt.2455000)go to 3
	      if(char4.eq.'nan '.or.char4.eq.' nan') go to 3
	      read(char4,*)se
	      see(i+1,nob)= se - 4
	   else
	      if(nform.eq.1)read(line,*)dat,qi,e
	      if(nform.eq.3)read(line,*)qi,e,dat,se
	      if(nform.eq.4)read(line,*)dat,qi,e,se
	      see(i+1,nob)=1.5
	   endif
	   if(.not.dave)then
	      if(nob.eq.nobmoa)e = e*1.
	      if(nob.eq.1)e = e*1.
	      if(nob.eq.2)e = e*1.
	      if(nob.eq.3)e = e*1.
	      if(nob.eq.4)e = e*1.
	      if(nob.eq.5)e = e*1
	      if(nob.eq.6)e = e*1
	      if(nob.eq.7)e = e*1
	   endif
	   i=i+1
	   bad(i,nob) = .false.
	   badd(i,nob) = .false.
	   if(dat.gt.2440000)dat = dat - 2450000.d0
	   if(nob.eq.nobmoa)then
	      flux(i,nob) = qi/6500. + offmoa
	      err(i,nob) = e/6500.
	      call hjdcor(hjd,dat,alpha,delta)
	      dat = hjd
	   else
	      if(rootname.eq.'pys')then
		 flux(i,nob) = 1-qi/10000.
		 err(i,nob) = e/10000.
	      else
c		 e = sqrt(e**2+ 0.006**2)
		 flux(i,nob) = 10.**(0.4*(18.-qi))
		 err(i,nob) = e*flux(i,nob)*log(10.)/2.5
	      endif
	   endif
	   date(i,nob)= dat
c	   write(6,*)nob,i,date(i,nob),flux(i,nob)
	   if(qi.eq.0.)bad(i,nob) = .true.
	   go to 3
 4	   continue
	   knt(nob) =i
 5	continue
	write(6,*)knt
c	read(5,*)xyz
c	write(6,*) a
c	write(6,*) date(10, 1), flux(10, 1), err(10, 1)
	if(.not.dave)then
c	   bad(1037,1) = .true. ! 3492.0041
c	   bad(1007,nobmoa) = .true. ! 3900.07
	endif
	if(nobbad.ne.-1)kdatbad = knt(nobbad)
	jloop = 0
 8	continue
	if(nobbad.ne.-1)then
	   kntfac = 100
	   fac = 0.05
	   jloop = jloop + 1
	   kntgood = kdatbad - jloop + 1
	   sigmax = sqrt(2*log(kntgood/3./sqrt(3.14159/2)))
	   sigmax = sqrt(2*log(kntgood/sigmax/sqrt(3.14159/2)))
	endif
	iloop = 0

	write(61, *)'# nob, k, t, f(4), amp*(b0-gamma*b1), db0, db1'
	write(62, *)'# nob, k, t, f(1),f(2),f(3),f(4), damp,dfdx,y,sig2'
	write(63, *)'# nob, k, t, x, x2, amp, b0, b1'
 6	continue
	iloop = iloop + 1
	do 10 i=1,ntop
	   d(i) = 0
	   do 10 j=1,ntop
	      b(i,j) = 0
 10	continue
	do 20 i=1,nobs
	   chi2(i) = 0
	   chi2max(i) = 0
	   ktot(i) = 0
 20	continue
	t0 = a(1)
	beta = a(2)
	te = a(3)
        rho= a(4)
        gammai=a(5)
	gammav=a(6)
	gammah=a(7)
	piex=a(8)
	piey=a(9)
	do 85 nob=1,nobs
	   fs = a(nparms-2+3*nob)
	   fb = a(nparms-1+3*nob)
	   fsee = a(nparms+3*nob)
	   do 25 i=nparms+1,ntop
	      f(i)=0
 25	   continue
	   do 84 k =1,knt(nob)
	      if(bad(k,nob))go to 84
	      if(badd(k,nob))go to 84
	      t = date(k,nob)
	      fl = flux(k,nob)
	      sig2 = err(k,nob)**2
 26	      format(3f12.6)
	      tau = (t-t0)/te
	      call geta(qn,qe,t,alpha,delta,t0par)
	      call gett(qnp,qep,t,qlat(nob),qlong(nob),alpha,delta)
	      qn = qn + qnp
	      qe = qe + qep
	      dtau =  piex*qn + piey*qe
	      dbeta= -piex*qe + piey*qn
	      dbeta = -dbeta
	      taup = tau + dtau
	      betap = beta + dbeta
	      x2 = betap**2 + taup**2
	      x = sqrt(x2)
              z = x/rho
	      call
     *        getb0p(b0,b1,db0,db1,z,b0tab,b1tab,db0tab,db1tab,ntab)      
	      amp = (x2+2)/x/sqrt(x2+4)
	      damp = -8*amp/x/(x2+2)/(x2+4)
c              f56 =  fs*amp*db1*x/rho**2
              f56 =  -fs*amp*b1
	      if(kgamma(nob).eq.1)then
		 gamma = gammai
		 f(5) = f56
		 f(6) = 0
		 f(7) = 0
	      endif
	      if(kgamma(nob).eq.2)then
		 gamma = gammav
		 f(5) = 0
		 f(6) = f56
		 f(7) = 0
	      endif
	      if(kgamma(nob).eq.3)then
		 gamma = gammah
		 f(5) = 0
		 f(6) = 0
		 f(7) = f56
	      endif
	      if(kgamma(nob).eq.4)then
		 gamma = (gammai+gammav)/2
		 f(5) = f56/2
		 f(6) = f56/2
		 f(7) = 0
	      endif
	      dxdtaup = taup/x
	      dxdbetap = betap/x
	      dtaupdt0= -1/te
	      dtaupdte = -tau/te
	      dfdx = fs*(damp*(b0-gamma*b1)+amp*(db0-gamma*db1)/rho)
	      f(1) = dfdx*dxdtaup*dtaupdt0
	      f(2) = dfdx*dxdbetap
	      f(3) = dfdx*dxdtaup*dtaupdte
          f(4) = -fs*amp*db0*x/rho**2
     *               +fs*amp*db1*x/rho**2*gamma
c             tau =  piex*qn + piey*qe
c	      dbeta= -piex*qe + piey*qn
c	      f(8) = dfdx*(dxdtaup*qn - dxdbetap*qe)
c	      f(9) = dfdx*(dxdtaup*qe + dxdbetap*qn)
	      f(8) = dfdx*(dxdtaup*qn + dxdbetap*qe)
	      f(9) = dfdx*(dxdtaup*qe - dxdbetap*qn)
	      f(nparms-2+3*nob) = amp*(b0-gamma*b1)
	      f(nparms-1+3*nob) = 1
	      f(nparms+3*nob) = see(k,nob)
	      flpre = fs*amp*(b0-gamma*b1) + fb + fsee*see(k,nob)
	      y = fl - flpre
	      chi2add = y**2/sig2
c	      write(6,*) k, nob, date(k,nob), y, chi2add
          if(iloop.eq.1) write(61, *) nob, k, t, f(4),
     +	          amp*(b0-gamma*b1), db0, db1
	      if(iloop.eq.1) write(62, *) nob,k,t,f(1),f(2),f(3),f(4),
     +            damp,dfdx,y,sig2
          if(iloop.eq.1) write(63, *) nob, k, t, x, x2, amp, b0, b1

	      do 80 i = 1,ntop
		 d(i) = d(i) + y*f(i)/sig2
		 do 80 j = 1,ntop
		    b(i,j) = b(i,j) + f(i)*f(j)/sig2
 80	      continue
c	if(nob.eq.3)write(6,*)b(14,14),b(14,15),b(15,15),f(14),f(15)
	      if(chi2add.gt.chi2max(nob))then
		 if(nob.ne.1.or.abs(t-t0par).gt.1)then
		    chi2max(nob) = chi2add
		    kmax(nob) = k
		 endif
	      endif
	      chi2(nob) = chi2(nob) + chi2add
 84	   continue
 85	continue
	do 86 i=1,ntop
	   if(kfix(i).eq.1)b(i,i)=b(i,i)+1e16
 86	continue
c	do i = 1, nparms
c	   do j = 1, i
c	      write(6,*) 'b', i, j, b(i,j)
c	   enddo
c	enddo
c	write(6,*)b(14,14),b(14,15),b(15,15)
	call inv(b,c,dum,ntop,ntop)
	do 90 i=1,ntop
	   da(i) = 0
	   do 90 j=1,ntop
	      da(i) = da(i) + c(i,j)*d(j)
 90	continue
	n1 = 9 + (2-1)*3 + 1
	n2 = 9 + (3-1)*3 +1
	ec1 = sqrt(c(n1,n1)/a(n1)**2 - 2*c(n1,n2)/a(n1)/a(n2)
     *                + c(n2,n2)/a(n2)**2)
	col1 = -2.5*log10(a(n2)/a(n1))
	qi1 = 18 - 2.5*log10(a(n1))
	qi1e = sqrt(c(n1,n1))/a(n1)
	n1 = 9 + (2-1)*3 + 2
	n2 = 9 + (3-1)*3 +2
	ec2 = sqrt(c(n1,n1)/a(n1)**2 - 2*c(n1,n2)/a(n1)/a(n2)
     *                + c(n2,n2)/a(n2)**2)
	col2 = -2.5*log10(a(n2)/a(n1))
	qi2 = 18 - 2.5*log10(a(n1))
	qi2e = sqrt(c(n1,n1))/a(n1)
	write(6,89)col1,ec1,qi1,qi1e,col2,ec2,qi2,qi2e,0.888
	c11 = c(8,8)
	c12 = c(8,9)
	c22 = c(9,9)
	det = c11*c22-c12**2
	a1 = a(8)
	a2 = a(9)
	trace = c11 + c22
	ec1 = sqrt((trace + sqrt(trace**2 - 4*det))/2)
	ec2 = sqrt((trace - sqrt(trace**2 - 4*det))/2)
	theta = atan(2*c12/(c11-c22))/2
	write(6,899)a1,a2,c11,c22,c12,ec1,ec2,theta*180/3.14159
 89	format(9f8.4)
 899	format(9f9.5)
	theta = atan(2*c12/(c11-c22))/2
	cth = cos(theta)
	sth = sin(theta)
	a1p = a1*cth + a2*sth
	a2p =-a1*sth + a2*cth
	err1 = sqrt(c11*cth**2 + 2*c12*sth*cth + c22*sth**2)
	err2 = sqrt(c11*sth**2 - 2*c12*sth*cth + c22*cth**2)
c	write(6,*)err1,err2,theta
	sig2 = (c22*a1**2-2*c12*a1*a2+c11*a2**2)/det
c	write(6,96)a1p,a2p,err1,err2,180*theta/3.14159,sqrt(sig2)
c	write(6,*)1
c	write(6,*)err1*err2,sqrt(det)
c	write(6,*)2
	do 91 i=1,ntop
	   s(i) = sqrt(c(i,i))
 91	continue
c	colerr = 13/12.*sqrt(c(13,13)/a(13)**2+c(16,16)/a(16)**2
c     *            -2*c(13,16)/a(13)/a(16))
c	colerr2 = 13/12.*sqrt(c(14,14)/a(14)**2+c(17,17)/a(17)**2
c     *            -2*c(14,17)/a(14)/a(17))
	write(6,*)colerr,colerr2

	write(6,*) 'iloop=', iloop
	if(iloop.le.3) then
	   ifile = 50+iloop
	   write(ifile,*) '# iloop'
	   write(ifile,*) iloop
	   write(ifile,*) '# kfix'
	   write(ifile,*) kfix
	   write(ifile,*) '# fac'
	   write(ifile,*) fac
	   write(ifile,*) '# chi2'
	   write(ifile,*) (chi2(nob), nob=1,nobs)
	   write(ifile,*) '# a'
	   write(ifile,96)a
	   write(ifile,*) '# da'
	   write(ifile,96)da
	   write(ifile,*) '# s'
	   write(ifile,96)s
	   write(ifile,*) '# b'
	   write(ifile,*) b
	   write(ifile,*) '# c'
	   write(ifile,*) c
	   write(ifile,*) '# d'
	   write(ifile,*) d
	   write(ifile,*) '# da'
	   write(ifile,*) da
	   flush(ifile)
	endif

	do 92 i=1,ntop
	   do 92 j=1,ntop
	      c(i,j) = c(i,j)/s(i)/s(j)
 92	continue
	
	write(6,96)a
	write(6,*)
	write(6,96)da
 96	format(f10.5,40f9.5)
	write(6,*)
	write(6,96)s
	write(6,*)
	do 97 i=1,ntop
c	   write(6,96)(c(i,j),j=1,ntop)
 97	continue
	write(6,96)c(8,9)
	chi2tot = 0
	knttot = 0
	do 971 nob=1,nobs
	   chi2tot =chi2tot + chi2(nob)
	   knttot = knttot + knt(nob)
 971	continue
	write(6,*)
	if(nobs.eq.1)then
	   write(6,981)chi2max,kmax,(date(kmax(nob),nob),nob=1,nobs)
	   write(6,991)chi2,chi2tot,knt,knttot
	endif
	if(nobs.eq.2)then
	   write(6,982)chi2max,kmax,(date(kmax(nob),nob),nob=1,nobs)
	   write(6,992)chi2,chi2tot,knt,knttot
	endif
	if(nobs.eq.3)then
	   write(6,983)chi2max,kmax,(date(kmax(nob),nob),nob=1,nobs)
	   write(6,993)chi2,chi2tot,knt,knttot
	endif
	if(nobs.eq.4)then
	   write(6,984)chi2max,kmax,(date(kmax(nob),nob),nob=1,nobs)
	   write(6,994)chi2,chi2tot,knt,knttot
	endif
	if(nobs.eq.5)then
	   write(6,985)chi2max,kmax,(date(kmax(nob),nob),nob=1,nobs)
	   write(6,995)chi2,chi2tot,knt,knttot
	endif
	if(nobs.eq.6)then
	   write(6,986)chi2max,kmax,(date(kmax(nob),nob),nob=1,nobs)
	   write(6,996)chi2,chi2tot,knt,knttot
	endif
	if(nobs.eq.7)then
	   write(6,987)chi2max,kmax,(date(kmax(nob),nob),nob=1,nobs)
	   write(6,997)chi2,chi2tot,knt,knttot
	endif
	if(nobs.eq.8)then
	   write(6,988)chi2max,kmax,(date(kmax(nob),nob),nob=1,nobs)
	   write(6,998)chi2,chi2tot,knt,knttot
	endif
	if(nobs.eq.9)then
	   write(6,989)chi2max,kmax,(date(kmax(nob),nob),nob=1,nobs)
	   write(6,999)chi2,chi2tot,knt,knttot
	endif
c 981	format(1f12.4,1i4,1f10.4)
c 991	format(2f12.4,2i4)
c 982	format(2f12.4,2i4,2f10.4)
c 992	format(3f12.4,3i4)
c 983	format(3f12.4,3i4,3f10.4)
c 993	format(4f12.4,4i4)
c 984	format(4f12.4,4i4,4f10.4)
c 994	format(5f12.4,5i4)
c 985	format(5f12.4,5i4,5f10.4)
c 995	format(6f12.4,6i4)
c 986	format(6f12.4,6i4,6f10.4)
c 996	format(7f12.4,7i4)
c 987	format(7f12.4,7i4,7f10.4)
c 997	format(8f12.4,8i4)
c 988	format(8f12.4,8i4,8f10.4)
c 998	format(9f12.4,9i4)
 981	format(1f12.4/,1i12/,1f12.4/)
 991	format(2f12.4/,2i12)
 982	format(2f12.4/,2i12/,2f12.4/)
 992	format(3f12.4/,3i12)
 983	format(3f12.4/,3i12/,3f12.4/)
 993	format(4f12.4/,4i12)
 984	format(4f12.4/,4i12/,4f12.4/)
 994	format(5f12.4/,5i12)
 985	format(5f12.4/,5i12/,5f12.4/)
 995	format(6f12.4/,6i12)
 986	format(6f12.4/,6i12/,6f12.4/)
 996	format(7f12.4/,7i12)
 987	format(7f12.4/,7i12/,7f12.4/)
 997	format(8f12.4/,8i12)
 988	format(8f12.4/,8i12/,8f12.4/)
 998	format(9f12.4/,9i12)
 989	format(9f12.4/,9i12/,9f12.4/)
 999	format(10f12.4/,10i12)
	facold = fac
	if(kntfac.le.1)then
	   if(nobbad.ne.-1)then
	      go to 200
	   else
	      read(5,*)fac
	   endif
	   if(fac.gt.6)kntfac = fac
	   if(fac.gt.2)fac=facold
	else
	   kntfac = kntfac-1
	endif
	if(fac.lt.-2)go to 200
	do 110 i=1,ntop
	   if(kfix(i).eq.1)go to 110
	   a(i) = a(i) + fac*da(i)
 110	continue
	go to 6
 200	continue
	ifile = 60
	write(ifile,*) '# chi2'
	write(ifile,*) (chi2(nob), nob=1,nobs)
	write(ifile,*) '# a'
	write(ifile,96) a
	write(ifile,*) '# s'
	write(ifile,96) s
	flush(ifile)
	
	if(nobbad.ne.-1)then
	   chi2bad = chi2(nobbad)
	   chi2maxbad = chi2max(nobbad)
	   if(chi2maxbad.lt.sigmax**2*chi2bad/kntgood)go to 300
	   badd(kmax(nobbad),nobbad) = .true.
	 write(6,*)kmax(nobbad),date(kmax(nobbad),nobbad),chi2bad,kntgood
	   if(nobbad.ne.nobmoa)then
	      write(81,301)kmax(nobbad),nobbad,date(kmax(nobbad),nobbad),
     *             chi2maxbad,sigmax,chi2bad/kntgood
 301	      format(6x,'bad(',i5,',',i2,')=.true. !',f10.4,3f8.3)
	   else
	      write(81,302)kmax(nobbad),date(kmax(nobbad),nobbad),
     *             chi2maxbad,sigmax,chi2bad/kntgood
 302	      format(6x,'bad(',i5,',nobmoa)=.true. !',f10.4,3f8.3)
	   endif
	   go to 8
 300	   continue
	   write(6,*)'suggested error factor',sqrt(chi2bad/(kntgood-2))
	   stop
	endif
	nob = 1
	fsogle = a(nparms-2+3*nob)
	fbogle = a(nparms-1+3*nob)
	do 220 nob=1,nobs
	   ifile = nob+36
	   fs = a(nparms-2+3*nob)
	   fb = a(nparms-1+3*nob)
	   fsee = a(nparms+3*nob)
	   if(nob.eq.1.or.nob.eq.2.or.
     *           nob.eq.4.or.nob.eq.5.or.nob.eq.6.or.nob.eq.7)then
	      gamma = a(5)
	   endif
	   if(nob.eq.3)then
	      gamma = a(6)
	   endif
	      if(kgamma(nob).eq.1)then
		 gamma = a(5)
	      endif
	      if(kgamma(nob).eq.2)then
		 gamma = a(6)
	      endif
	      if(kgamma(nob).eq.3)then
		 gamma = a(7)
	      endif
	      if(kgamma(nob).eq.4)then
		 gamma = (a(5)+a(6))/2
	      endif
	   do 210 k =1,knt(nob)
	      if(bad(k,nob))go to 210
	      dat = date(k,nob)
	      fl = flux(k,nob)
	      e = err(k,nob)
c	      if(nob.eq.7)write(6,*)e,1
	      fl = (fl-fb-fsee*see(k,nob))/fs*fsogle + fbogle
	      if(fl.le.0)go to 210
	      e = e/fs*fsogle
c	      if(nob.eq.7)write(6,*)e,2
	      e = e/fl*2.5/log(10.)
c	      if(nob.eq.7)write(6,*)e,3,fl
c	      if(nob.eq.7)read(5,*)xyz
	      fl = 18 - 2.5*log10(fl)
c	      fl = 2.5*log10(fl/(fsogle+fbogle))
	      t = date(k,nob)
	      tau = (t-t0)/te
	      call geta(qn,qe,t,alpha,delta,t0par)
	      call gett(qnp,qep,t,qlat(nob),qlong(nob),alpha,delta)
	      qn = qn + qnp
	      qe = qe + qep
c	      call geta(qn,qe,t)
	      qw = qe*cos(0.2) + qn*sin(0.2)
	      dtau =  piex*qn + piey*qe
	      dbeta= -piex*qe + piey*qn
	      dbeta = -dbeta
	      taup = tau + dtau
	      betap = beta + dbeta
	      x2 = betap**2 + taup**2
	      x = sqrt(x2)
              z = x/rho
	      call 
     *        getb0p(b0,b1,db0,db1,z,b0tab,b1tab,db0tab,db1tab,ntab)      
	      amp = (x2+2)/x/sqrt(x2+4)
	      flpre = fsogle*amp*(b0-gamma*b1) + fbogle
	      flpre =  18 - 2.5*log10(flpre)
	      call geta(qn,qe,t,alpha,delta,t0par)
	      call gett(qnp,qep,t,qlat(1),qlong(1),alpha,delta)
	      qn = qn + qnp
	      qe = qe + qep
	      dtau =  piex*qn + piey*qe
	      dbeta= -piex*qe + piey*qn
	      dbeta = -dbeta
	      taup = tau + dtau
	      betap = beta + dbeta
	      x2 = betap**2 + taup**2
	      x = sqrt(x2)
              z = x/rho
	      call 
     *        getb0p(b0,b1,db0,db1,z,b0tab,b1tab,db0tab,db1tab,ntab)      
	      amp = (x2+2)/x/sqrt(x2+4)
	      flpreogle = fsogle*amp*(b0-gammai*b1) + fbogle
	      flpreogle =  18 - 2.5*log10(flpreogle)
c	      flpre =  2.5*log10(flpre/(fsogle+fbogle))
	      res = (fl-flpre)
	      if(nob.eq.6)write(55,208)dat,z,log10(z)
 208	      format(f10.5,f9.4,f8.4)
	      write(ifile,209)dat,fl,e,res,res/e,flpre
 209	      format(f10.4,3f9.4,f9.3,f9.4,2f8.4)
 210	   continue
 220	continue
	gamma = a(5)
	t = t0 - 2.5*te
	k = 0
 230	continue
	dt = te/100
	if(abs(t-t0)/te.lt.0.1)dt =te/10000
	t = t+dt
	k = k + 1
	if(t.gt.t0+2*te)go to 240
	tau = (t-t0)/te
	call geta(qn,qe,t,alpha,delta,t0par)
	call gett(qnp,qep,t,qlat(nob),qlong(nob),alpha,delta)
	qn = qn + qnp
	qe = qe + qep
c	call geta(qn,qe,t)
	dtau =  piex*qn + piey*qe
	dbeta= -piex*qe + piey*qn
	dbeta = -dbeta
	taup = tau + dtau
	betap = beta + dbeta
	x2 = betap**2 + taup**2
	x = sqrt(x2)
        z = x/rho
	call getb0p(b0,b1,db0,db1,z,b0tab,b1tab,db0tab,db1tab,ntab)      
	amp = (x2+2)/x/sqrt(x2+4)
	fl = fsogle*amp*(b0-gamma*b1) + fbogle
	fl =  18 - 2.5*log10(fl)
	flp=  18 - 2.5*log10(fsogle*amp + fbogle)
c	flp=  2.5*log10((fsogle*amp + fbogle)/(fsogle+fbogle))
 211	format(f10.4,f8.4,f12.4,2f8.4,3f9.4,4f8.4)
	write(35,211)t,fl,z,b1,flp,amp,b0,b1,qn,qe,dtau,dbeta
	if(mod(k,5).eq.0)write(31,212)taup,betap,t
 212	format(2f8.3,f9.2)
c	write(6,*)t,fl,z,b1,flp,amp,b0,b1
c	read(5,*)xyz
	go to 230
 240	continue
        write(6,241)(a(i),i=1,5)
	do 243 j=2,ntop/5 + 1
	   k1 = j*5-4
	   k2 = j*5
	   if(k2.gt.ntop)k2=ntop
        write(6,242)(a(i),i=k1,k2)
 241    format('        data a/',f11.6,',',4(f10.6,','))
 242    format('     *  ',5(f10.6,','))
 243	continue
c        write(6,241)a
c 241    format(f11.6,',',40(f10.6,','))
	do 250 nz=1,700
	   z = nz/100.
	call  getb0p(b0,b1,db0,db1,z,b0tab,b1tab,db0tab,db1tab,ntab)      
	write(34,249)z,b1,b1*0.5*2.5/log(10.)
 249	format(3f8.4)
 250	continue
	close(61)
	close(62)
	deltat = a(3)*sqrt(a(4)**2-a(2)**2)
	write(6,251)a(1)-deltat,a(1)+deltat
 251	format(2f12.5)
	stop
	end
	

c
	subroutine inv(a,ainv,adum,n,nmat)
	real*8 a(nmat,nmat),adum(nmat,nmat),ainv(nmat,nmat)
	real*8 hdum(1000),hinv(1000),Q
	data tol/0.00000001/
	if(n.gt.1000)then
		write(6,*)n,' bigger than 1000'
		stop
	endif
	do 10 i=1,n
	do 10 j=1,n
	adum(i,j) = a(i,j)
	ainv(i,j) = 0
   10	continue
	do 15 i=1,n
	ainv(i,i)=1
   15	continue
	do 100 j=n,1,-1
	if(abs(adum(j,j)).gt.tol) go to 40
	do 20 j1=j-1,1,-1
	if(abs(adum(j,j1)).gt.tol) go to 21
   20	continue
	write(6,*)j, 'degenerate matrix'
	stop
   21	continue
	do 25 i=1,n
	hdum(i) = adum(i,j)
	hinv(i) = ainv(i,j)
   25	continue
	do 30 i=1,n
	adum(i,j) = adum(i,j1)
	ainv(i,j) = ainv(i,j1)
   30	continue
	do 35 i=1,n
	adum(i,j1) = hdum(i)
	ainv(i,j1) = hinv(i)
   35	continue
   40	continue
	q = 1/adum(j,j)
	do 45 i=1,j
	adum(i,j) = adum(i,j)*q
	hdum(i) = -adum(i,j)
   45	continue
	do 50 i=1,n
	ainv(i,j) = ainv(i,j)*q
	hinv(i) = -ainv(i,j)
   50	continue
	do 70 j1 = 1,j-1
	q = adum(j,j1)
	if(q.eq.0.0)go to 69
	do 55 i=1,j
	adum(i,j1) = adum(i,j1) + q*hdum(i)
   55	continue
	do 60 i=1,n
	ainv(i,j1) = ainv(i,j1) + q*hinv(i)
   60	continue
   69	continue
   70	continue
  100	continue
	do 200 j=1,n-1
	do 110 i=1,n
	hinv(i) = ainv(i,j)
  110	continue
	do 130 j1=j+1,n
	q = -adum(j,j1)
	if(q.eq.0.0)go to 129
	do 120 i=1,n
	ainv(i,j1) = ainv(i,j1) + q*hinv(i)
  120	continue
  129	continue
  130	continue
  200	continue
	do 220 i=1,n
	do 220 j=1,n
	adum(i,j)=0
  	do 220 k=1,n
	adum(i,j) = adum(i,j) + a(i,k)*ainv(k,j)
  220	continue
	do 230 i=1,n
	if(abs(adum(i,i)-1).gt.tol)write(6,*)i,i,adum(i,i)
	do 230 j=1,n
	if(i.eq.j)go to 230
	if(abs(adum(i,j)).gt.tol)write(6,*)i,j,adum(i,j)
  230	continue
	return
	end
c
      subroutine mkgvec(epslog,gvec,dgvec,ndiv,epstop,epsbot)
      implicit real*8 (a-h,o-z)
      dimension epslog(ndiv),gvec(ndiv),dgvec(ndiv)
      data pi/3.141592653589793/
      range = log10(epstop/epsbot)
      epsdiv = range/(ndiv-1)
      do 20 i=1,ndiv
         epslg = log10(epsbot) + (i-1)*epsdiv
         epslog(i) = epslg
         eps = 10.**epslg
         mdiv = nint(pi*10/eps)
         if(mdiv.lt.100)mdiv = 100
         if(eps.lt.2)then
            app = 2*log(8/eps)
         else
            app = 2*pi/eps*(1-1/eps**2)
         endif
         tot = 0
         totder=0
         do 10 j=1,mdiv
            theta = (j-0.5)/mdiv*pi
            den = sqrt(eps**2 + 4*(sin(theta/2))**2)
            tot = tot + 2*pi/mdiv/den
            totder=totder - 2*pi/mdiv/den**3*eps
 10      continue
         gvec(i) = tot
         dgvec(i)= totder
c         write(16,11)i,eps,tot,app,tot/app,totder
 11      format(i4,5f12.6)
 20   continue
      return
      end
c
	subroutine getpsi(psi,phi,ecc)
	implicit real*8 (a-h,o-z)
	pi = 3.14159265d0
	psi= phi
	do 10 i=1,4
	   fun = psi - ecc*sin(psi)
	   dif = phi - fun
	   der = 1 - ecc*cos(psi)
	   psi = psi + dif/der
 10	continue
	return
	end

c
	subroutine geta(qn,qe,hjd,alpha,delta,t0)
	implicit real*8 (a-h,o-z)
	real*8 sun(3),xpos(3),ypos(3),rad(3),north(3),east(3)
	real*8 spring(3),summer(3)
	data spring/1.,0.,0./
	data summer/0.,0.9174,0.3971/
	data pi/3.14159265/
	ecc = 0.0167
c	vernal = 2719.0
	vernal = 2719.55
c	vernal = 2719.0 + 1000.
	offset = 75
c	offset = 0
	peri   = vernal - offset
	phi = (1 - offset/365.25)*2*pi
c	phi = (1 - 91/365.25)*2*pi
c	phi = 0
	call getpsi(psi,phi,ecc)
c	write(6,*)psi,phi
	costh = (cos(psi) - ecc)/(1-ecc*cos(psi))
	sinth = -sqrt(1-costh**2)
	do 3 i = 1,3
c	   xpos(i) = spring(i)*cos(psi) + summer(i)*sin(psi)
c	   ypos(i) =-spring(i)*sin(psi) + summer(i)*cos(psi)
	   xpos(i) = spring(i)*costh + summer(i)*sinth
	   ypos(i) =-spring(i)*sinth + summer(i)*costh
 3	continue
c	write(6,4)xpos
c	write(6,4)ypos
c	read(5,*)xyz
 4	format(3f10.4)
	north(1) = 0
	north(2) = 0
	north(3) = 1
	radian = 180/3.14159265
c	alpha = (18 +  05./60 + 40.00/3600.)*15
c	delta = -(32 + 56./60 + 08.6/3600.)
c	t0 = 2771
c	t0 = t0-1
	rad(1) = cos(alpha/radian)*cos(delta/radian)
	rad(2) = sin(alpha/radian)*cos(delta/radian)
	rad(3) = sin(delta/radian)
	call cross(east,north,rad)
	call dot(e2,east,east)
	do 5 i=1,3
	   east(i) = east(i)/sqrt(e2)
 5	continue
	call cross(north,rad,east)
 6	format(3f7.3)
c	theta = (t0+1 - peri)/365.25*360.
	phi   = (t0+1 - peri)/365.25*2*pi
	call getpsi(psi,phi,ecc)
	qn2 = 0
	qe2 = 0
	do 10 i=1,3
c	   sun(i) = xpos(i)*cos(theta/radian)+ypos(i)*sin(theta/radian)
	   sun(i) = xpos(i)*(cos(psi)-ecc) +
     *              ypos(i)*sin(psi)*sqrt(1-ecc**2)
	   qn2 = qn2 + sun(i)*north(i)
	   qe2 = qe2 + sun(i)*east(i)
 10	continue
c	theta = (t0-1 - peri)/365.25*360.
	phi   = (t0-1 - peri)/365.25*2*pi
	call getpsi(psi,phi,ecc)
	qn1 = 0
	qe1 = 0
	do 20 i=1,3
c	   sun(i) = xpos(i)*cos(theta/radian)+ypos(i)*sin(theta/radian)
	   sun(i) = xpos(i)*(cos(psi)-ecc) +
     *              ypos(i)*sin(psi)*sqrt(1-ecc**2)
	   qn1 = qn1 + sun(i)*north(i)
	   qe1 = qe1 + sun(i)*east(i)
 20	continue
c	theta = (t0 - peri)/365.25*360.
	phi   = (t0 - peri)/365.25*2*pi
	call getpsi(psi,phi,ecc)
	qn0 = 0
	qe0 = 0
	do 30 i=1,3
c	   sun(i) = xpos(i)*cos(theta/radian)+ypos(i)*sin(theta/radian)
	   sun(i) = xpos(i)*(cos(psi)-ecc) +
     *              ypos(i)*sin(psi)*sqrt(1-ecc**2)
	   qn0 = qn0 + sun(i)*north(i)
	   qe0 = qe0 + sun(i)*east(i)
 30	continue
	vn0 = (qn2-qn1)/2
	ve0 = (qe2-qe1)/2
	factor = 365.25*4.74
c	write(6,*)qn0,qe0,vn0*factor,ve0*factor
c -0.0927258367 -0.502471007  2.01547714 -25.4080095
c	read(5,*)xyz
	t = hjd
c	theta = (t - peri)/365.25*360.
	phi   = (t - peri)/365.25*2*pi
	call getpsi(psi,phi,ecc)
	qn = -qn0 - vn0*(t-t0)
	qe = -qe0 - ve0*(t-t0)
	do 40 i=1,3
c	   sun(i) = xpos(i)*cos(theta/radian)+ypos(i)*sin(theta/radian)
	   sun(i) = xpos(i)*(cos(psi)-ecc) +
     *              ypos(i)*sin(psi)*sqrt(1-ecc**2)
	   qn = qn + sun(i)*north(i)
	   qe = qe + sun(i)*east(i)
 40	continue
c	write(31,11)njd,qn,qe
 11	format(i6,2f9.5)
 100	continue
	return
	end
c       
	subroutine cross(c,a,b)
	implicit real*8 (a-h,o-z)
	dimension a(3),b(3),c(3)
	c(1) = a(2)*b(3) - b(2)*a(3)
	c(2) = a(3)*b(1) - b(3)*a(1)
	c(3) = a(1)*b(2) - b(1)*a(2)
	return
	end
c
      subroutine dot(c,a,b)
      implicit real*8 (a-h,o-z)
      dimension a(3),b(3)
      c = a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
      return
      end

C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.

      FUNCTION rc(x,y)
      REAL*8 rc,x,y,ERRTOL,TINY,SQRTNY,BIG,TNBG,COMP1,COMP2,THIRD,C1,C2,
     *C3,C4
      PARAMETER (ERRTOL=.04d0,TINY=1.69d-38,SQRTNY=1.3d-19,BIG=3d37,
     *TNBG=TINY*BIG,COMP1=2.236d0/SQRTNY,COMP2=TNBG*TNBG/25d0,THIRD=
     *1d0/3d0,C1=.3d0,C2=1d0/7d0,C3=.375d0,C4=9d0/22d0)
      REAL*8 alamb,ave,s,w,xt,yt
      if(x.lt.0d0.or.y.eq.0d0.or.(x+dabs(y)).lt.TINY.or.(x+
     *dabs(y)).gt.BIG.or.(y.lt.-COMP1.and.x.gt.0d0.and.x.lt.COMP2))pause 
     *'invalid arguments in rc'
      if(y.gt.0d0)then
        xt=x
        yt=y
        w=1d0
      else
        xt=x-y
        yt=-y
        w=dsqrt(x)/dsqrt(xt)
      endif
1     continue
        alamb=2d0*dsqrt(xt)*dsqrt(yt)+yt
        xt=.25d0*(xt+alamb)
        yt=.25d0*(yt+alamb)
        ave=THIRD*(xt+yt+yt)
        s=(yt-ave)/ave
      if(dabs(s).gt.ERRTOL)goto 1
      rc=w*(1d0+s*s*(C1+s*(C2+s*(C3+s*C4))))/dsqrt(ave)
      return
      END
c
	subroutine 
     *    getb0p(b0,b1,db0,db1,z,b0tab,b1tab,db0tab,db1tab,ntab)      
	implicit real*8 (a-h,o-z)
	dimension b0tab(ntab),b1tab(ntab),db0tab(ntab),db1tab(ntab)
	if(z.lt.0.001)then
	   b0 = 2*z
	   db0= 2
	   b1 = -5./14.*z
	   db1 = -5./14
	   return
	endif
	if(z.gt.10)then
	   b0 = 1 + 1/8./z**2
	   db0= -1/4./z**3
	   b1 = 1./40./z**2
	   db1 = -1./20./z**3
	   return
	endif
	k = nint(1000*z-0.5)
	if(k.lt.1)k=1
	if(k.gt.ntab)k=ntab
	x = z*1000 - k
	b0 = x*b0tab(k+1) +(1-x)*b0tab(k)
	db0 = x*db0tab(k+1) +(1-x)*db0tab(k)
	b1 = x*b1tab(k+1) +(1-x)*b1tab(k)
	db1 = x*db1tab(k+1) +(1-x)*db1tab(k)
	return
	end
c
	subroutine hjdcor(hjd,date,alpha,delta)
	implicit real*8 (a-h,o-z)
	real*8 sun(3),xpos(3),ypos(3),rad(3)
	real*8 spring(3),summer(3)
	data spring/1.,0.,0./
	data summer/0.,0.9174,0.3971/
	data pi/3.14159265d0/
	radian = 180/pi
	ecc = 0.0167
	vernal = 2719.55
	offset = 75
	peri   = vernal - offset
	phi = (1 - offset/365.25)*2*pi
	call getpsi(psi,phi,ecc)
	costh = (cos(psi) - ecc)/(1-ecc*cos(psi))
	sinth = -sqrt(1-costh**2)
	do 3 i = 1,3
	   xpos(i) = spring(i)*costh + summer(i)*sinth
	   ypos(i) =-spring(i)*sinth + summer(i)*costh
 3	continue
	rad(1) = cos(alpha/radian)*cos(delta/radian)
	rad(2) = sin(alpha/radian)*cos(delta/radian)
	rad(3) = sin(delta/radian)
	phi   = (date - peri)/365.25*2*pi
	call getpsi(psi,phi,ecc)
	qr0 = 0
	do 30 i=1,3
	   sun(i) = xpos(i)*(cos(psi)-ecc) +
     *              ypos(i)*sin(psi)*sqrt(1-ecc**2)
	   qr0 = qr0 - sun(i)*rad(i)*8.31/1440.
 30	continue
	hjd = date + qr0
 100	continue
	return
	end
c
c       to de-implement terrestrial parallax: change reearth (/100)
c       
c
	subroutine gett(qn,qe,hjd,qlat,qlong,alpha,delta)
	implicit real*8 (a-h,o-z)
	real*8 sun(3),xpos(3),ypos(3),rad(3),north(3),east(3)
	real*8 spring(3),summer(3)
	data pi/3.14159265d0/
	radian = 180/pi
        rearth = 20000/pi/1.5e8 !/100   (remove "/100" to implement)
	vernal = 2719.0
	north(1) = 0
	north(2) = 0
	north(3) = 1
c	alpha = (18 +  05./60 + 41.05/3600.)*15
c	delta = -(28 + 45./60 + 36.2/3600.)
c	t0 = 4233.66
	rad(1) = cos(alpha/radian)*cos(delta/radian)
	rad(2) = sin(alpha/radian)*cos(delta/radian)
	rad(3) = sin(delta/radian)
	call cross(east,north,rad)
	call dot(e2,east,east)
	do 5 i=1,3
	   east(i) = east(i)/sqrt(e2)
 5	continue
	call cross(north,rad,east)
c
	qn = 0
	qe = 0
	qr = 0
        phase = (hjd-vernal)*366.25/365.25 + qlong/360
c        write(6,*)phase
        sun(1) = -cos(phase*2*pi)*cos(qlat/radian)
        sun(2) = -sin(phase*2*pi)*cos(qlat/radian)
        sun(3) = -sin(qlat/radian)
	do 30 i=1,3
	   qn = qn + sun(i)*north(i)*rearth
	   qe = qe + sun(i)*east(i)*rearth
	   qr = qr + sun(i)*rad(i)*rearth
 30	continue
	return
	end
c
c
	subroutine getlatlon(qlat,qlong,nobs)
	implicit real*8 (a-h,o-z)
	dimension qlat(nobs),qlong(nobs)
	do 10 nob=1,nobs
c ---   ogle
	if(nob.eq.1)then    
	   qlong(nob) = -70.702
	   qlat(nob) = -29.0083
	endif
c ---   moa 170 27.9East 43 59.2south  (MOA)
	if(nob.eq.4)then
	   qlong(nob) = +(170 + 27.9/60.)
	   qlat(nob) =  -(43 + 59.2/60.)
	endif
c ---   ctio
	if(nob.eq.2.or.nob.eq.3)then
	   qlong(nob) =  -70.805
	   qlat(nob)  =  -30.168
	endif
c ---   TASMANIA
	if(nob.eq.5)then
	   qlong(nob) =  +(147 +26/60. + 21/3600.)
	   qlat(nob)  =  -(42 +48/60. +18/3600.)
	endif
c ---   Mt Lemmon
	if(nob.eq.6)then
	   qlong(nob) =  -110.789
	   qlat(nob)  =  +32.443
	endif
cc ---   farmcove 174:53:37E 36:53:37S  ???
c ---   lt
	if(nob.eq.7)then
c	   qlong(nob) = +(174 + 53/60. + 37/3600.)
c	   qlat(nob) =   -(36 + 53/60. + 43/3600.)
	   qlong(nob) = -(17  + 52/60. + 45/3600.)
	   qlat(nob) =  +(28 + 45/60. + 44.8/3600.)
	endif
c ---   bronberg
	if(nob.eq.8)then
	   qlong(nob) = +(28 + 26/60. + 18/3600.)
	   qlat(nob)  = -(25 + 54/60. + 32/3600.)
	endif
 10	continue
	return
	end
c

c     
      subroutine getrootband(qlat,qlong,bandname,rootname,filename)
      implicit real*8 (a-h,o-z)
      character*1 filechar(21)
      character*21 filename
      character*3  rootname
      character*1  bandname
      character*8 word
      logical found

      qlong = -70.702
      qlat  = -29.0083
      bandname = 'I'
      return	

c       Skip everything else
      read(filename,1)filechar
 1    format(21a1)
      do 10 i=20,2,-1
         if(filechar(i).eq.'.')go to 12
 10   continue
      write(6,*)'bad filename ',filename
      stop
 12   continue
      bandname = filechar(i-1)
      write(rootname,13)filechar(i+1),filechar(i+2),filechar(i+3)
 13   format(3a1)
c     SURVEY TEAMS
c  --   kmtc (ctio)
      word = 'KMTC'
      length=4
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = -70.815
         qlat  = -30.165
         bandname = 'I'
         return
      endif	
c  --   kmta (SSO)
      word = 'KMTA'
      length=4
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = 149.050
         qlat =  -31.267
         bandname = 'I'
         return
      endif	
c  --   kmts (saao)
      word = 'KMTS'
      length=4
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = +20.8105
         qlat  = -32.3783
         bandname = 'I'
         return
      endif	
c ---   ogle
      word = 'phot.dat'
      length = 8
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = -70.702
         qlat  = -29.0083
         bandname = 'I'
         return
      endif
c ---   moa
      word = 'moa'
      length = 3
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = +(170 + 27.9/60.)
         qlat  = -(43 + 59.2/60.)
         bandname = 'I'
         return
      endif
      if(rootname.ne.'pho'.and.rootname.ne.'dia')go to 20
c     uFUN normal
c ---   ctio
      word = 'CT13'
      length = 4
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = -70.815
         qlat  = -30.165
         return
      endif
c ---   Mt Lemmon
      word = 'LOAO'
      length = 4
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = -110.76
         qlat  = +32.43
         return
      endif
c ---  Farmcove
      word = 'FCO'
      length = 3
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = +(174 + 53/60. + 37/3600.)
         qlat  =  -(36 + 53/60. + 43/3600.)
         return
      endif
c ---  Auckland
      word = 'AO'
      length = 2
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = +(174 + 46/60. + 37/3600.)
         qlat  =  -(36 + 54/60. + 22/3600.)
         return
      endif
c ---   CAO
      word = 'CAO'
      length = 3
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = -(68 + 10/60. + 49/3600.)
         qlat  =  -(22 +57/60. + 10/3600.)
         return
      endif
c ---   Bronberg
      word = 'Bron'
      length = 4
      call vetword(found,word,length,filechar,0)
 8    continue
      if(found)then
         qlong = +(28 + 26/60. + 44/3600.)
         qlat  = -(25 + 54/60. + 48/3600.)
         return
      endif
c ---   Hunters Hill
      word = 'HHO'
      length = 3
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = +(149 + 06/60. + 36/3600.)
         qlat  = -(35 + 09/60. + 45/3600.)
         return
      endif
c ---   Kumeu
      word = 'Kumeu'
      length = 5
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = +(174 + 31/60. + 29/3600.)
         qlat  = -(36 + 48/60. + 23/3600.)
         return
      endif
c ---   Hereford
      word = 'HAO'
      length = 3
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = -(110 + 14/60. + 15/3600.)
         qlat  = +(31 + 27/60. + 08/3600.)
         return
      endif
c ---   Molehile
      word = 'MAO'
      length = 3
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = +(174 + 43/60. + 43/3600.)
         qlat  = -(36 + 47/60. + 55/3600.)
         return
      endif
c ---   MDM
      word = 'MDM'
      length = 3
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = -(111 + 36/60. + 56/3600.)
         qlat  = +(31 + 57/60. + 05/3600.)
         return
      endif
c ---   Palomar
      word = 'PAL'
      length = 3
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = -(116 + 51/60. + 36/3600.)
         qlat  = +(33 + 21/60. + 26/3600.)
         return
      endif
c ---   Perth
      word = 'Perth'
      length = 5
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = +(115 + 49/60. + 00/3600.)
         qlat  = -(31 + 58/60. + 00/3600.)
         return
      endif
c ---   Possum
      word = 'Pos'
      length = 3
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = +(177 + 53/60. + 29/3600.)
         qlat  = -(38 + 37/60. + 26/3600.)
         return
      endif
c ---   Southern Stars
      word = 'SSO'
      length = 3
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = -(149 + 35/60. + 15/3600.)
         qlat  = -(17 + 33/60. + 04/3600.)
         return
      endif
c ---   Vintage Lane
      word = 'VLO'
      length = 3
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = +(173 + 50/60. + 21/3600.)
         qlat  = -(41 + 29/60. + 30/3600.)
         return
      endif
c ---   Wise
      word = 'Wise'
      length = 3
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = +(34 + 45/60. + 44/3600.)
         qlat  = +(30 + 35/60. + 50/3600.)
         return
      endif
      go to 120
 20   continue
c ---  Bronberg (alternate)      
      word = 'bron'
      length = 4
      call vetword(found,word,length,filechar,1)
      if(found)then
         bandname='U'
         go to 8
      endif
      if(rootname.ne.'raw'.and.rootname.ne.'pys')go to 40
c ---   Tasmania
      word = 'U'
      length = 1
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = +(147 + 26/60. + 21/3600.)
         qlat  = -(42 + 48/60. + 18/3600.)
         return
      endif
c ---   SAAO
      word = 'A'
      length = 1
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = +20.8105
         qlat  = -32.3783
         return
      endif
c ---   FTN
      word = 'H'
      length = 1
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = -(156 + 15/60. + 21/3600.)
         qlat  = +(20 + 42/60. + 27/3600.)
         return
      endif
c ---   Liverpool (cararies
      word = 'L'
      length = 1
      call vetword(found,word,length,filechar,0)
      if(found)then
         qlong = -(17 + 53/60. + 45/3600.)
         qlat  = +(28 + 45/60. + 45/3600.)
         return
      endif
 40   continue
      go to 120
 120  continue
      write(6,121)filechar
 121  format(21a1,1x,'lat-long not found')
c      read(5,*)xyz
      qlat = 0
      qlong= 0
      return
      end
c
      subroutine vetword(found,word,length,filechar,iloop)
      character*8 word
      character*1 wordchar(8),filechar(21)
      logical found
      found = .false.
      read(word,1)wordchar
 1    format(8a1)
      do 10 i=0,(21-length)*iloop
         do 5 j=1,length
            if(wordchar(j).ne.filechar(i+j))go to 10
 5       continue
         found = .true.
         return
 10   continue
      return
      end
c
      subroutine vetword2(found,word,length,linechar,iloop)
      character*8 word
      character*1 wordchar(8),linechar(80)
      logical found
      found = .false.
      read(word,1)wordchar
 1    format(8a1)
      do 10 i=0,(80-length)*iloop
         do 5 j=1,length
            if(wordchar(j).ne.linechar(i+j))go to 10
 5       continue
         found = .true.
         return
 10   continue
      return
      end
