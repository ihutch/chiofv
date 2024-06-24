! chigauss is aimed to be a stripped down version of chiofv (but not much).
! With the added facility of giving electron f(v) as sum of gaussians. 

! Velocities are input and output in the same sqrt(T0/m_i) units.
! The range of velocities plotted is 2vw+vdmax-vdmin, where vw is 
! propto sqrt(max(gpar(3,1),gpar(4,1)). So first Gaussian sets plot width.
      program chiofv
      integer nv,npts
      real vw
! Number of velocities in integrals nv, number of Re(vp) points npts
! These are two different velocity arrays.
! Number of Im(vp) points: nim. 
      parameter (nv=400,npts=200,nimax=60,nim2=2*nimax+1,nmmax=-nimax)
! Input to the analysis routines on nv-mesh distrib, velocity, f^prime.
      real f(0:nv+2),v(0:nv+2),ww(0:nv+2)
! Storage for extra wider plot of fe(v) 
      real fep(1:nv),vep(1:nv)
! Intermediate analysis output on npts mesh
      real vp(npts),fv(npts),fvex(npts),fvcomb(npts)
! Final chi(vp) for complex grid.
      real cvreal(npts,-nimax:nimax),cvimag(npts,-nimax:nimax)
      real ceimag(npts,-nimax:nimax),cereal(npts,-nimax:nimax)
! Contouring work array:
      character cworka(npts,-nimax:nimax)
      character ilb
! The velocity step and start.
      real dv,v0,v1
      real Te
      real rmitome
!
! Things describing the input distribution as sum of Gaussians. 
      integer npar
      parameter (npar=4,ngmax=6)
! Parameters; density, velocity shift, vtperp, vtparallel.      
      real gpar(npar,ngmax),gepar(npar,ngmax)
      real fstat(ngmax)
      real theta
      complex cvint,vpec,chihat,chie
      external cvint,chihat
      integer nim,nmm,nvin
      logical laspect,lcolor,lextra,lgrowth,lcoptic,ltwotone,lthresh
      logical lcombplot,lelectrons
      real vpimax,vpi(-nimax:nimax)
      character*50 cfilename

      data nim/nimax/nmm/nmmax/laspect/.false./lcolor/.true./
      data lextra/.false./lgrowth/.true./ltwotone/.true./
      data lcombplot/.false./lelectrons/.true./
      data vpimax/.4/vpi/nim2*0/
      common/ilabel/ilb

! Default: Plot to screen
      call pfset(3)
      call charsize(.02,.02)
      ilb='i'
      fflat=0.
! Silence warnings.
      cworka(1,1)=char(0)
! Extent of integration factor
      vwfac=6.
! Default gaussian data:
      ng=1
      gpar(:,1)=[1.,0.,1.,1.]
      nge=1
      gepar(:,1)=[1.,0.,1.,1.]
! And other defaults
      theta=0.
      Te=1.
      rmitome=1836
      vdmin=0.01
      vdmax=0.01
      vw=vwfac*sqrt(max(gpar(3,1),gpar(4,1)))
      vrgfac=1.
! Default electron landau damping weight is what a maxwellian of
! temperature Te will give.
! Needs to account for the mass ratio 1836 by default.
      eldweight=sqrt(3.141593/(2.*rmitome))
! Smoothing of analytic continuation default setting.
      smfac=1.
! Noise amplitude:
      amp=.00
! Default no coptic data
      lcoptic=.false.
!---------------------------------------------------------------
! Parse command line for changes to defaults
      call parsecmdline(gpar,npar,ngmax,ng,vdmin,vdmax,vw,Te,theta       &
     &     ,vwfac,vpimax,nim,nmm,nimax,smfac,laspect,lcolor,eldweight    &
     &     ,vflat,wflat,fflat,lextra,amp,nw,lgrowth,lcoptic,cfilename    &
     &     ,ltwotone,vrgfac,lthresh,rmitome,omegac,lcombplot,gepar,nge)
!---------------------------------------------------------------
      if(.not.lcoptic)then
! Construct the test distribution function as a sum of Gaussians.
! Scale to sound-speed units assuming inputs relative to Ti=1:
         if(vdmax.eq.0.01)vdmax=-100.
         if(vdmin.eq.0.01)vdmin=100.
         if(Te.eq.99.)then  ! No electrons.
            Te=1.
            lelectrons=.false.
            ltwotone=.false.
            lcombplot=.true.
            ilb='e'
 ! Avoid label clash.
            vpimax=(real(int(vpimax*10.))+.5)/10.
         endif
! Decide vp-range-width vw.
         ct=cos(theta)
         do i=1,ng
            gpar(2,i)=gpar(2,i)/sqrt(Te)
            gpar(3,i)=gpar(3,i)/(Te)
            gpar(4,i)=gpar(4,i)/(Te)
            if(gpar(2,i)*ct.gt.vdmax)vdmax=gpar(2,i)*ct
            if(gpar(2,i)*ct.lt.vdmin)vdmin=gpar(2,i)*ct
         enddo
!         vw=vwfac*sqrt(max(gpar(3,1),gpar(4,1)))
         vw=vwfac*sqrt(maxval(gpar(3,1:npar)))
         dv=(vdmax-vdmin+2.*vw)/nv
         v0=vdmin-vw
!         write(*,*)'vwfac,vdmin,vdmax',vwfac,vdmin,vdmax

! Calculate total ion density
         pne=0.
         do i=1,ng
            pne=pne+gpar(1,i)
         enddo
! Set the distribution function to be sum of maxwellians.
!      write(*,*)'dv,v0',dv,v0
         call fvinit(nv,dv,v0,f,ng,gpar,theta)
! Normalize it and initialize velocity.
         do i=0,nv+2
            f(i)=f(i)/pne
            v(i)=(v0+i*dv)
         enddo
         
         if(fflat.ne.0.)then
            iv=int((vflat-v0)/dv)
            iw=int(wflat/dv)
            call flatspot(nv,dv,f,iv,iw,fflat)
            write(*,'(a,3f6.2)')'Flatspot vflat,wflat,fflat'             &
     &           ,vflat,wflat,fflat
         endif
         call COPTICverif(nv,dv,v0,f)
! Distribution created.
      else
! External data from coptic
         open(13,file=cfilename, form='formatted',status='old',err=101)
         read(13,*) nvin
         read(13,*) dv
         read(13,*) v1
         v0=v1-dv
         if(nvin.gt.nv)then
           write(*,*) 'Input f longer than expected.','nvin=',nvin
           call exit(1)
         endif
! There's a big problem here. nv ought to be in the file.
! Also the velocity has to be a uniform grid. Really only v0 and dv
! ought to be specified.
         read(13,'(1g14.6)',end=102)(f(i),i=1,nvin)
! Also there must not be non-zero values at the ends. So we ought not
! to read the ends we should just put them to zero.
! Also it ought to be normalized. ...
         f(0)=0.
         f(1)=0.
! Set the end of the distribution function to zero.
         do i=nvin+1,nv+2
            f(i)=0.
         enddo
         vavei=0
         fwt=0
      endif  ! End of ion distribution initialization. 
! Find the average ion velocity.
      do i=0,nv+2
         v(i)=(v0+i*dv)
         vavei=vavei+f(i)*v(i)
         fwt=fwt+f(i)
      enddo
      vavei=vavei/fwt
! Increment Smoothing of analytic continuation based on ranges.
      smfac=smfac+ 10.*(npts/(3.*nv))*dv/vpimax

! Find electron susceptibility weight for propagation angle theta
      if(cos(theta).gt.0.)then 
         eldweight=eldweight/cos(theta)
      else
         write(*,*)'CosTheta is non-positive',theta,cos(theta)           &
     &        ,' no angle factors.'
      endif
!-----------------------------------------------------------------
! Perhaps add noise to ion distribution ?
      call addnoise(nv,f,amp)
! Perhaps smooth it :
      if(nw.gt.0)call trismooth(nv-3,nw,f(3),ww,f(3))
!-----------------------------------------------------------------
! Prevent using twotone plot of chi_i if ion shift is large
      vbig=10.
      if(ltwotone.and.(abs(v(0)).gt.vbig.or.abs(v(nv)).gt.vbig))then
         write(*,*)'WARNING: LARGE SHIFT. ltwotone set to false'
         ltwotone=.false.
      endif
!-----------------------------------------------------------------
      if(lthresh)then
! I don't recall what this code section is for. But it often terminates
! the execution. 
! Calculate the chi along just the real axis for thresholds
         call stationary(nv,f,ngmax,fstat,nstat)
         if(nstat.eq.3)then
! Find fractional depth assuming that the middle one is the local minimum.
            fmax=max(fstat(1),fstat(3))
            fmid=min(fstat(1),fstat(3))
            depth=(fmid-fstat(2))/fmid
            write(*,*)nstat,(fstat(i),i=1,nstat),depth
         else
            write(*,*)'Incorrect number of stationaries',nstat
         endif
         call chiion(v0,dv,f,nv,fv,vp,npts,fmax                       &
     &     ,nimax,0,0,vpi,vpimax,smfac,cvreal,cvimag)
!         write(*,'(i5,3f8.4)')(i,vp(i),cvreal(i,0),cvimag(i,0),i=1,npts)
         do i=npts/10,npts-1  ! Test if this calculation makes sense.
            if(cvimag(i,0)*cvimag(i+1,0).le.0.and.cvreal(i,0).le.0)then
               ri=cvimag(i,0)/(cvimag(i,0)-cvimag(i+1,0))
               cvr0=ri*cvreal(i+1,0)+(1-ri)*cvreal(i,0)
               cvi0=ri*cvimag(i+1,0)+(1-ri)*cvimag(i,0)
               vp0=ri*vp(i+1)+(1-ri)*vp(i)
               write(*,*)depth,abs(gpar(2,1)-gpar(2,2))/2.,-Te/cvr0,     &
     &              gpar(1,2)/(gpar(1,1)+gpar(1,2)),gpar(3,2)/gpar(3,1)
               stop
            endif
         enddo
      endif
!-----------------------------------------------------------------
! Calculate chi-ion from it. This is the main calculation of chi arrays
      call chiion(v0,dv,f,nv,fv,vp,npts,fmax                          &
     &     ,nimax,nim,nmm,vpi,vpimax,smfac,cvreal,cvimag)
! This is where we should put the multigaussian ion chi calculation.


      vte2=2.*rmitome !The assumption is that Te=T0 giving v-units.
      fvexmax=-1.e-6
      if(lelectrons)then
! Get chi electron and add on real and imag parts for original plot.
         ct=cos(theta)
         if(ct.gt..01)then
            vefactor=1./(ct*sqrt(vte2))
         else
            write(*,*)'Dangerously low cos theta',ct,theta
            stop
         endif
! Also construct the effective electron distribution fvex.
!         write(*,*)'nge=',nge
         fvex=0.
         do j=1,npts
            do k=nmm,nim
               chie=0.
               do l=1,nge  ! Multigaussian calculation
                  vpec=cmplx(vp(j)-gepar(2,l),vpi(k))*gepar(1,l)     &
                       &    /(ct*sqrt(2.*rmitome*gepar(4,l)))
                  chie=chie+chihat(vpec)
               enddo
               if(.not.ltwotone)then
                  cvimag(j,k)=cvimag(j,k)+imag(chie)
                  cvreal(j,k)=cvreal(j,k)+real(chie)
               endif
               cereal(j,k)=real(chie)
               ceimag(j,k)=imag(chie)
            enddo
            do l=1,nge
               fvex(j)=fvex(j)+gepar(1,l)*rmitome/sqrt(3.141593*2.*rmitome)&
               & *exp(-((vp(j)-gepar(2,l))/ct)**2/(2.*rmitome*gepar(4,l)))
            enddo
         enddo
      endif
! Calculate the effective distribution function shape.
      do j=1,npts
         fvcomb(j)=fv(j)+fvex(j)
      enddo
!---------------------------------------------------------------

! Usually unused plots:
      if(lextra)then
! Plot the wide electron distribution.
         vgtmax=maxval(gepar(4,1:nge))*sqrt(rmitome)
         vgemin=minval(gepar(2,1:nge))-2*vgtmax
         vgemax=maxval(gepar(2,1:nge))+2*vgtmax
         do j=1,nv
            vep(j)=vgemin+(vgemax-vgemin)*(j-1)/(nv-1)
            do l=1,nge
               fep(j)=fep(j)+gepar(1,l)*rmitome/sqrt(3.141593*2.*rmitome)&
               & *exp(-((vep(j)-gepar(2,l))/ct)**2/(2.*rmitome*gepar(4,l)))
            enddo
         enddo
         call autoplot(vep,fep,nv)
         call axlabels('v/(T!d0!d/m!di!d)','(m!di!d/m!de!d)f!de!d')
         call pltend
! Penrose style real vp plots.
         call realvpplot(vp,cvreal(1,0),cvimag(1,0),npts,nimax+1)
         call argdiagplot(cvreal(1,0),cvimag(1,0),npts,Te)
!         call vertslice(vpi(-nim),npts,nim, cvreal(1,-nim))
      endif

! Contour the susceptibility for stability analysis
      call multiframe(2,1,0)
      call ticnumset(7)
      koff=int(min(max(1.-vrgfac,0.)*npts/2.,npts-1.))
      nrmin=1+koff
      nrmax=npts-koff
      if(lcombplot)then            ! Combined plot
         call minmax(fvcomb(nrmin),nrmax-nrmin+1,fvmin,fvmax)
         call pltinit(vp(nrmin),vp(nrmax),fvmin,fvmax)
         call axis
         if(.not.lelectrons)then
            call axlabels('','f!de!d')
         else
            call axlabels('','f!di!d+f!de!d(m!di!d/m!de!d)')
         endif
         if(lelectrons)then
            call dashset(1)
            call color(3)
            call legendline(.02,.95,0,' f!de!d(m!di!d/m!de!d)')
            call polyline(vp(nrmin),fvex(nrmin),nrmax-nrmin+1)
            call dashset(0)
            call color(15)
         endif
         call polyline(vp(nrmin),fvcomb(nrmin),nrmax-nrmin+1)
      else ! Older plot options
         call plotfofv(vp(nrmin),fv(nrmin),nrmax-nrmin+1,fmax)
         call dashset(2)
         call color(3)
         ff=-fmax
         fvexmax=maxval(fvex)
         if((fvex(nrmin)-fvexmax).lt.ff                             &
              .or.(fvex(nrmax)-fvexmax).lt.ff.or.fvexmax.lt.0)then 
!  Do not use suppressed zero
            write(*,*)fvex(1),fvex(npts),fmax,fvexmax
            do i=1,npts
               fvex(i)=(fmax/fvexmax)*fvex(i)
            enddo
            call legendline(.5,.95,0,                                    &
     &           ' (f!dimax!d/f!demax!d).f!de!d(v)')
         else      !  Suppressed zero
            do i=1,npts
               fvex(i)=(fvex(i)-fvexmax)+fmax
               if(fvex(i).gt.fvexmax)fvexmax=fvex(i)
            enddo
            call legendline(.5,.95,0,                                    &
     &           ' (f!de!d(v)-f!de!d(0))m!di!d/m!de!d+f!dimax!d')
         endif
         call polyline(vp(nrmin),fvex(nrmin),nrmax-nrmin+1)
         call dashset(0)
         call color(15)
      endif
      call chicontour(vp,vpi(nmm),npts,nrmin,nrmax,nim,nmm,cvreal(1,nmm) &
     &     ,cvimag(1,nmm),ceimag(1,nmm),cworka,laspect,lcolor,ltwotone)
      call multiframe(0,0,0)
! Make the cvreal the total chi if it is not already (for acpplot).
      if(ltwotone)cvreal=cvreal+cereal
      if(ltwotone)cvimag=cvimag+ceimag
! Calculate parameters along contours and plot:
      if(lgrowth)call acpplot(vp,vpi(nmm),vavei,npts,nim,nmm,cvreal(1    &
     &     ,nmm),cvimag(1,nmm),cworka,laspect,lcolor)
      call exit(0)
 101  write(*,*)'Could not open coptic file',cfilename
      call exit(1)
 102  write(*,*)'Premature coptic file end'
      end
! End of main program.
!*****************************************************************
!*****************************************************************
! Calculate the real and imaginary parts of chi for a specified
! ion distribution function. Ion chi only, here.
      subroutine chiion(v0,dv,f,nv, fv,vp,npts,fmax                      &
     &     ,nimax,nim,nmm,vpi,vpimax,smfac,cvreal,cvimag)

! On entry: v0, dv are the origin and step of the velocity for f(v),
!           f(v) is the distribution of length 0:nv+2.
!           nimax defines the imaginary velocity allocated length
!           nim,nmm define its actual array length.
!           vpimax defines the imag-velocity range,
!           smfac the smoothing factor.
!
! On exit:
!          Arrays of real velocity of length npts are in vp, 
!          fv is the distribution, 
!          vpi contains the imaginary velocity array. 
!          cvreal, cvimag contain the real and imaginary part of
!          chi.(k\lambda_{De})^2.
!
! Input to the analysis routines on nv-mesh distrib, velocity, f^prime.
      real f(0:nv+2)
!      real v(0:nv+2) is not used, uniform grid v0,dv instead
! Intermediate analysis output on npts mesh
      real vp(npts),fv(npts)
      real cvreal(npts,-nimax:nimax),cvimag(npts,-nimax:nimax)

! Fourier transform variables
      parameter (nptsmax=200)
      parameter (efold2=10.)
      integer lensav,nfftwork
      complex ftfp(nptsmax),ftbk(nptsmax)
! Has to be at least 2*N + INT(LOG(REAL(N))/LOG(2.)) + 4. , (N=npts)
      parameter (lensav=3*nptsmax,nfftwork=2*nptsmax)
      real wsave(lensav),fftwork(nfftwork)
      real vpimax,vpi(-nimax:nimax)
      complex cvint,cvi,vpc
      complex chihat
!      logical lchaprox,lchieinc
      external cvint,chihat

!      write(*,*)nv,npts,nim,nmm,nimax,v0

      if(npts.gt.nptsmax)stop 'npts .gt. nptsmax'
!-----------------------------------------------------------------
      fmax=0.
      do j=1,npts
         vp(j)=v0+dv+dv*(nv-1)*(j-1.)/(npts-1.)
         vi=1.+(nv-1)*(j-1.)/(npts-1.)
         i=int(vi)
         xi=vi-i
         if(i.gt.nv)write(*,*)'Too large i',i,j
! Interpolated values of f(v) and f'(v)
         fv(j)=(f(i)*(1.-xi)+f(i+1)*xi)
         if(fv(j).gt.fmax)fmax=fv(j)
! Derivative indices
         ixd=nint(xi)
         id=i+ixd
         xid=xi+0.5-ixd
!--------------------
! Real values of vp.
! For a range of vp's, calculate the pvintegral and imaginary part.
         cvimag(j,0)=-3.1415926*(                                        &
     &        (f(id)-f(id-1))/(dv)*(1.-xid)                              &
     &        +(f(id+1)-f(id))/(dv)*xid)
! The ion principal value integral
         cvreal(j,0)=-pvint(v0,dv,f,nv,vp(j))
!--------------------
! Positive values of imaginary part of vp.
! Integration along the real axis gives chi integral.
         do k=1,nim
            vpi(k)=vpimax*k/nim
! Complex v_p.
            vpc=cmplx(vp(j),vpi(k))
! Integral (which is chi.(k\lambda_{De})^2
            cvi=-cvint(v0,dv,f,nv,vpc)
            cvreal(j,k)=real(cvi)
            cvimag(j,k)=imag(cvi)
         enddo
!--------------------
! Also get chi for negative imaginary part.
         do k=-1,nmm,-1
!            vpi(k)=vpimax*k/abs(nmm)  ! superceded with
            vpi(k)=vpimax*k/abs(nim)
! Complex v_p.
            vpc=cmplx(vp(j),vpi(k))
! Integral (which is chi.(k\lambda_{De})^2 plus the electron terms.
            cvi=-cvint(v0,dv,f,nv,vpc)
            cvreal(j,k)=real(cvi)
            cvimag(j,k)=imag(cvi)
         enddo
      enddo
!---------------------------------------------------------------
! Analytically continue fp into the lower half plane to give pole 
! contribution proportional to df/dv.
! Copy f(v) data. 
!      write(*,*)'Analytic Continuation'
      do j=1,npts
         ftfp(j)=fv(j)
      enddo
!      write(*,*)'Calling cfft1i',npts,lensav
      call cfft1i(npts,wsave,lensav,ier)
      if(ier.ne.0)stop 'cfft1i error'
!      write(*,*)'Calling cfft1f',nv+3,lensav,ier
      call cfft1f(npts,1,ftfp,npts,wsave,lensav,fftwork,nfftwork,ier)
      if(ier.ne.0)stop 'cfft1f error'      
! Now ftfp contains the fourier transform of f(). v<->t p<->omega
      dt=(nv*dv)/npts    ! Spacing of the (npts) v-values: delta v in tex.
! 1/dt is the maximum p 
      dodt=2.*3.1415926/npts   ! delta omega * delta v ->  dp*dt
      do i=0,nmm,-1        ! For each imaginary vp 
! Bandwidth limit
         do j=1,npts       ! Along the real direction
            kmod=mod(j-1+npts/2,npts)-npts/2    ! Reordered real index
            p=kmod*dodt    ! p*dt
!  Differentiate (once only)
            if(i.eq.0)ftfp(j)=ftfp(j)*cmplx(0.,1./dt)*sin(p)
            pv=p*vpi(i)/dt
! smfac is approximately 1. The max pzi should be about 5 (e-foldings).
! The max of the following expression is  efold2/(2*smfac) which is
! ~5 when efold2=10.
            pzi=pv/(1.+(pv*smfac/efold2)**2)
            ftbk(j)=ftfp(j)     *exp(pzi)
         enddo
! For testing only:
!            ftbk(j)=ftfp(j)
         call cfft1b(npts,1,ftbk,npts,wsave,lensav,fftwork,nfftwork,ier)
! then ftbk contains the analytic continuation of f' to the
! corresponding value of vpi, (as a fn of vpr) which we use as the pole
! contribution to add (times 2\pi i) to the negative vpi cases.
! Check of (last) transform. 
!       write(*,'(5f10.4)')(fv(j),ftbk(j),ftfp(j),j=1,npts)
         do j=1,npts
            if(i.eq.0)then 
            else
! Add appropriately weighted to cv
               cvimag(j,i)=cvimag(j,i)-2.*3.1415926*real(ftbk(j))
               cvreal(j,i)=cvreal(j,i)-2.*3.1415926*imag(ftbk(j))
            endif
         enddo
      enddo
!      write(*,*)'Analytic Continuation Completed'
      end
!**********************************************************************
!*********************************************************************
!      real function pvint(v0,dv,f,nv,vp)
! Integrate a distribution function f to obtain the susceptibility
! principal value integral defined as I(v_p)= \int {df(v)\over
! dv}/(v-v_p) dv. Here v_p is the phase velocity omega/k. f is the 1-d
! distribution function along v, and the resultant must be multiplied by
! (\omega_{pi}/k)^2/n_i to yield susceptibility.
! 
! The principal value integral is implemented by integrating by parts
! twice to obtain the form I = \int (v-v_p)(ln|v-v_p|-1)
! f'''(v)dv. (Which requires that f', f'' are zero at the extrema of the
! integration.  Then the value of f''' is evaluated by natural finite
! differences on a uniform mesh:
!
!  f'''_{i+1/2} = ([f_{i+2}+f_i-2f_{i+1}]- [f_{i+1}+f_{i-1}-2f_i])/dv^3
! 
! This is then integrated with I = \sum dv x value at i+1/2.  The
! limitations arise from the high order difference and rounding.  But
! the velocity function is non-singular once integrated, so no big
! difficulty arises in evaluating the integral through the former pole.
! f(v) ought to be zero for at least two values at the ends of its
! range.  The integral is done only over the non-zero values. The
! velocity grid must be uniform, and only v_0,dv is passed. 2:nv is the
! non-zero range. Convergence is second order in dv, i.e. in 1/nv.
!
! If velocities are in units of sqrt(2T_i/m_i), and the total i-th species
! density is given by n_i = \int f dv, then the i-susceptibility is
! \chi_i(\omega,k)= [Z_i^2 T_e/T_i] [-I(v_p)/2]/(n_e k^2\lambda_{De}^2)
!
! The imaginary part of the susceptibility in these units comes from the
! above equation with I(v_p) replaced by \pi df/dv giving half the residue.
! 
! Alternatively, if velocities are in units of sqrt(T_e/m_i), then
! \chi_i(\omega,k)= [Z_i^2] [-I(v_p)]/(n_e k^2\lambda_{De}^2)
!
! In either case, if f_i is normalized, one must multiply by n_i.
!
! Return pvint(=I) for the distribution function f, which is an array of
! values on a velocity grid starting at v0 (at index 0), spaced by dv,
! and extending to index nv+2. The first two and last two values would
! normally be zero, so as to satisfy boundary conditions.
! 
      real function pvint(v0,dv,f,nv,vp)
      integer nv
      real v0,vp,dv
      real f(0:nv+2)
      parameter (small=1.e-20)
      pvint=0.
      do i=1,nv
         v=v0+(i+0.5)*dv
         vd=v-vp
         if(abs(vd).lt.small)then 
            vfn=0.
         else
            vfn=vd*(alog(abs(vd))-1)
         endif
         pvint=pvint+vfn*((f(i+2)-f(i-1))-3.*(f(i+1)-f(i)))
      enddo
      pvint=pvint/dv**2
      end
!*********************************************************************
!      complex function cvint(v0,dv,f,nv,vp)
! Perform the real and imaginary susceptibility integrations for a complex
! \omega/k =v_p with non-zero imaginary part, along the real v-axis.
! Denoting the real and imaginary parts of vp by vr, and vi respectively,
! the integral is
! J = \int (df/dv) (v-vr)/[(v-vr)^2+vi^2] df
!         + i \int (df/dv) vi /[(v-vr)^2 + vi^2] df
!   = -\int (d^2f/dv^2)[ln[{(v-vr)/vi}^2+1]/2 + i arctan{(v-vr)/vi}] dv.
!
! Interpolation to f''(i) = (f(i-1)-2*f(i)+f(i+1))/dv^2 
! So f''(i+1/2)=(f(i-1)-f(i)-f(i+1)+f(i+2))/2.dv^2
!
! f must go to zero at the ends of the integration to make this form valid.

      complex function cvint(v0,dv,f,nv,vp)
      complex vp
      integer nv
      real v0,dv
      real f(0:nv+2)
      parameter (small=1.e-6)
      real cr,ci

      vr=real(vp)
      vi=imag(vp)
      if(abs(vi).lt.small)then
         write(*,*)'cvint error: vi too small',vi
         cvint=cmplx(0.,0.)
         return
      endif
      cr=0.
      ci=0.
      do i=1,nv
!         v=v0+(i+0.5)*dv
! Four-point formula f''= this combination /2dv^2
!         fpp=(f(i-1)-f(i)-f(i+1)+f(i+2))

         v=v0+(i)*dv
!c Three-point formula. f''= this combination/dv^2
         fpp=(f(i-1)-2.*f(i)+f(i+1))
         va=(v-vr)/vi
         cr=cr-alog(va**2+1.)*fpp
         ci=ci-atan(va)*fpp
      enddo
! Four-point:
!      cr=cr*0.25/dv
!      ci=ci*0.5/dv
! Three-point:
      cr=cr*0.5/dv
      ci=ci*1./dv
      
      cvint=cmplx(cr,ci)
      end

!********************************************************************
! Return the 1-D distribution function in array f(0:nv+2) for nv+3
! velocity points from v0 to v0+(nv+2)*dv being the sum of ng gaussians
! whose parameters are in gpar(*,ig)=(den, vel, vtperp, vtparallel) 
! At an angle of velocity component to drift velocity, theta.

! gpar(3,4) here have already been scaled dividing by Te and velocities
! dv, v0, gpar(2) by sqrt(T0).  Consequently the maxwellian is
! exp(-v^2/(2*vt2)), vt2=Ti/T0 and the velocities v0 dv etc can be
! considered to be in units of sqrt(T0/mi).
      subroutine fvinit(nv,dv,v0,f,ng,gpar,theta)
      integer nv
      real dv,v0
      real f(0:nv+2)
      integer npar
      parameter (npar=4)
! Parameters; den, vel, vtperp, vtparallel.      
      real gpar(npar,ng)
      real theta

      do i=0,nv+2
         f(i)=0.
      enddo
      ct=cos(theta)
      st=sin(theta)
      do k=1,ng
         vt2=gpar(3,k)*st**2+gpar(4,k)*ct**2
         fnorm=1/(sqrt(2.*3.141593*vt2))
         do i=0,nv+2
            v=v0+i*dv
            f(i)=f(i)+gpar(1,k)*fnorm*exp(-(v-gpar(2,k)*ct)**2           &
     &           /(2.*vt2))
         enddo
      enddo
      end
!*********************************************************************
      subroutine realvpplot(vp,cvreal,cvimag,npts,nimax)
      integer npts,nimax
      real vp(npts),cvreal(npts,0:nimax),cvimag(npts,0:nimax)
      call autoplot(vp,cvreal(1,0),npts)
      call polymark(vp,cvreal(1,0),npts,1)
      call vecw(vp(1),0.,0)
      call vecw(vp(npts),0.,1)
!      call axlabels('v!dp!d/(2Ti/m_i)!u1/2!u',
!     $     '!Ax!@ .T!di!d(k!Al!@!dD!d)!u2!u/T!de!d')
      call axlabels('v!dp!d/(T!de!d/m!di!d)!u1/2!u',                     &
     &     '!Ax!@ .(k!Al!@!dD!d)!u2!u')
      call legendline(.7,.9,0,'!AR!@(!Ax!@!di!d)')
      call boxtitle                                                      &
     &     ('Susceptibility along the real !Aw!@/k=v!dp!d axis.')
      call winset(.true.)
      call dashset(2)
      call color(2)
      call polyline(vp,cvimag(1,0),npts)
      call legendline(.7,.85,0,'!AI!@(!Ax!@!di!d)')
! cv lines:
      call color(6)
      call dashset(3)
      call polyline(vp,cvreal(1,1),npts)
      call dashset(4)
      call polyline(vp,cvimag(1,1),npts)
      call color(15)
      call dashset(0)
      call pltend()
      end
!******************************************************************
      subroutine argdiagplot(pv,fvp,npts,Te)
      integer npts
      real pv(npts),fvp(npts),Te
! Argand diagram plot.
      call autoplot(pv,fvp,npts)
      call axlabels('Real','Imaginary')
      call winset(.true.)
! Set the space between labels to .05 normal units and line drawing off.
      call labeline([.05],[0.],0,'',-99)
      call labeline(pv,fvp,npts,'>',1) 
      call labeline([.005],[0.],0,'',-99)
      call labeline(pv,fvp,npts,'"',1)
      call labeline([0.],[0.],0,'',-99)
! Draw axes through origin
      call vecw(-Te,0.,0)
      call vecw(Te,0.,1)
      call vecw(0.,-Te,0)
      call vecw(0.,Te,1)
      call pltend()
      end
!*******************************************************************
! Contour the real and imaginary parts of chi on 2D grid of vp
! vp(1:npts) vpi(nmm:nim) cworka is a work array at least as big.
! laspect if true says plot conserving aspect ratio.  
      subroutine chicontour(vp,vpi,npts,nrmin,nrmax,nim,nmm,cvreal       &
     &     ,cvimag,ceimag,cworka,laspect,lcolor,ltwotone)
      integer npts,nim,nmm,nrmin,nrmax
      real vp(npts),vpi(nmm:nim),cvreal(npts,nmm:nim),cvimag(npts,       &
     &     nmm:nim),ceimag(npts,nmm:nim)
      character cworka(npts,nmm:nim),ilb
      logical laspect,lcolor
      logical ltwotone
      common/ilabel/ilb

      include 'acpathcom.f'

      parameter (npixel=240)
      integer ired(npixel),igreen(npixel),iblue(npixel)
      parameter (ncl=9)
      real zcl(-ncl:ncl)
! Cunning contours so as to be able to do coloring
      data zcl/-120.,-50.,-20.,-10.,-5.,-2.,-1.,-.5                      &
     &     ,-.2,0.,.2,.5,1.,2.,5.,10.,20.,50.,120./
      if(laspect)then
         call pltinaspect(vp(nrmin),vp(nrmax),vpi(nmm),vpi(nim))
      else
         call pltinit(vp(nrmin),vp(nrmax),vpi(nmm),vpi(nim))
      endif

      if(lcolor)then
         do k=1,npixel 
            if(k.le.npixel/2-1.or.(k.le.npixel/2.and..not.ltwotone))then
               ired(k)=65535
               iblue(k)=65535
               igreen(k)=65535
            elseif(k.le.npixel/2)then
               ired(k)=58000
               iblue(k)=58000
               igreen(k)=58000
            else
               ired(k)=50000
               iblue(k)=50000
               igreen(k)=50000
            endif
         enddo
         call accisgradset(ired,igreen,iblue,npixel)
         icl=2*ncl+1
! Vector, Color, Omit lines.
         icsw=1+16+32+64
         call contourl(cvreal(nrmin,nmm),cworka,npts,nrmax-nrmin+1,nim   &
     &        -nmm+1,zcl(-ncl),icl,vp(nrmin),vpi,icsw)
      endif

! Draw zero line
      call color(12)
      call vecw(vp(nrmin),0.,0)
      call vecw(vp(nrmax),0.,1)
      call color(3)
! Vector, negative real part contours.
      icsw=1
      icl=ncl+1
      call contourl(cvreal(nrmin,nmm),cworka,npts,nrmax-nrmin+1,nim-nmm  &
     &     +1,zcl(-ncl),icl,vp(nrmin),vpi,icsw)
      call color(15)
      call axis()
      call axis2()
! Correct the cvreal
!      if(ltwotone)cvreal=cvreal+1.
! Zero imaginary part contour,
      icsw=1
      icl=1
      if(.not.ltwotone)iacpsw=1    ! Path documentation
      call contourl(cvimag(nrmin,nmm),cworka,npts,nrmax-nrmin+1,nim-nmm  &
     &     +1,zcl(0),icl,vp(nrmin),vpi,icsw)
      if(ltwotone)then
         call axlabels('!AR!@(v!dp!d) /(T!d0!d/m!di!d)!u1/2!u',          &
     &        '  !AI!@(v!dp!d) /(T!d0!d/m!di!d)!u1/2!u')
         call legendline(.6,.9,258,                                      &
     &        '(k!Al!@!dD!d)!u2!u!Ax!@!di!d contours')
         do k=nmm,nim
            do j=1,npts
               cvimag(j,k)=cvimag(j,k)+ceimag(j,k)
            enddo
         enddo
! Second imaginary contour
         call color(1)
         icl=-1
! Turn on path documentation.
         iacpsw=1
         iacpcp=1
         iacpcon(iacpcp)=0
         iacppt=0
         call contourl(cvimag(nrmin,nmm),cworka,npts,nrmax-nrmin+1,nim   &
     &        -nmm+1,zcl(0),icl,vp(nrmin),vpi,icsw)
         call legendline(.77,.65,0,'imag(!Ax!@)')
         iacpsw=0
         call color(15)
      else
         call axlabels('!AR!@(v!dp!d) /(T!de!d/m!d'//ilb//'!d)!u1/2!u',  &
     &        '  !AI!@(v!dp!d) /(T!de!d/m!d'//ilb//'!d)!u1/2!u')
         call legendline(.6,.9,258,'(k!al!@!dDe!d)!u2!u!Ax!@ contours')
      endif
      call legendline(.77,.75,0,'imag')
      call color(3)
      call legendline(.77,.85,0,'real')
      call color(15)
      call pltend()
      end
!*******************************************************************
!*******************************************************************
      subroutine plotfofv(vp,fv,npts,fmax)
      integer npts
      real vp(npts),fv(npts),fmax
      yt=1.15*fmax
      yb=-0.06*yt
!      yb=-0.18*yt
      call pltinit(vp(1),vp(npts),yb,yt)
      call axptset(0.,abs(yb)/(yt-yb))
      call axis()
      call axis2()
      call polyline(vp,fv,npts)
      call axlabels('','f(v)')
      end
!********************************************************************
! Calculate and plot the result of the path/contour following.
      subroutine acpplot(vp,vpi,vavei,npts,nim,nmm,cvreal,cvimag,cworka  &
     &     ,laspect,lcolor)
      integer npts,nim,nmm
      real vp(npts),vpi(nmm:nim),cvreal(npts,nmm:nim),cvimag(npts,       &
     &     nmm:nim)
      character cworka(npts,nmm:nim)
      logical laspect,lcolor

! iacpcon is the pointers to starts of contours which are 
! xacp(iacpcon(i)+1),yacp(iacpcon(i)+1) 
! iacpcon(i)+1 is the start index of each contour segment.
! iacpcon(i+1)-iacpcon(i) is the length of contour i.
! iacpcp is the current contour number, and when completed
! is 1+ the total number of contours. iacpcon(iacpcp) points to 
! the end of the last contour.
      include 'acpathcom.f'
      real curk(iacpmax),curfr(iacpmax),curfi(iacpmax)
      character*20 label

      iplotinit=0
! Do over contour (segments)
      curfimx=-1.e6
      ncur=0
      xmax=3.2
      xmin=0.
      wimfac=1.
      do i=1,iacpcp-1
         icur=0
         id=0
         vpimx=-1.e6
         cfimin=1.e6
         cfimax=-1.e6
         cfrmin=1.e6
         cfrmax=-1.e6
         do j=1,iacpcon(i+1)-iacpcon(i)
! Convert x,y values to fractional indexes then interpolate for chi
            xf=(1.+(xacp(iacpcon(i)+j)-vp(1))/(vp(npts)-vp(1))*(npts-1))
            yf=(1.+(yacp(iacpcon(i)+j)-vpi(nmm))/(vpi(nim)-vpi(nmm))     &
     &           *(nim-nmm))
            chirept=amultilin2(cvreal(1,nmm),npts,npts,nim-nmm+1,        &
     &           xf,yf,ierr)
! Add only if in the negative chi-real region
            if(chirept.lt.0.)then
! Add point to curve
               icur=icur+1
               curk(icur)=sqrt(-chirept)
! Real part of this phase velocity
!               vprm=(vp(1) +(xacp(iacpcon(i)+j)-1)*(vp(npts)-vp(1))
!     $              /(npts-1))
!               vprm=xacp(iacpcon(i)+j) ! Old v relative to electrons
               vprm=xacp(iacpcon(i)+j)-vavei  ! relative to vavei
! Hence real part of frequency k.v_p:
               curfr(icur)=curk(icur)*vprm
! Imaginary part of phase velocity
               vpim=yacp(iacpcon(i)+j)
! Imaginary part of frequency, not now divided by scaling factor.
               curfi(icur)=curk(icur)*vpim
!       write(*,'(i5,5f10.4)')icur,curk(icur),curfi(icur),vprm,vpim
               if(curk(icur).lt.xmax)then ! Adust imag limits
                  if(curfi(icur).gt.cfimax)cfimax=curfi(icur)
                  if(curfi(icur).lt.cfimin)cfimin=curfi(icur)
                  if(curfr(icur).gt.cfrmax)cfrmax=curfr(icur)
                  if(curfr(icur).lt.cfrmin)cfrmin=curfr(icur)
!       write(*,*)cfimin,cfimax,cfrmin,cfrmax
               endif
               if(vpim.gt.vpimx)vpimx=vpim
               if(curfi(icur).gt.curfimx)curfimx=curfi(icur)
! Decide the position of the line's label:
               xpos=xmax*(.6-.05*i)
               if(curk(icur).lt.xpos .and. id.eq.0)id=icur
!               write(*,'(2i4,7f9.4)')i,j,xacp(iacpcon(i)+j)
!     $              ,yacp(iacpcon(i)+j),curk(icur),curfr(icur)
!     $              ,curfi(icur),vprm,vpim
            endif
         enddo
! Segment completed
!         write(*,*)id,icur,vpimx,cfimin,cfimax
!         if(id.gt.0.and.icur.gt.1 .and. vpimx.gt.-0.05)then
         if(id.gt.0.and.icur.gt.1)then
            if(iplotinit.eq.0)then !Initialize plot
!         write(*,'(10f8.3)')(curfr(k),k=1,icur-1)
            ymax=cfimax
            ymin=cfimin
            wimfac=(cfimax-cfimin)/(cfrmax-cfrmin)
            wimfaclog=nint(log10(wimfac))
            wimfacl5=wimfaclog-.5
            wimfac=10.**(min(wimfaclog,wimfacl5))
!            write(*,*)'wimfac',wimfac
            ymin=min(-.2,-max(abs(cfimin),abs(cfimax)))
            ymax=max(.2,max(abs(cfimin),abs(cfimax)))
            call pltinit(xmin,xmax,ymin,ymax)
            call axis()
            call axlabels('k!Al!@!dD!d','!Aw!@/!Aw!@!dpi!d')
            call axptset(1.,0.)
            call ticrev()
            call altyaxis(1./wimfac,1./wimfac)
            call axlabels('','!AI!@(!Aw!@)/!Aw!@!dpi!d')
            call ticrev()
            call axptset(0.,0.)
            call winset(.true.)
            call vecw(0.,0.,0)
            call vecw(5.,0.,1)
            call legendline(.63,.95,0,'!AR!@(!aw!@)-k!p!o_!o!qv!di!d')
            call dashset(2)
            call legendline(.02,.95,0,'!AI!@(!Aw!@)')
            call dashset(0)
            iplotinit=1
            endif
            ncur=ncur+1
! There is a curve to plot. Plot it.
            call dashset(0)
            call color(ncur)
            call iwrite(ncur,iwdth,label)
!            call labeline(curk,curfr,icur,label,10)
!            id=0.7*icur
            call polyline(curk,curfr*wimfac,icur)
            call jdrwstr(wx2nx(curk(id))+0.02,wy2ny(curfr(id)*wimfac)    &
     &           ,label,0.)
            call dashset(2)
            call polyline(curk,curfi,icur)
            call jdrwstr(wx2nx(curk(id)),wy2ny(curfi(id)),label,0.)
!            write(*,*)'   k       Re(o)       Im(o)  Seg',i
!     $           ,j
!            do k=1,icur
!               write(*,'(3f10.4)')curk(k),curfr(k),curfi(k)
!            enddo
         endif
      enddo
      if(iplotinit.gt.0)call pltend()
      end
!******************************************************************
! Use interp at two columns to interpolate values; then interpolate
! those two values. Hence giving a multilinear interpolation in 2D.
      real function amultilin2(data,ld,nx,ny,xf,yf,ierr)
      real data(ld,ny)
      real xf,yf
      ierr=0
      ix=int(xf)
      iy=int(yf)
      yff=yf-iy
      xff=xf-ix
! Allow values to the upper limit
      if(ix.eq.nx)then
         ix=ix-1
         xff=1.
      endif
      if(iy.eq.ny)then
         iy=iy-1
         yff=1.
      endif
! But not beyond.
      if(ix.lt.1.or.ix.ge.nx.or.iy.lt.1.or.iy.ge.ny)then
         ierr=1
         write(*,*)'amultlin2 index error',xf,nx,yf,ny
      endif
      f11=data(ix,iy)
      f12=data(ix,iy+1)
      f21=data(ix+1,iy)
      f22=data(ix+1,iy+1)
      d1=f11*(1-xff)+f21*xff
      d2=f12*(1-xff)+f22*xff
      amultilin2=d1*(1.-yff)+d2*yff
      end
!********************************************************************
!******************************************************************
      subroutine parsecmdline(gpar,npar,ngmax,ng,vdmin,vdmax,vw,Te,theta &
     &     ,vwfac,vpimax,nim,nmm,nimax,smfac,laspect,lcolor,eldweight    &
     &     ,vflat,wflat,fflat,lextra,amp,nw,lgrowth,lcoptic,cfilename    &
     &     ,ltwotone,vrgfac,lthresh,rmitome,omegac,lcombplot,gepar,nge)
      integer npar,ngmax
      real gpar(npar,ngmax),gepar(npar,ngmax)
      real vdmin,vdmax,vw,Te,theta,vwfac,amp
      logical laspect,lcolor,lextra,lgrowth,lcoptic,ltwotone,lthresh
      logical lcombplot
      character*50 cfilename

      character*50 argument
      logical lgset,lgeset
      data lgset/.false./lgeset/.false./
      argument=' '
      eld=1.
      nw=0
      ng=1
      nge=1

! Deal with arguments.
      do ia=1,iargc()
         call getarg(ia,argument)
         write(*,*)argument(1:20)
         if(argument(1:1).ne.'-')then
            if(.not.lcoptic)then
!        Read the input file.
               write(*,*)'Attempting to open 12',iargc()
               open(13,file=argument,form='formatted',status='old',err   &
     &              =101)
               read(13,*,end=102)theta,Te
               write(*,'(a,f8.3,a,f8.3)')'theta=',theta,'  Te=',Te
               ct=cos(theta)
               if(ct.le.0)then
                  write(*,*)'CosTheta is negative. Not allowed.',theta
                  stop
               endif
               read(13,*,end=102)ng
               if(ng.gt.ngmax)then
                  write(*,*)'Too many gaussians'
                  stop
               endif
               vdmin=0.
               vdmax=0.
               write(*,'(a)')'          Den   Velocity   Tperp   Tpara'
               do i=1,ng
                  read(13,*,end=102)(gpar(k,i),k=1,4)
                  write(*,'(a,i1,a,4f8.3)')'gpar(',i,')=',(gpar(k,i),k=1 &
     &                 ,4)
                  if(gpar(2,i).lt.vdmin)vdmin=ct*gpar(2,i)
                  if(gpar(2,i).gt.vdmax)vdmax=ct*gpar(2,i)
                  vwhere=vwfac*sqrt(sqrt(1-ct**2)*gpar(3,i)              &
     &                 +ct**2*gpar(4,i))
                  if(vwhere.gt.vw)vw=vwhere
               enddo
            else
!          File name of coptic data.
               cfilename=argument
            endif
         else
! Switch arguments.
            if(argument(1:2).eq.'-i')read(argument(3:),*,err=103)vpimax
            if(argument(1:2).eq.'-n')read(argument(3:),*,err=103)nim
            if(argument(1:2).eq.'-m')read(argument(3:),*,err=103)nmm
            if(argument(1:2).eq.'-s')read(argument(3:),*,err=103)smfac
            if(argument(1:2).eq.'-w')read(argument(3:),*,err=103)vwfac
            if(argument(1:2).eq.'-u')read(argument(3:),*,err=103)amp
            if(argument(1:2).eq.'-x')read(argument(3:),*,err=103)vrgfac
            if(argument(1:2).eq.'-p')read(argument(3:),*,err=103)nw
            if(argument(1:2).eq.'-T')read(argument(3:),*,err=103)Te
            if(argument(1:2).eq.'-O')read(argument(3:),*,err=103)omegac
            if(argument(1:6).eq.'-theta')then 
               read(argument(7:),*,err =103)theta  !degrees
               theta=3.1415926/180.*theta          !radians
            elseif(argument(1:2).eq.'-t')then
               ltwotone=.not.ltwotone
            endif
            if(argument(1:2).eq.'-r')read(argument(3:),*,err=103)rmitome
            if(argument(1:3).eq.'-ve')then
               if(lgeset)nge=nge+1
               read(argument(4:),*,err=103,end=105)(gepar(k,nge),k=1,4)
               lgeset=.true.
            elseif(argument(1:2).eq.'-v')then
               if(lgset)ng=ng+1
               read(argument(3:),*,err=103,end=105)(gpar(k,ng),k=1,4)
               lgset=.true.
            endif
            if(argument(1:3).eq.'-Vi')read(argument(4:),*,err=103)vdmin
            if(argument(1:3).eq.'-Va')read(argument(4:),*,err=103)vdmax
            if(argument(1:2).eq.'-h')goto 104
            if(argument(1:2).eq.'-?')goto 104
            if(argument(1:2).eq.'-a')laspect=.not.laspect
            if(argument(1:2).eq.'-c')lcolor=.not.lcolor
            if(argument(1:2).eq.'-b')lextra=.not.lextra
            if(argument(1:2).eq.'-g')lgrowth=.not.lgrowth
            if(argument(1:2).eq.'-S')lthresh=.not.lthresh
            if(argument(1:2).eq.'-e')read(argument(3:),*,err=103)eld
            if(argument(1:2).eq.'-d')call pfset(-3)
            if(argument(1:7).eq.'-COPTIC')lcoptic=.not.lcoptic
            if(argument(1:3).eq.'-CP')lcombplot=.not.lcombplot
            if(argument(1:2).eq.'-f')then
               if(fflat.eq.0.)then
                  fflat=1.
                  vflat=.5
                  wflat=.1
               endif
               if(argument(1:3).eq.'-fv')read(argument(4:),*,err         &
     &              =103)vflat
               if(argument(1:3).eq.'-fw')read(argument(4:),*,err         &
     &              =103)wflat
               if(argument(1:3).eq.'-ff')read(argument(4:),*,err         &
     &              =103)fflat
!               write(*,*)fflat,vflat,wflat
            endif
         endif
 105     continue
      enddo

      eldweight=eldweight*eld
      if(iargc().le.0)then
         write(*,*)'No input file or switches specified. Using defaults'
      endif
      if(nim.gt.nimax)nim=nimax
      if(abs(nmm).gt.nimax)nmm=nimax
      if(nmm.gt.0)nmm=-nmm

      return

 101  write(*,*)'Failed to open file: ',argument
 103  write(*,*)'Error parsing argument: ',argument
 104  write(*,*)'Usage: chiofv [switches] [filename]'
      write(*,'(a,f7.4,a,f7.4,a)')' -theta k-angle to v,B (deg)[',theta, &
     &     ']rad;  -O omegac_e/omega_pe [',omegac,']'
      write(*,'(a,f7.3,a)')' -T Te                      [',Te,']'
      write(*,'(a,f7.1,a)')' -r mi/me                   [',rmitome,']'
      write(*,'(a,f7.4,a)')' -i max imaginary vp        [',vpimax,']'
      write(*,'(2a,2i4,a)')' -n -m imag mesh numbers  '                  &
     &     ,'  [',nim,nmm,'] +ve,-ve. Negative range=n/m*vpimax'
      write(*,'(a,f7.4,a)')' -u uncertainty (noise)     [',amp,']'
      write(*,'(a,f7.2,a,f7.2,a)')' -Va maximum vreal of plot  [',       &
     &     vdmax,'] -Vi minimum vreal of plot  [',vdmin,']'
      write(*,'(a,f7.4,a)')' -x vrange fraction plot    [',vrgfac,']'
      write(*,'(a,f7.3,a,f7.3,a)')' -s smoothing factor        [',       &
     &     smfac,'] -w integration width factor[',vwfac,']'
      write(*,*)'-a toggle aspect ratio      [',laspect                  &
     &     ,'] -c toggle coloring          [',lcolor,']'
      write(*,*)'-b toggle extra detail plots[',laspect,']',             &
     &     ' -g toggle growth rate plot  [',lgrowth,']'
      write(*,*)'-t toggle two tone plot     [',ltwotone,']',            &
     &     ' -S stability threshold      [',lthresh,']'
      write(*,*)'-CP toggle combined f plot  [',ltwotone,']'
      write(*,*)'-d disable display of plots'
      write(*,*)'-e ELD weight factor'
      write(*,*)'-p integer width of smoothing triangle box',nw
      write(*,'(a,4f7.2,a)')' -ve<n>,<v>,<Tpp>,<Tpl> Electron Gaussian ['  &
     &     ,(gpar(i,ng),i=1,4),']'
      write(*,'(a,4f7.2,a)')' -v<n>,<v>,<Tpp>,<Tpl>  Ion Gaussian      ['  &
     &     ,(gpar(i,ng),i=1,4),']'
      write(*,*)'   First -v[e] overwrites default, succeeding adds new'
      write(*,*)'-f flatspot: -fv center-velocity -fw width'             &
     &     ,' -ff slope-factor'
      write(*,*)'-COPTIC : toggle external distribution reading lcoptic'
      write(*,*)'   subsequent filename refers to external data.'
      call exit
 102  write(*,*)'File ended prematurely',ng,k    
      end

!*******************************************************************
! Create a flatspot by adding a zero-integral function at 
! position iv, of width iw
      subroutine flatspot(nv,dv,f,iv,iw,fflat)
      real f(0:nv+2)
      integer nv,iv,iw

      if(iw.le.1)return
      if((iv+2*iw).gt.(nv+2).or.((iv-2*iw).lt.0))then
         write(*,*)'Flat spot overlaps ends',iw,iv,nv
         return
      endif

      slope=(f(iv+1)-f(iv-1))/(2.*dv)

! Triangle shaped add-on has analyticity problems.
      if(.false.)then
      do i=iv-2*iw,iv+2*iw
         ir=i-iv
         if(abs(ir).gt.iw)ir=sign(2*iw-abs(ir),ir)
!         write(*,*)i,ir
         f(i)=f(i)-fflat*slope*dv*ir
      enddo
      endif

! Small gaussian derivative.
      do i=0,nv+2
         ir=i-iv
         f(i)=f(i)-fflat*slope*dv*ir*exp(-(float(ir)/iw)**2)
      enddo

      end
!******************************************************************
! This replaces a dummy routine in the accis library. 
! Stores the just completed contour of length imax and value cv.
!      subroutine acpathdoc(z,cv,L,imax,xc,yc,i)
      subroutine acpathdoc(imax,xc,yc,cv)
      integer imax
      real cv,xc(imax),yc(imax)
! common block for communications includes iacp...
      include 'acpathcom.f'
      if(iacpsw.ne.0)then
! Switched on:
         if(imax.ne.0.)then
! Store onto the contour stack.
            do i=1,imax
               if(iacppt.lt.iacpmax)then
                  iacppt=iacppt+1
                  xacp(iacppt)=xc(i)
                  yacp(iacppt)=yc(i)
!                  if(i.eq.1)write(*,'(a,3i4,3f9.3)')'acpathdoc',
!     $        iacppt,iacpcp,iacpcon(iacpcp),cv,xacp(iacppt),yacp(iacppt)
               else
                  write(*,*)'Exhausted iacpmax length',iacpmax
                  return
               endif
            enddo
! Terminate a contour and store its end.
            if(iacpcp.lt.iacpconmax)then
               iacpcp=iacpcp+1
               iacpcon(iacpcp)=iacppt
            else
               write(*,*)'Exhausted iacpconmax length',iacpconmax
               return
            endif
         endif
!         write(*,'(a,3i4,3f9.3)')'acpathdoc',
!     $        iacppt,iacpcp,iacpcon(iacpcp),cv,xacp(iacppt),yacp(iacppt)
      else ! Reset call
         iacpcp=1
         iacppt=0
      endif
! iacpcon is the pointer to start of contour which is 
! xacp(iacpcon(i)+1),yacp(iacpcon(i)+1) 
! iacpcon(i)+1 is the start index of each contour segment.
! iacpcon(i+1)-iacpcon(id) is the length of contour i.
! iacpcp is the current contour number, and when completed
! is 1+ the total number of contours. iacpcon(iacpcp) points to 
! the end of the last.
      end
!**********************************************************************
      subroutine addnoise(nv,f,amp)
      integer nv
      real f(0:nv+2)

! Initialize random number generator.
      df=ran1(-1)
      do i=3,nv-1
         df=ran1(0)
!         write(*,*)i,' ran1=',df
         f(i)=f(i)+df*amp
      enddo

      end
!**********************************************************************
      FUNCTION ran1(IDUM)
! Returns a uniform random deviate between 0 and 1.
! Reentrant version of ran1 makes the state visible in common.
! Get explicit 
      implicit none
      real ran1
      integer idum
! Common here requires separation of processes for reentrancy.
! So this is not thread safe.
!      include 'ran1com.f'
! Internal state of the fortran random number generator ran1 is made
! visible through this common block. It can then be saved by just
! writing ranstate(1:103) and restored by reading it back.
! This requires integers to be no longer than reals.
      real ranstate(103)
! Internal state of ran1:
      real rrnd(97)
      integer irx1,irx2,irx3,jrc
      equivalence (rrnd,ranstate)
      equivalence (irx1,ranstate(98)),(irx2,ranstate(99)),               &
     &     (irx3,ranstate(100)),(jrc,ranstate(101))
! Internal state of the gaussian random generator gasdev:
      integer gd_iset
      real gd_gset
      equivalence (gd_iset,ranstate(102)),(gd_gset,ranstate(103))
! The whole thing:
      common /ran1com/ranstate

      real RM1,RM2
      integer M1,M2,M3,ic1,ic2,ic3,ia1,ia2,ia3,j
      save
!      DIMENSION Rrnd(97)
      PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
      PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
      PARAMETER (M3=243000,IA3=4561,IC3=51349)
! Can't do auto initialize. Have to do it by hand.
!      DATA IFF /0/
!      IF (IDUM.LT.0.OR.IFF.EQ.0) THEN
      IF (IDUM.LT.0) THEN
!        IFF=1
        IrX1=MOD(IC1-IDUM,M1)
        IrX1=MOD(IA1*IrX1+IC1,M1)
        IrX2=MOD(IrX1,M2)
        IrX1=MOD(IA1*IrX1+IC1,M1)
        IrX3=MOD(IrX1,M3)
        do J=1,97
           IrX1=MOD(IA1*IrX1+IC1,M1)
           IrX2=MOD(IA2*IrX2+IC2,M2)
           Rrnd(J)=(FLOAT(IrX1)+FLOAT(IrX2)*RM2)*RM1
        enddo
! Tell gasdev to reset.
        gd_iset=0
      ENDIF
      IrX1=MOD(IA1*IrX1+IC1,M1)
      IrX2=MOD(IA2*IrX2+IC2,M2)
      IrX3=MOD(IA3*IrX3+IC3,M3)
      JRC=1+(97*IrX3)/M3
      if(JRC.GT.97.OR.JRC.LT.1)then
         write(*,*)'RAN1 Error!'
      endif
      RAN1=Rrnd(JRC)
      Rrnd(JRC)=(FLOAT(IrX1)+FLOAT(IrX2)*RM2)*RM1
      RETURN
      END
!*****************************************************************
      subroutine trismooth(n,nw,v,ww,vs)
! Smooth a vector using a triangular box.
! So vs= 1/(nw+1) sum_{-nw}^{+nw} v_j*(1-j/(nw+1))
! On entry
!   n is the length of v
!   nw is the smoothing length in samples (0 => do nothing)
!   v is the input vector
!   ww is a work space of length n
! On exit
!   vs is the smoothed result. 
! The arrays v and vs or the arrays v and ww, but not ww and vs, 
! can be the same storage locations.
! If the vector is on a non-uniform grid, then this smoothing is not
! area conserving, because it knows nothing about the grid spacing. 
      integer n,nw
      real v(n),vs(n),ww(n)
      do i=1,n
         ww(i)=v(i)
      enddo
      p=nw+1.
      do i=1,n
         vi=0.
         do j=-nw,nw
! Deal with the ends...
            k=min(n,max(1,i+j))
            vi=vi+ww(k)*(1.-abs(float(j))/p)
         enddo
         vs(i)=vi/p
      enddo
      end
!*******************************************************************
! Test the trismooth routine on a step function.
      subroutine trismoothtest(nw)
      integer nw,nt
      parameter (nt=200)
      real v(nt),vs(nt),ww(nt),x(nt)

      do i=1,nt
         x(i)=i
         if(i.lt.nt/2)then
            v(i)=1.
         else
            v(i)=0.
         endif
      enddo
      call trismooth(nt,nw,v,ww,vs)
!      write(*,*)v,vs
      call autoplot(x,v,nt)
      call color(1)
      call polyline(x,vs,nt)
      call pltend()
      end
!**********************************************************************
! Test the COPTIC input data reading
      subroutine COPTICverif(nv,dv,v0,f)
      integer nv
      real dv,v0,f(0:nv+2)
      open(14,file='COPTICverif.dat',status='unknown',err=101)
      close(14,status='delete')
      open(14,file='COPTICverif.dat',status='new',err=101)
          write(14,*) nv
          write(14,*) dv
          write(14,*) v0
          write(14,'(1g14.6)')(f(i),i=0,nv+2)
      close(14)
      goto 102
 101  write(*,*)'Error opening file'
      close(14,status='delete')
      call exit(1)
 102  continue
      end
!*********************************************************************
      subroutine stationary(nv,f,ngmax,fstat,nstat)
      real f(nv),fstat(ngmax)
! Finds stationary values of f and return nstat values in fstat
      nstat=0
      do i=2,nv-1
      if((f(i+1)-f(i))*(f(i)-f(i-1)).lt.0..or.f(i+1)-f(i).eq.0)then
            nstat=nstat+1
            fstat(nstat)=f(i)
         endif
         if(nstat.eq.ngmax)return
      enddo
      end
