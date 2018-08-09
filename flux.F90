      PROGRAM flux
!===============================================================================
!
!     PURPOSE:
!      to process FLUX data (obtained from raw data)
!         1. Coordinate rotation
!         2. WPL correction
!         3. Calculation of Potential evaporation
!         4. Calculation of Priestely-Taylor coefficient
!         5. Calculation of conductances
!         6. Calculation of Decoupling factor
!
!     Record of revisions:
!     Date        Programmer         Description of change
!     ====        ==========         =====================
!  Feb/28/2002    Jinkyu Hong            Original code
!  Mar/09/2002    Jinkyu Hong           Add WPL correction
!  Mar/14/2002    Jinkyu Hong           Add quality control
!  Mar/21/2002    Jinkyu Hong     Add the calculation of conductances
!  Aug/17/2005    Jinkyu Hong     Modification of WPL and PF rotation routines
!  Dec/23/2008    Jinkyu Hong     Update for new data/modify for pgf compiler
!  Jan/07/2009    Jinkyu Hong     Update for handling bad data points
!  Jan/13/2009    Jinkyu Hong     Modification for a new data logger (CR23X)
!
!     Contact information: jhong@yonsei.kr / hong.jinkyu@gmail.com
!
!===============================================================================
!
!                                CONDITIONS OF USE
!
!   This program is available for scientific research. The FORTRAN code and 
! documentation are available subject to the following conditions. I assume 
! that you accept the following conditions of use if you use the program.
!
!   1. It shall not be used for any commercial applications.
!
!   2. This program has been developed for the purpose of Biometeorology Lab, 
! Yonsei University and no warranty is given as to its suitability for other 
! users.
!
!   3. Jinkyu Hong can not guarantee that this program will work correctly in all 
! circumstances. Neither does Jinkyu Hong and Biometeorology lab, Yonsei University
! accept any liability whatsoever for any error or omission in the program.
!
!   4. Jinkyu Hong retains the proprietary rights and copyright of the program.
!
!   5. Suitable acknowledgement of Jinkyu Hong and Biometeorology Lab, Yonsei 
! University shall be made in any published work involving this program.
!
!   6. You can share this program with all persons and groups without permission 
! of Biometeorology Lab, Yonsei University, only when you accept the above 
! CONDITIONS OF USE.
!===============================================================================
!  
!     Description of Variables
!      nt = total number of lines to be processed
!      nc = total number of variables (column) to be processed
!
!      um(3) = mean wind speed
!      su(4) = standard deviation of wind speed and temperature
!      xcov(4,3) = covariances
!      ta = Air temperature measured by slow response sensor like HMP45C
!      ea = Water vapor pressure measured by slow response sensor (hPa)
!      co2 = CO2 concentration (mg/m3)
!      rn = net radiation (W/m2)
!      rflx(1) = downward shortwave radiation (W/m2)
!      rflx(2) = upward shortwave radiation (W/m2)
!      rflx(3) = downward longwave radiation (W/m2)
!      rflx(4) = upward longwave radiation (W/m2) 
!      g = soil heat flux (W/m2)
!      prec = precipitation rate (mm)
!      pa = atmospheric pressure (kPa) 
!      vdf = vapor pressure deficit (Pa)
!      ga = aerodynamic conductance (m/s)
!      gc = surface conductance (m/s)
!      gi = climatological conductance (m/s)
!      gr = radiative conductance (m/s) 
!   
!     NOTE:
!      1. Before calling subroutine "rotation1" and "rotation2", 
!                                    you must call subroutine "average".
!      2. Before using this program, YOU HAVE TO CHECK UP diagnostic values first.
!      3. "nc" should be more than 3.
!      4. Be careful to confirm the "Unit" of each variables
!
!
!     Reference
!
!Hong, J. and J. Kim, 2002, On processing raw data from micrometeorological field
!         experiments, Korean J. Agric. For. Meteorol., 4, 119-126
!Foken, T. and B. wichura, 1996, Tools for quality assessment of surface-based
!         flux measurements, Boundary-Layer Meteorology, 78, 83-105
!Webb, E. K., G. I. Pearman, R. Leuning, 1980, Correction of flux measurements 
!         for density effects due to heat and water vapor transfer, Quart. J. Met.
!         Soc., 106, 85-100 
!Wilczak, J. M, S. P. Oncley, S. A. Stage, 2001,Sonic anemometer tilt correction
!         algorithms, Boundary-Layer Meteorology, 99, 127-150
! 
!
!===============================================================================

      implicit none

#undef gap_filling

!===============================================================================
! User should define parameters correctly before running this program
!-------------------------------------------------------------------------------
      integer, parameter :: atime = 30              !averaging time as minutes
      integer, parameter :: ndoy = 366              !number of days in one year
      integer, parameter :: nt = (60/atime)*24*366  !number of data in one year
!      integer, parameter :: nc = 51                 !number of variables (CR5000)
!      integer, parameter :: nc2 = 5                !number of variables (CR23X)
      integer, parameter :: nc = 29                 !number of variables (CR5000)
      integer, parameter :: nc2 = 29                !number of variables (CR23X)
      integer, parameter :: nc3 = 64                !number of variables (ST data)
      integer, parameter :: nh = 3                  !number of header files 
      integer, parameter :: nq = 8                  !wind direction class for PFR
      integer, parameter :: nday = 56               !processing every "nday"
!      integer, parameter, dimension(nm) :: nmonth =  &
!         (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)

      real, parameter :: hm = 40.0             !measurement height
      real, parameter :: d0 = 14.0             !zero-plane displacement height
      real, parameter :: zm = hm-d0            !effective height (=hm-d0)
      real, parameter :: pi = 3.141592654
      real, parameter, dimension(nq+1) :: wdc = &!wind direction class for PFR
         (/ 0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0, 360.0 /)
 
      character(len=40), parameter ::  &
                            input = './L0/'
      character(len=40), parameter ::  &
                            input2 = './L0s/'
      character(len=40), parameter ::  &
                            input3 = './L0a/'
      character(len=40), parameter ::  &
                            output = './L1/'

      character(len=25), parameter, dimension(nh) :: hfile = &   !header files
                                 (/ './header_cr5000.dat',  &
                                    './header_cr23x.dat',   &
                                    './header_st.dat'    /)

      logical, parameter :: FLAG2 = .false. !Flag for the third rotation
      logical, parameter :: FLAG4 = .false.  !Flag for PFR method
      logical, parameter :: flag5 = .false. !Flag for gap-filling of meteo data

!===============================================================================

      integer :: i, j 
      integer :: idummy
      integer :: is, js, ks, ierror, nf, nn, n2, mm, nm

      real :: x(nc) 
      real :: um(3)        
      real :: su(4)
      real :: xcov(4,3)    
      real :: ta, ea, co2, rn, rflx(4), g, press, prec 
      real :: rh39m, airt39m, wind39m, qt
      real :: vdf, pa
      real :: ga, gc, gi, gr, omega
      real :: rho_a, rho_v, rho
      real :: uw, vw, uv, ut, vt, wt, uq, vq, wq, uc, vc, wc, wd
      real :: p(nq,3,3)            !P Matrix for PF method
      real :: b(nq,3), c(nq,3)     !B matrix for PF method
      real :: a(3)                 !rotation angle
      real :: wbar                 !vertical wind by WPL correction
      real :: mean                 !User-defined function
      real :: zeta                 !Monin-Obukhov lenght (z/L)
      real :: ustar                !friction velocity
      real :: ntot                 !total number of data for calculating statistics
      real :: agc                  !AGC value of LI7500
      real, allocatable, dimension(:,:) :: ax     !information for gap-filling of meteo data

      character(len=21) :: date
      character(len=21), dimension(nc+1) :: header
      character(len=21), dimension(nc2+1) :: header2
      character(len=21), dimension(nc3+1) :: header3
      character(len=60) :: fname(10), fname2(10), fname3(10), path
      character(len=60) :: path_1, path_2

      logical :: FLAG1                      !Flag to indicate the end of file
      logical :: FLAG3                      !Qualtity control flag

      path_2 = trim(trim(output)//'L1.dat')

      open(51,file=path_2,status='unknown',action='write')


!===============================================================================
!     Read the information on the directory and file name to be processed
!===============================================================================

      call fileinfo(input,fname,nf)
      call fileinfo(input2,fname2,idummy)
#ifdef gap_filling
      call fileinfo(input3,fname3,nf)
#endif

!===============================================================================
!    reading header information of data file
!===============================================================================

       call headerinfo(nc,hfile(1),header)
       call headerinfo(nc2,hfile(2),header2)
#ifdef gap_filling
       call headerinfo(nc3,hfile(3),header3)
#endif

!===============================================================================
!    End of reading the information on the directory and files to be processed
!===============================================================================
      print *, 'nf =   ', nf;pause

      DO 10 ks = 1, nf           !how many files will be processed?

       path = trim(input)//trim(fname(ks))
       open( 1,file=path,status='old',action='read',iostat=ierror)

       path = trim(input2)//trim(fname2(ks))
       open( 2,file=path,status='old',action='read',iostat=ierror)

       FLAG1 = .false.  

       if(flag5==.true.) then   ! gap-filling of meteo data
          path = trim(input3)//trim(fname3(ks))
          open( 3,file=path,status='old',action='read',iostat=ierror)

          call count(nn)
          allocate(ax(3,nn))
          call gapfilling(atime,nday,nn,nc,nc2,nc3, &
                   header,header2,header3,ax)
       endif

!-------------------------------------------------------------------------------

       IF(FLAG4==.true.) THEN   !Apply PF method

        nm = int(ndoy/nday)

        do 20 i = 1, nm          !divide the data every segment (e.g.,14days) 

           mm = 60/atime*24*nday
           if(i==nm) mm = 2*mm - 1
           nn = 0

           b(:,:) = 0.0            !initialization 
           c(:,:) = 0.0
           p(:,:,:) = 0.0
           um(:) = 0.0

           do 30 js = 1, mm

             call read(1,nc,nc2,nn,date,x,um,su,xcov,ta,ea,co2,rn,rflx,&
                   header,header2,g,prec,press,wd,ntot,agc,FLAG1,FLAG3)

             if(FLAG1==.true.) then
               mm = nn 
               backspace(1);backspace(2)
               exit
             endif

             if(FLAG3==.false.) then
                call pfr_wd(is,nq,wd,wdc)
             
                call pmatrix(um,p(is,:,:),c(is,:))
             endif

30         continue
            
           print *, 'mmm = ', mm, nn

           do js = 1, nq
              if(p(js,1,1)<=3.0) then
                 print *, 'Warning: small number of data in p matrix'
                 print *, js, p(js,1,1)
                 stop
              endif
           enddo

           do js = 1, nq
              call pfmethod(p(js,:,:),b(js,:),c(js,:)) 
           enddo

           do js = 1, mm
              backspace(1);backspace(2)
           enddo

           do 40 is = 1, mm
              call init(nc,x,um,su,xcov,uw,vw,uv,ut,vt,wt, &
                        uq,vq,wq,uc,vc,wc,a)
        
              call read(2,nc,nc2,n2,date,x,um,su,xcov,ta,ea,co2,rn,rflx,&
                    header,header2,g,prec,press,wd,ntot,agc,FLAG1,FLAG3)

              call pfr_wd(js,nq,wd,wdc)

              call pfrotation(um,su,xcov,uw,vw,uv,ut,vt,wt,uq,vq,wq,&
                     uc,vc,wc,ustar,a,p(js,:,:),b(js,:))

              call wpl(wt,wq,wc,ta,ea,co2,press,rho_a,rho_v,rho,wbar)

              call qcontrol(um,su,zm,ta,ea,uw,wt,wq,wc,zeta,ustar,FLAG3)

              call conductance(um,rn,g,ta,ea,press,uw,wt,wq,rho,vdf, &
                    ustar,ga,gc,gi,gr,omega,FLAG3)

              call outer(date,um,su,ta,ea,co2,rn,rflx,g,prec,press, &
                    ustar,uw,wt,wq,wc,ntot,agc,ga,gc,gi,gr,rho,FLAG3)

40         continue

20     continue

       ELSE                     !Apply Wesley's rotation

        DO 50


         call init(nc,x,um,su,xcov,uw,vw,uv,ut,vt,wt, &
                     uq,vq,wq,uc,vc,wc,a)

         call read(2,nc,nc2,nn,date,x,um,su,xcov,ta,ea,co2,rn,rflx,  &
               header,header2,g,prec,press,wd,ntot,agc,FLAG1,FLAG3)
         if(FLAG1==.true.) then
           exit
         endif

!
!    Application of 1st and 2nd rotation simultaneously
!
         call rotation12(um,su,xcov,uw,vw,uv,vt,wt,vq,wq,vc,wc,ustar,a) 
!
!    Application of 3rd rotation when FLAG2 is TRUE
!
         if(FLAG2==.true.) then
          call rotation3(um,su,xcov,uw,vw,uv,vt,wt,vq,wq,vc,wc,a) 
         endif

         call wpl(wt,wq,wc,ta,ea,co2,press,rho_a,rho_v,rho,wbar)

         call qcontrol(um,su,zm,ta,ea,uw,wt,wq,wc,zeta,ustar,FLAG3)
         
         call conductance(um,rn,g,ta,ea,press,uw,wt,wq,rho,vdf, &
                ustar,ga,gc,gi,gr,omega,FLAG3)
!        call potentialET

         call outer(date,um,su,ta,ea,co2,rn,rflx,g,prec,press, &
               ustar,uw,wt,wq,wc,ntot,agc,ga,gc,gi,gr,rho,FLAG3)
    
50     CONTINUE 

       ENDIF


      close(1); close(2)
10    CONTINUE

      END PROGRAM flux



!===============================================================================
!
!                            SUBROUTINEs
!
!===============================================================================

      subroutine fileinfo(input,fname,nf)

      implicit none

      integer :: ierror
      integer, intent(out) :: nf

      character(len=40), intent(in) :: input
      character(len=60), intent(out) :: fname(10)

      call system('\ls '//trim(input)//' > ./ls.dat')

      open( 3,file='./ls.dat',status='unknown',action='read')

       nf = 0 !Calculation of the number of file at each directory

       do
         read(3,'(a)',iostat=ierror) fname(nf+1)
         if(ierror==-1) exit
         nf = nf + 1
       enddo
       close(3)

!      print *, 'nf in subroutine =   ', nf
      call system('rm -f ls.dat')  !Remove "ls.dat" file.    

      return

      end subroutine fileinfo
 
!===============================================================================
      subroutine headerinfo(nc,hfile,header)

      implicit none

      integer :: i
      integer, intent(in) :: nc
      character(len=25), intent(in) :: hfile
      character(len=21), dimension(nc+1), intent(out) :: header

      open(3,file=hfile,status='old',action='read')

      do i = 1, nc+1
         read(3,*) header(i)
         header(i) = trim(header(i))
      enddo

      close(3)

      return

      end subroutine headerinfo

!===============================================================================
      subroutine count(nn)

      implicit none

      integer, intent(out) :: nn
      integer :: istat

      character(len=21) :: date


      nn = 0
      do 10
         read(1,*,iostat=istat) date
         if(istat==0) nn = nn + 1
         if(istat==-1) exit
10    continue

      rewind(1)
      return

      end subroutine count

!===============================================================================

      subroutine init(nc,x,um,su,xcov,uw,vw,uv,ut,vt,wt, &
                     uq,vq,wq,uc,vc,wc)
      implicit none  
      integer :: nc, is, js
      real :: x(nc)
      real :: um(3),su(4),xcov(4,3)
      real :: uw,vw,uv,ut,vt,wt,uq,vq,wq,uc,vc,wc
      real :: a(3)                   !rotation angles

       do is = 1, nc
        x(is) = 0.0
       enddo

       um(:)=0.0; su(:)=0.0
 
       do is = 1, 3
        a(is) = 0.0
        do js = 1, 4
         xcov(js,is) = 0.0
        enddo
       enddo
 
       uw=0.0;vw=0.0;uv=0.0
       ut=0.0;vt=0.0;wt=0.0
       uq=0.0;vq=0.0;wq=0.0
       uc=0.0;vc=0.0;wc=0.0

      RETURN
      END   


!===============================================================================

      subroutine pfr_wd(js,nq,wd,wdc)
      implicit none

      integer, intent(out) :: js
      integer, intent(in) :: nq
      real, intent(in) :: wd
      real, dimension(nq+1), intent(in) :: wdc

       if( (wd>=wdc(1)).and.(wd<wdc(2)) ) then
         js = 1
       else if( (wd>=wdc(2)).and.(wd<wdc(3)) ) then
         js = 2
       else if( (wd>=wdc(3)).and.(wd<wdc(4)) ) then
         js = 3
       else if( (wd>=wdc(4)).and.(wd<wdc(5)) ) then
         js = 4
       else if( (wd>=wdc(5)).and.(wd<wdc(6)) ) then
         js = 5
       else if( (wd>=wdc(6)).and.(wd<wdc(7)) ) then
         js = 6
       else if( (wd>=wdc(7)).and.(wd<wdc(8)) ) then
         js = 7
       else if( (wd>=wdc(8)).and.(wd<wdc(9)) ) THEN
         js = 8
       end if

      return

      end subroutine pfr_wd


!===============================================================================

      SUBROUTINE outer(date,um,su,ta,ea,co2,rn,rflx,g,prec,press, &
                  ustar,uw,wt,wq,wc,ntot,agc,ga,gc,gi,gr,rho,FLAG3)
      implicit none
       integer :: i
       real, intent(inout) :: um(3), su(4)
       real, intent(in) :: ta, ea, co2
       real, intent(in) :: rn, rflx(4), g, prec, press
       real, intent(in) :: ustar, uw, wt, wq, wc, ntot, agc
       real, intent(in) :: ga, gc, gi, gr
       real :: cp, rho, lambda, gamma 
       real :: H, LE, F_co2
       character (len=21), intent(in) :: date
       logical, intent(in) :: FLAG3

!-----------------------      NOTE     -----------------------------------------
!
! This program uses the sonic temperature instead of actual air temperature. 
!   Therefore, in calculating H and LE, heat capacity of dry air is applied 
!   with the sonic temperature (~ virtual temperature)
!
!-------------------------------------------------------------------------------

       if(FLAG3==.true.) then
          um(:) = -9999.9; su(:) = -9999.9
          H = -9999.9; LE = -9999.9; F_co2 = -9999.9
       else
          cp = 1004.67                          !Heat capacity of DRY air [K/kg/K] 
          lambda  = (2501300. - 2366.*ta)/1000. !Latent heat of varporization [J/g]

          H = rho*cp*wt
          LE = lambda*wq
          F_co2 = wc
       endif

       if(F_co2==0) then
       print *, 'date, H, LE, F_co2' 
       print *, date, H, LE, F_co2;pause
       endif

       write(51,100) date,(um(i),i=1,3),(su(i),i=1,4),ta,ea,co2, &
                      rn,(rflx(i),i=1,4),g,prec,press,ga,gc,gi,gr,&
                      H,LE,F_co2,ustar,ntot,agc,wt,wq
100    format(a21,1x,40(f16.4,1x)) 


      RETURN
      END

!===============================================================================

      SUBROUTINE read(nw,nc,nc2,nn,date,x,um,su,xcov,ta,ea,co2,rn,rflx,&
                  header,header2,g,prec,press,wd,ntot,agc,FLAG1,FLAG3)
      implicit none  
      integer, intent(in) :: nw
      integer, intent(in) :: nc, nc2
      integer, intent(out) :: nn
      integer :: nt
      integer :: is, js, ks, istat, istat2
      integer :: nuid, nsid(4), ncid(12), ntid(3)
      integer :: npid, nrn, npre(2), nwd, n_tot, n_agc, nrflx(4)

      real :: doy, x(nc), x2(nc2)
      real, intent(out) :: um(3), su(4), xcov(4,3)
      real, intent(out) :: ta, ea, co2, rn, rflx(4), g, prec, press, wd
      real, intent(out) :: ntot, agc
      real, parameter :: rv = 461.51       ! gas constant for dry air (J/kg/K)

      character(len=20), intent(out) :: date
      character(len=20) :: date2
      character(len=21), dimension(nc+1), intent(in) :: header
      character(len=21), dimension(nc2+1), intent(in) :: header2
      character(len=8) :: time

      logical, intent(out) :: FLAG1
      logical, intent(out) :: FLAG3
 

!-----------------------------------------------------------------------
!     Analysing header information
!-----------------------------------------------------------------------

       do 10 is = 1, nc+1
          if(header(is)(1:8)=='Ux_1_Avg') nuid = is-1
          if(header(is)(1:8)=='cv_Ux_Ux') nsid(1) = is-1
          if(header(is)(1:8)=='cv_Uy_Uy') nsid(2) = is-1
          if(header(is)(1:8)=='cv_Uz_Uz') nsid(3) = is-1
          if(header(is)(1:8)=='cv_Ts_Ts') nsid(4) = is-1
          if(header(is)(1:8)=='cv_Uz_Ux') ncid(1) = is-1
          if(header(is)(1:8)=='cv_Ux_Uy') ncid(2) = is-1
          if(header(is)(1:8)=='cv_Uy_Uz') ncid(3) = is-1
          if(header(is)(1:8)=='cv_Ux_Ts') ncid(4) = is-1
          if(header(is)(1:8)=='cv_Uy_Ts') ncid(5) = is-1
          if(header(is)(1:8)=='cv_Uz_Ts') ncid(6) = is-1
          if(header(is)(1:9)=='cv_Ux_co2') ncid(7) = is-1
          if(header(is)(1:9)=='cv_Uy_co2') ncid(8) = is-1
          if(header(is)(1:9)=='cv_Uz_co2') ncid(9) = is-1
          if(header(is)(1:9)=='cv_Ux_h2o') ncid(10) = is-1
          if(header(is)(1:9)=='cv_Uy_h2o') ncid(11) = is-1
          if(header(is)(1:9)=='cv_Uz_h2o') ncid(12) = is-1
          if(header(is)(1:8)=='Ts_1_Avg') ntid(1) = is-1
          if(header(is)(1:9)=='co2_1_Avg') ntid(2) = is-1
          if(header(is)(1:9)=='h2o_1_Avg') ntid(3) = is-1
          if(header(is)(1:11)=='press_1_Avg') npid = is-1
          if(header(is)(1:5)=='rain1') npre(1) = is-1
          if(header(is)(1:5)=='rain2') npre(2) = is-1
          if(header(is)(1:12)=='wnd_dir_comp') nwd = is-1
          if(header(is)(1:5)=='n_Tot') n_tot = is-1
          if(header(is)(1:3)=='agc') n_agc = is-1
10     continue
       
       do 11 is = 1, nc2+1
          if(header2(is)(1:6)=='NetRad') nrn = is-1
          if(header2(is)(1:11)=='cm3up_w_STD') nrflx(1) = is
          if(header2(is)(1:5)=='cm3dn') nrflx(2) = is-1
          if(header2(is)(1:5)=='cg3up') nrflx(3) = is-1
          if(header2(is)(1:5)=='cg3dn') nrflx(4) = is-1
11     continue
!-----------------------------------------------------------------------

       FLAG1= .false.
       FLAG3= .false.

       read(1,*,iostat=istat) date,time,doy,(x(is),is=1,nc) 
       print *, date, time, doy,istat
       date = trim(date)//' '//trim(time)
       read(2,*,iostat=istat2) date2,time,doy,(x2(is),is=1,nc2) 
       print *, date2, time, doy, istat2
       date2 = trim(date2)//' '//trim(time)
       print *, (x(is),is=1,nc);pause
       print *, (x2(is),is=1,nc2)
       
       if(istat==-1) then           !-1 in Linux and alpha machine
         print *, nn, 'End of data'
!         nn = nn + 1
         FLAG1 = .true.
!       else if (date/=date2) then
!         print *, time
!         print *, date, date2, 'inconsistent timestamp in two files'
!         pause
!         stop
       else if (istat>0) then
         print *, 'Error in reading input file 1'
         stop
       else if (istat2>0) then
         print *, 'Error in reading input file 2'
         stop
       else if((nw==1).and.(abs(x(nuid))>9999)) then
         print *, 'skip the data due to spike in sonic data'
         nn = nn + 1
         FLAG3 = .true.
         print *, date, nn, x(nuid);pause
!       else if((nw==1).and.(abs(x(ntid(2)))>9999)) then
!         print *, 'skip the data due to spike in IRGA data'
!         nn = nn + 1
!         FLAG3 = .true.
!         print *, date, nn, x(ntid(2))
       else
         nn = nn + 1
         print *, 'nn =   ', nn

         do is = 1, nc
            if(abs(x(is))>=99999) x(is)=9999.9
         enddo 

         if(abs(x(nuid))>=9999) FLAG3 = .true.
         if(abs(x(n_tot))<1000) FLAG3 = .true.

         do is = 1, 3 
           um(is)=x(nuid-1+is)
         enddo

         do is = 1, 4 
           su(is)=sqrt( x(nsid(is)) )
         enddo

         ks = 1
         DO is = 1, 4
          DO js = 1, 3 
            xcov(is,js) = x(ncid(ks))
            ks = ks + 1
          ENDDO
         ENDDO

         co2 = x(ntid(2))
         press = x(npid)*10.  ! unit conversion from "kPa" to "hPa"
         if(FLAG3==.false.) then
           ea = (x(ntid(3))/1000.)*rv*(ta+273.15)/100.   !unit is "hPa"
           ta = x(ntid(1))      ! unit is degreeC 
         else
           ea = 9999.9
           ta = 9999.9
         endif

         prec = (x(npre(1))+x(npre(2)))/2.0    !rainfall
         wd = x(nwd)                !wind direction 
         ntot = x(n_tot)            !total number of data used in calculting statistics
         agc = x(n_agc)             !AGC value from LI7500

         rn = x2(nrn)               !net radiation
         do is = 1, 4  
          rflx(is) = x2(nrflx(is))
         enddo
         g = x2(nrn)/10.            !This must be modified.!------------------------------!
       endif
      

       if(nn==0) then
        print *, 'ZERO nn:File exists, but there is no data in a file.'
        stop
       endif 
  
      RETURN
      END

!===============================================================================
      subroutine gapfilling(atime,nday,nn,nc,nc2,nc3, &
                   header,header2,header3,ax)

      implicit none

      integer, intent(in) :: atime, nday, nn, nc, nc2, nc3
      integer :: mm, nm, mm2
      integer :: istat
      integer :: is, js, ks
      integer :: ntid(3), ntid2(3), nrn, nrn2

      real, dimension(3,nn), intent(out) :: ax
      real :: doy
      real, dimension(nc,nn) :: x
      real, dimension(nc+2,nn) :: x2
      real, dimension(nc3,nn) :: x3
      real, allocatable, dimension(:) :: xx, yy
      real, allocatable, dimension(:,:) :: a, b

      character(len=21), dimension(nc+1), intent(in) :: header
      character(len=21), dimension(nc2+1), intent(in) :: header2
      character(len=21), dimension(nc3+1), intent(in) :: header3
      character(len=21) :: date
      character(len=5) :: time

      mm = (60/atime)*24*nday
      nm = int(nn/mm)
      allocate(a(3,nm)); allocate(b(3,nm))
      a = 0.0; b = 0.0

!      print *, 'mm =  ', mm, nm;pause

      do 10 js = 1, nn
         read(1,*,iostat=istat) date,time,doy,(x(is,js),is=1,nc)
         if(istat>0) then
           print *, 'reading error in subroutine gapfilling'
           stop
         endif

         read(2,*,iostat=istat) date,time,doy,(x2(is,js),is=1,nc2)
         if(istat>0) then
           print *, 'reading error in subroutine gapfilling'
           stop
         endif

         read(3,*,iostat=istat) date,time,doy,(x3(is,js),is=1,nc3)
         if(istat>0) then
           print *, 'reading error in subroutine gapfilling'
           stop
         endif
10    continue

      do is = 1, 3
         rewind(is)
      enddo

      do 20 is = 1, nc+1
         if(header(is)(1:8)=='Ts_1_Avg') ntid(1) = is-1
         if(header(is)(1:9)=='co2_1_Avg') ntid(2) = is-1
         if(header(is)(1:9)=='h2o_1_Avg') ntid(3) = is-1
20    continue

      do 21 is = 1, nc2+1
         if(header2(is)(1:6)=='NetRad') nrn = is-1
21    continue


      do 22 is = 1, nc3+1
         if(header3(is)(1:11)=='Rn_CNR1_Avg') nrn2 = is-1
         if(header3(is)(1:7)=='Ts_mean') ntid2(1) = is-1
         if(header3(is)(1:8)=='co2_mean') ntid2(2) = is-1
         if(header3(is)(1:7)=='h2o_Avg') ntid2(3) = is-1
22    continue

      do 30 is = 1, nm
         if(is==nm) then
            mm2 = nn - (is-1)*mm 
         else
            mm2 = mm
         endif

         print *, 'mm2 =  ', is,mm2, js+(is-1)*mm
         ks = 0
         do 40 js = 1, mm2
            if((abs(x(ntid(1),js+(is-1)*mm))<9999.9).and. &
               (abs(x3(ntid2(1),js+(is-1)*mm))<9999.9)) then
              ks = ks + 1
            endif  
40       continue             

         allocate(xx(ks)); allocate(yy(ks))

         ks = 0
         do 41 js = 1, mm2
            if((abs(x(ntid(1),js+(is-1)*mm))<9999.9).and. &
               (abs(x3(ntid2(1),js+(is-1)*mm))<9999.9)) then
              ks = ks + 1
              xx(ks) = x3(ntid2(1),js+(is-1)*mm)
              yy(ks) = x(ntid(1),js+(is-1)*mm)
            endif
41       continue

         print *, 'ks =  ', ks
         call lregress(ks,xx,yy,a(1,is),b(1,is))
          
!         print *, 'a, b ==  ', is, a(1,is), b(1,is);pause
30    continue

      deallocate(a); deallocate(b)
      deallocate(xx); deallocate(yy)

      return 

      end subroutine gapfilling

!===============================================================================

      subroutine lregress(n,x,y,a,b)

      implicit none
      integer, intent(in) :: n
      integer :: is
      real, intent(in) :: x(n),y(n)
      real :: t1,t2,y1,y2
      real, intent(out) :: a,b

      a=0.0
      b=0.0
      t1=0.0
      t2=0.0
      y1=0.0
      y2=0.0

      do 20 is=1,n
         t1=t1+x(is)
         t2=t2+x(is)*x(is)
         y1=y1+y(is)
         y2=y2+y(is)*x(is)
20    continue

      a=(y2-t1*y1/float(n))/(t2-t1*t1/float(n))
      b=(y1-a*t1)/float(n)

      return

      end subroutine lregress

!===============================================================================

      SUBROUTINE pmatrix(um,p,c,uid)
      IMPLICIT NONE
       INTEGER :: i    
       INTEGER :: uid
       REAL :: um(3), p(3,3), c(3)

       p(1,1) = p(1,1) + 1       !Calculation of P matrix components(3x3 matrix)
       p(1,2) = p(1,2) + um(1)
       p(1,3) = p(1,3) + um(2)
       p(2,1) = p(1,2)
       p(2,2) = p(2,2) + um(1)**2
       p(2,3) = p(2,3) + um(1)*um(2)
       p(3,1) = p(1,3)
       p(3,2) = p(2,3)
       p(3,3) = p(3,3) + um(2)**2
   
       c(1) = c(1) + um(3)   !Calculation of C matrix components(1x3 matrix)
       c(2) = c(2) + um(1)*um(3)
       c(3) = c(3) + um(2)*um(3)

!       PRINT *, 'p matrix in pfrotation in pmatrix'
!       PRINT *, (p(1,i),i=1,3)
!       PRINT *, (p(2,i),i=1,3)
!       PRINT *, (p(3,i),i=1,3)
!       PRINT *, (c(i),i=1,3), '  cmatrix '

      RETURN
      END

!===============================================================================
      SUBROUTINE pfmethod(p,b,c)
      IMPLICIT NONE
       INTEGER :: i
       INTEGER, PARAMETER :: np = 3     !P is (np x np) matrix
       INTEGER :: indx(np)
       REAL, INTENT(inout) ::p(np,np)
       REAL :: b(3), c(3), d
       REAL :: sb, cb, sa, ca

       b(:) = c(:)                 !For preservation of c values

       CALL ludcmp(p,np,np,indx,d) !LU Decomposition for the least square method
       CALL lubksb(p,np,np,indx,b) !LU Backsubstitution for getting "b" values

       p(3,1) = -b(2)/SQRT(b(2)**2 + b(3)**2 + 1) !Eqn.(42) in Wilczak et al.(2001)
       p(3,2) = -b(3)/SQRT(b(2)**2 + b(3)**2 + 1)
       p(3,3) = 1/SQRT(b(2)**2 + b(3)**2 + 1) 

       sb = -p(3,2)/SQRT(p(3,2)**2 + p(3,3)**2)
       cb =  p(3,3)/SQRT(p(3,2)**2 + p(3,3)**2) !Eqn.(44) in Wilczak et al.(2001)
       sa =  p(3,1)  
       ca = SQRT(p(3,2)**2 + p(3,3)**2) 

       p(1,1) = ca; p(1,2) = sa*sb; p(1,3) = -sa*cb ! Calculation of P matrix
       p(2,1) = 0.; p(2,2) = cb   ; p(2,3) = sb     ! P = (CD)^T Eqn.(36) 

!       PRINT *, 'p matrix in pfmethod III'
!       PRINT *, (b(i),i=1,3)
!       PRINT *, (p(1,i),i=1,3)
!       PRINT *, (p(2,i),i=1,3)
!       PRINT *, (p(3,i),i=1,3)

      RETURN  
      END     

!===============================================================================

      SUBROUTINE pfrotation(um,su,xcov,uw,vw,uv,ut,vt,wt, &
                            uq,vq,wq,uc,vc,wc,ustar,a,p,b)
      IMPLICIT NONE
       INTEGER :: is, js
       REAL, INTENT(inout) :: um(3), su(4)
       REAL, INTENT(in) :: xcov(4,3), a(3)
       REAL, INTENT(out) :: uw, vw, uv, ut, vt, wt
       REAL, INTENT(out) :: uq, vq, wq, uc, vc, wc, ustar
       REAL, INTENT(in) :: p(3,3), b(3)
       REAL :: ui, vi, wi, v(3)

!     "eta" is angle between the mean stream line and each-run stream line.

!
!      Coordinate rotation by Planar Fit Method
!      
       ui=0; vi=0; wi=0

!       ui = p(1,1)*um(1) + p(1,2)*um(2) + p(1,3)*(um(3)-b(1))
!       vi = p(2,1)*um(1) + p(2,2)*um(2) + p(2,3)*(um(3)-b(1))
!       wi = p(3,1)*um(1) + p(3,2)*um(2) + p(3,3)*(um(3)-b(1))

       ui = p(1,1)*um(1) + p(1,2)*um(2) + p(1,3)*um(3)
       vi = p(2,1)*um(1) + p(2,2)*um(2) + p(2,3)*um(3)
       wi = p(3,1)*um(1) + p(3,2)*um(2) + p(3,3)*um(3)

       um(1)=ui; um(2)=vi; um(3)=wi    !rotation of mean wind speed

       v(1) = (p(1,1)**2)*(su(1)**2) + (p(1,2)**2)*(su(2)**2)     &
            + (p(1,3)**2)*(su(3)**2)                              &
            + 2.*(p(1,1)*p(1,2)*xcov(1,2)+p(1,2)*p(1,3)*xcov(1,3) &
            + p(1,1)*p(1,3)*xcov(1,1)) + (p(1,3)**2)*(b(1)**2)

       v(2) = (p(2,1)**2)*(su(1)**2) + (p(2,2)**2)*(su(2)**2)     &
            + (p(2,3)**2)*(su(3)**2)                              &
            + 2.*(p(2,1)*p(2,2)*xcov(1,2)+p(2,2)*p(2,3)*xcov(1,3) &
            + p(2,1)*p(1,3)*xcov(1,1)) + (p(2,3)**2)*(b(1)**2)

       v(3) = (p(3,1)**2)*(su(1)**2) + (p(3,2)**2)*(su(2)**2)     &
            + (p(3,3)**2)*(su(3)**2)                              &
            + 2.*(p(3,1)*p(3,2)*xcov(1,2)+p(3,2)*p(3,3)*xcov(1,3) &
            + p(3,1)*p(3,3)*xcov(1,1)) + (p(3,3)**2)*(b(1)**2)

!       print *, 'b1, b2, b3  ', b(1), b(2), b(3)
!       print *, 'su1, su2, su3  ', su(1), su(2), su(3)
!       print *, 'p33 = ', p(3,3)
!       print *, v(1), v(2), v(3)
!       print *, p(3,1)*p(3,2)*xcov(1,2)
!       print *, p(3,2)*p(3,3)*xcov(1,3)
!       print *, p(3,1)*p(3,3)*xcov(1,1)
!
!      Rotation of vertival turbulent fluxes
! 
       uw = p(1,1)*p(3,1)*su(1)**2+p(1,2)*p(3,2)*su(2)**2 &
         + p(1,3)*p(3,3)*su(3)**2 &
         + (p(1,1)*p(3,2)+p(1,2)*p(3,1))*xcov(1,2)+(p(1,1)*p(3,3)      &
         + p(1,3)*p(3,1))*xcov(1,1)                 &
         + (p(1,2)*p(3,3)+p(1,3)*p(3,2))*xcov(1,3)  &
         + p(1,3)*p(3,3)*(b(3)**2)
       vw = p(2,1)*p(3,1)*su(1)**2+p(2,2)*p(3,2)*su(2)**2  &
         + p(2,3)*p(3,3)*su(3)**2&
         + (p(2,1)*p(3,2)+p(2,2)*p(3,1))*xcov(1,2) &
         + (p(2,1)*p(3,3)+p(2,3)*p(3,1))*xcov(1,1) &
         + (p(2,2)*p(3,3)+p(2,3)*p(3,2))*xcov(1,3) &
         + p(2,3)*p(3,3)*(b(3)**2)

       wt = p(3,1)*xcov(2,1)+p(3,2)*xcov(2,2)+p(3,3)*xcov(2,3)
       wq = p(3,1)*xcov(4,1)+p(3,2)*xcov(4,2)+p(3,3)*xcov(4,3)
       wc = p(3,1)*xcov(3,1)+p(3,2)*xcov(3,2)+p(3,3)*xcov(3,3)

       su(1)=SQRT(v(1)); su(2)=SQRT(v(2)); su(3)=SQRT(v(3)) 

       ustar = sqrt( sqrt( uw**2 + vw**2 ) )

!       print *, p(3,1), p(3,2), p(3,3)
!       print *, 'aaaa', wc, xcov(3,1),xcov(3,3), wq, xcov(4,3)
!       print *, xcov(4,1),xcov(4,2), xcov(4,3)
!       print *, 'aaaa', wc, xcov(3,1),xcov(3,3), 2450*wq, 2450*xcov(4,3);pause

       call rotation1(um,su,uw,vw,uv,ut,vt,wt,uq,vq,wq,uc,vc,wc,a)

      RETURN
      END


!===============================================================================
   
      SUBROUTINE rotation1(um,su,uw,vw,uv, &
                           ut,vt,wt,uq,vq,wq,uc,vc,wc,a)
      IMPLICIT NONE
       INTEGER :: is, js
       REAL, INTENT(inout) :: um(3), su(4)
       REAL :: xcov(4,3)
       REAL, INTENT(inout) :: uw, vw, uv
       REAL, INTENT(inout) :: ut, vt, wt, uq, vq, wq, uc, vc, wc, a(3)
       REAL :: uw0, vw0, uv0, ut0, vt0, uq0, vq0, uc0, vc0, su0(4)
       REAL :: ce, se
       REAL :: ce2, se2
       REAL :: u(3)
       REAL :: v(3)

       
       a(1) = atan( um(2)/um(1) )
       ce=um(1)/(sqrt(um(1)**2 + um(2)**2))
       se=um(2)/(sqrt(um(1)**2 + um(2)**2))

       ce2=ce**2; se2=se**2

!     Calculation of variance from standard deviation
       su0(1)=su(1)**2; su0(2)=su(2)**2
       uw0=uw; vw0=vw; uv0=uv; ut0=ut; vt0=vt
       uq0=uq; vq0=vq; uc0=uc; vc0=vc
!       
       u(1) =  um(1)*ce + um(2)*se 
       u(2) = -um(1)*se + um(2)*ce

       v(1) = su0(1)*ce2 + su0(2)*se2 +  2.*uv0*ce*se
       v(2) = su0(1)*se2 + su0(2)*ce2 - 2.*uv0*ce*se

       uw = uw0*ce + vw0*se 
       uv = uv0*(ce2-se2) - su0(1)*ce*se + su0(2)*ce*se
       vw = vw0*ce - uw0*se

       vt = vt0*ce - ut0*se 
       vc = vc0*ce - uc0*se 
       vq = vq0*ce - uq0*se 

       um(1)=u(1); um(2)=u(2)

!       print *, 'um, vm =', um(1), um(2)
!       print *, 'su0, sv0 = ', su0(1), su0(2)
!       print *, 'su, sv = ', su(1), su(2)
!       print *, 'ce, se = ', ce, se
!       print *, 'ce2, se2 = ', ce2, se2
!       print *, 'uv= ', uv
!       print *, 'v1, v2 = ', v(1), v(2)
       su(1)=SQRT(v(1)); su(2)=SQRT(v(2))
            
      RETURN
      END

!===============================================================================
   
      SUBROUTINE rotation12(um,su,xcov,uw,vw,uv,vt,wt,vq,wq, &
                           vc,wc,ustar,a)
      IMPLICIT NONE
      INTEGER :: is, js
      REAL :: um(3), su(4)
      REAL :: xcov(4,3)
      REAL,INTENT(OUT) :: uw, vw, uv, wt, wq, wc
      REAL,INTENT(OUT) :: vt, vq, vc, ustar, a(3)  
      REAL :: ce, se, ct, st
      REAL :: ce2, se2, ct2, st2
      REAL :: u(3)
      REAL :: v(3)

       a(1) = atan( um(2)/um(1) )
       a(2) = atan( um(3)/SQRT(um(1)**2 + um(2)**2) )
       ce=um(1)/(sqrt(um(1)**2 + um(2)**2))
       se=um(2)/(sqrt(um(1)**2 + um(2)**2))
       ct=(sqrt(um(1)**2 + um(2)**2)) &
             /(sqrt(um(1)**2 + um(2)**2 + um(3)**2))
       st=um(3)/(sqrt(um(1)**2 + um(2)**2 + um(3)**2))

       ce2=ce**2; se2=se**2; ct2=ct**2; st2=st**2

!     Calculation of variance from standard deviation
       su(1)=su(1)**2; su(2)=su(2)**2; su(3)=su(3)**2 
!       
       u(1) =  um(1)*ct*ce + um(2)*ct*se + um(3)*st
       u(2) = -um(1)*se    + um(2)*ce
       u(3) = -um(1)*st*ce - um(2)*st*se + um(3)*ct

       v(1) = su(1)*ct2*ce2 + su(2)*ct2*se2 + su(3)*st2 &
          + 2.*xcov(1,2)*ct2*ce*se&
          + 2.*xcov(1,1)*ct*st*ce + 2.*xcov(1,3)*ct*st*se 
       v(2) = su(1)*se2 + su(2)*ce2 - 2.*xcov(1,2)*ce*se
       v(3) = su(1)*st2*ce2 + su(2)*st2*se2 + su(3)*ct2  &
          + 2.*xcov(1,2)*st2*ce*se - 2.*xcov(1,1)*ct*st*ce &
          - 2.*xcov(1,3)*ct*st*se 
       uw = xcov(1,1)*ce*(ct2-st2) - 2.*xcov(1,2)*ct*st*ce*se&
          + xcov(1,3)*se*(ct2-st2) - su(1)*ct*st*ce2 - su(2)*ct*st*se2&
          + su(3)*ct*st
       uv = xcov(1,2)*ct*(ce2-se2) + xcov(1,3)*st*ce - xcov(1,3)*st*se&
          - su(1)*ct*ce*se + su(2)*ct*ce*se
       vw = xcov(1,3)*ct*ce - xcov(1,1)*ct*se-xcov(1,2)*st*(ce2-se2)&
          + su(1)*st*ce*se - su(2)*st*ce*se

       wt = xcov(2,3)*ct - xcov(2,1)*st*ce - xcov(2,2)*st*se
       wc = xcov(3,3)*ct - xcov(3,1)*st*ce - xcov(3,2)*st*se
       wq = xcov(4,3)*ct - xcov(4,1)*st*ce - xcov(4,2)*st*se 

       vt = xcov(2,2)*ce - xcov(2,1)*se 
       vc = xcov(3,2)*ce - xcov(3,1)*se 
       vq = xcov(4,2)*ce - xcov(4,1)*se 

       um(1)=u(1); um(2)=u(2); um(3)=u(3)
       su(1)=SQRT(v(1)); su(2)=SQRT(v(2)); su(3)=SQRT(v(3))
       ustar = sqrt( sqrt(uw**2+vw**2) )
            
      RETURN
      END

!===============================================================================

      SUBROUTINE rotation3(um,su,xcov,uw,vw,uv,vt,wt,vq,wq,vc,wc,a) 
      IMPLICIT NONE
      INTEGER :: is, js
      REAL :: um(3), su(4), xcov(4,3)
      REAL,INTENT(INOUT) :: uw, vw, uv, wt, wq, wc
      REAL,INTENT(INOUT) :: vt, vq, vc, a(3)
      REAL :: sp, cp
      REAL :: uw_t, vw_t
 
       a(3) = 1./2.*atan( 2.*vw/(su(2)**2 - su(3)**2) )

       sp = sin(a(3)); cp = cos(a(3))

       uw_t = -uv*sp + uw*cp
       vw_t = (su(3)**2-su(2)**2)*cp*sp + vw*(cp**2-sp**2)
       uw = uw_t; vw = vw_t
       wt = -vt*sp + wt*cp
       wc = -vc*sp + wc*cp
       wq = -vq*sp + wq*cp

      RETURN
      END
 
!===============================================================================

      SUBROUTINE wpl(wt,wq,wc,ta,ea,co2,press,rho_a,rho_v,rho,wbar)
      IMPLICIT NONE
      INTEGER:: is
       REAL :: wt, wq, wc, ta, ea, co2, wbar
       REAL :: tak                    !Air temperature (unit is Kelvin)
       REAL :: t(3)                   !Temporary variables for WPL correction
       REAL :: Pd                     !Dry air pressure
       REAL :: Press                  !Air pressure
!
!     If there are pressure and temperature measurements, 
!      we can get the information on air density accurately from ideal gas law.
!
       REAL, PARAMETER :: ma = 28.964   !mean molecular weight (g/ mol) 
       REAL, PARAMETER :: mv = 18.015   !molecular weight of water vapor (g/mol)
       REAL, PARAMETER :: Rd = 287.04   !gas constant for dry air (J/kg/K)
       REAL, PARAMETER :: Rv = 461.51   !gas constant for water vapor (J/kg/K)
       REAL, PARAMETER :: mu = ma/mv

       REAL :: rho_a     !dry air density  = Pa/(Rd Ta)
       REAL :: rho_v     !water vapor density = ea/(Rv Ta)
       REAL :: rho_c     !CO2 density from observation (mg/m3)
       REAL :: rho       !air density (=rho_a + rho_v)(Kg/m3)
       REAL :: sigma
      
       tak = ta + 273.15           !Unit conversion from degree to Kelvin
       Pd = Press - ea             !Dry air pressure (hPa)
       rho_c = co2/1000./1000.     !Unit conversion from (mg/m3) to (kg/m3) 
       rho_a = (Pd*100.)/(Rd*tak)
       rho_v = ea*100./(Rv*tak)
       rho = Press*100./Rd/tak
       sigma = rho_v/rho_a

!    
!     Unit conversions for WPL correction
!
       wq = wq/1000.                     !from g/m2/s to kg/m2/s
       wc = wc/1000./1000.               !from mg/m2/s to kg/m2/s

!       print *, 'wq= ', 2400*wq*1000
!       print *, 'wc= ', wc*1000*1000
!
!
!     WPL correction for latent heat flux and CO2 flux
!
       wbar = mu/(1.+mu*sigma)*wq/rho_a+wt*tak       !Eqn. (43) in Webb et al.(1980) 
       wq = (1.+mu*sigma)*(wq + (rho_v/tak)*wt)      !Eqn. (25) 
       wc = wc+mu*(rho_c/rho_a)*wq + (1.+mu*sigma)*(rho_c/tak)*wt !Eqn.(24)

!       print *, 'rho_v/tak*wt= ', (wq+(rho_v/tak)*wt)*1000
!       print *, 'sigma= ', sigma
!       print *, 'ta= ', ta
!       print *, 'ea= ', ea
!       print *, 'rho_c= ', rho_c
!       print *, 'rho_a= ', rho_a
!       print *, 'rho_v= ', rho_v
!       print *, 'wt= ', wt
!       print *, 'wq= ', 2400*wq*1000
!       print *, 'wc= ', wc*1000*1000
!       pause

!     Unit recovery
       wc = wc * (10.**6)                 !from kg/m2/s to mg/m2/s          
       wq = wq * (10.**3)                 !from kg/m2/s to  g/m2/s          

      RETURN
      END
  
!===============================================================================

      SUBROUTINE conductance(um,rn,g,ta,ea,press,uw,wt,wq,rho,vdf, &
                   ustar,ga,gc,gi,gr,omega,FLAG3)
      IMPLICIT NONE
      real, intent(in) :: um(3), rn, g, uw, wt, wq, rho, ustar
      real, intent(out) :: vdf, ga, gc, gi, gr
      REAL :: av        !Available energy (=H+LE)
      REAL :: ta, ea, press
      real :: es        !saturated vapor pressure
      real :: delta, epsilon 
      REAL :: a(4)
      REAL :: tr
      REAL :: omega     !Decoupling factor   
      REAL :: cp        !specific heat capacity of air
      REAL :: lambda    !heat of vaporization
      REAL :: gamma     !psychrometric constant
      REAL :: H         !sensible heat flux (W/m2)
      REAL :: LE        !latent heat flux (w/m2)
      REAL :: q         !specific humidity
      real :: tak
      REAL, PARAMETER :: stb = 5.67E-08  !Stefan-Boltzmann constant
      REAL, PARAMETER :: ep = 0.622
      logical, intent(in) :: FLAG3

      DATA a / 13.3185, 1.9760, 0.6445, 0.1299/ !coefficients for calculating es
!
!  Formular for lambda is same with one in SiB2 code.
!
!  Careful treatment of units of each variable
!
       q = ep*(ea*100.)/((ep-1.)*(ea*100.)+press)  !calculation of mixing ratio
       cp = 1004.67*(1+0.87*q)                     !heat capacity
       lambda  = (2501300. - 2366.*ta)             !unit is J/kg
       gamma =  (cp*press)/(ep*lambda) 
       H = rho*cp*wt
       LE = rho*(lambda/1000.)*wq
       av = H + LE                                 !available energy

!  Calculation of saturation vapor pressure and each change to temperature
!  See the book of Brutsaert (Evaporation into the Atmosphere)
!
       tak = ta + 273.15
       tr = 1.0-(373.15/tak)
       es = 1013.25*EXP( a(1)*tr-a(2)*(tr**2)-(a(3)*tr**3)-a(4)*(tr**4))
       delta = 373.15*es/(tak**2)*&
             ( a(1) - 2.*a(2)*tr - 3.*a(3)*(tr**2) - 4.*a(4)*(tr**3) )
       vdf = (es - ea)
       epsilon = delta/gamma
!       print *, ea/es*100. 
!       print *, 'vdf= ', vdf;pause 

!
!  Calculation of conductances
!
       if(FLAG3==.true.) then                            !Data quality check-up
         gi=9999.9; ga=9999.9; gc=9999.9; gr=9999.9;omega=9999.9
       else 
         gi = (gamma*av)/(rho*cp*vdf)
         ga = 1./( um(1)/ustar**2 + 4.626/ustar )       !aerodynamic conductance      
!        gc = 1./( rho*vdf/LE*lambda & 
!           + 1./ga*(epsilon*(H+LE)/LE-(epsilon+1.)) )
         gc = ( (1.+H/LE)*(epsilon+gi/ga)-epsilon-1. )/ ga
         gc = 1./gc                                     !canopy conductance 
         omega = (1.+epsilon)/(1.+epsilon+ga/gc)
         gr = 4.*stb*(ta**3)/rho/cp

!        print *, 'omega = ', omega;pause
       endif
 
      RETURN
      END

!===============================================================================

      SUBROUTINE qcontrol(um,su,zm,ta,ea,uw,wt,wq,wc,zeta,ustar,FLAG3)
      IMPLICIT NONE
      real, intent(in):: um(3), su(4), zm, ta, ea
      real, intent(in) :: uw, wt,wq, wc  
      real, intent(in) :: ustar
      real, intent(out) :: zeta       !=z/L
      real :: L                       !Obukhov length
      real :: itc, itc_obs            !Integral turbulence characteristics
      real :: c(2), difference
      real :: tak
      real, parameter :: gc = 9.81    !gravity constant
      real, parameter :: kc = 0.4     !von Karman constant
      real, parameter :: crv = 1.0    !Critical value of I.T.C.
      logical, intent(inout) :: FLAG3

!       IF(uw < 0.0) THEN
         L = -ustar**3/( kc*gc/tak*wt )
         zeta = zm/L
!       ELSE
!         PRINT *, 'Positive momentum flux happen!'
!         zeta = 9999.9
!         FLAG3 = .TRUE.
!         PRINT *, 'bad quality  ', 'positive momentum flux'
!         return
!       ENDIF

!       PRINT *, 'uw, wt, wq, wt, wc, zeta'
!       PRINT *, uw, wt, wq, wt, wc, zeta; pause

!       IF(um(1) < 0.05) THEN
!        FLAG3 = .TRUE.
!        PRINT *, 'Low wind speed:  < 1.0 m/s  ' 
!        return
!       ENDIF

!       IF(zeta>-0.4) THEN
!        c(1) = 1.4; c(2) = 0.0
!       ELSE
!        c(1) = 1.9; c(2) = 1.0/3.0
!       ENDIF

!       IF(zeta < 1.0) THEN
!         itc = c(1)*(1.0-zeta)**c(2)
!       ELSE
!         itc = c(1)
!       ENDIF

!       itc_obs = su(3)/SQRT(-uw)
!       difference = ABS((itc-itc_obs)/itc)

!       IF(difference > crv) THEN
!        FLAG3 = .true.
!        PRINT *, 'bad quality  ', 'out of range of I.T.C of w'
!        PRINT *, 'itc = ', itc, 'itc_obs = ', itc_obs
!        return
!       ENDIF


      RETURN
      END subroutine qcontrol
  
!=========================================================================

      SUBROUTINE ludcmp(a,n,np,indx,d)
      IMPLICIT NONE
      INTEGER :: n, np, indx(n)
      INTEGER :: i, j, k, imax
      INTEGER, PARAMETER :: NMAX = 500

      REAL :: aamax, dum, sum, vv(NMAX)
      REAL :: d, a(np,np)
      REAL, PARAMETER :: TINY = 1.0E-20

       d = 1.0

       DO i = 1, n
         aamax = 0.0
         DO j = 1, n
           IF(ABS(a(i,j)).gt.aamax) THEN
            aamax = ABS(a(i,j))
           ENDIF
         ENDDO
         IF(aamax==0.) THEN
          PRINT *, 'singular matrix in LU decomposition'
          pause
         ENDIF
       ENDDO

       DO j=1,n
        DO i=1,j-1
          sum=a(i,j)
          DO k=1,i-1
            sum=sum-a(i,k)*a(k,j)
         ENDDO
          a(i,j)=sum
        ENDDO

        aamax=0.
        DO i=j,n
          sum=a(i,j)
          DO k=1,j-1
            sum=sum-a(i,k)*a(k,j)
          ENDDO
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          IF(dum.ge.aamax) THEN
            imax=i
            aamax=dum
          ENDIF
        ENDDO
        IF(j.ne.imax) THEN
          DO k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
          ENDDO
          d=-d
          vv(imax)=vv(j)
        ENDIF
        indx(j)=imax
        IF(a(j,j)==0.) a(j,j)=TINY
        IF(j.ne.n) THEN
          dum=1./a(j,j)
          DO i=j+1,n
            a(i,j)=a(i,j)*dum
         ENDDO
        ENDIF
       ENDDO

      RETURN
      END


!=========================================================================

      SUBROUTINE lubksb(a,n,np,indx,b)
      IMPLICIT NONE
       INTEGER n,np,indx(n)
       REAL a(np,np),b(n)
       INTEGER i,ii,j,ll
       REAL sum

       ii=0
       DO i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        IF(ii.ne.0) THEN
          DO j=ii,i-1
            sum=sum-a(i,j)*b(j)
          ENDDO
        ELSE IF(sum.ne.0.) THEN
          ii=i
        ENDIF
        b(i)=sum
       ENDDO

       DO i=n,1,-1
        sum=b(i)
        DO j=i+1,n
          sum=sum-a(i,j)*b(j)
        ENDDO
        b(i)=sum/a(i,i)
       ENDDO

      RETURN
      END


!=========================================================================
!=========================================================================
!
!     User-defined Function
!
!=========================================================================
!=========================================================================


!=========================================================================
      REAL FUNCTION  mean(x,N)
!
!     PURPOSE:
!     To evaluate mean value
!
!=========================================================================
!
!    x = variable
!    N = the total number of data
!
!=========================================================================
!
      IMPLICIT NONE
!
!=========================================================================

      INTEGER :: is
      INTEGER, INTENT(IN) :: N
      REAL, DIMENSION(N), INTENT(IN) :: x
      REAL :: sum

      sum = 0.0

      DO is = 1, N
       sum = sum + x(is)
      END DO

      mean = sum / N

      END FUNCTION mean

!=========================================================================








