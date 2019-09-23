! FITTER_MC, Neil Drummond, 18/4/06.

! This is a general-purpose utility for using nl2sno to fit a model to a set
! of data.  The name of the data file can be specified as a command-line
! argument.  Error bars should be supplied, as the data will be Monte-Carlo
! simulated in order to estimate error bars on the fitted parameters and
! properties expressed in terms of the parameters.

! If you wish to change the model, you only need to alter the three subroutines
! containing the comment

! YOU NEED TO CHANGE THIS SUBROUTINE IF YOU HAVE A NEW MODEL.

! The program uses the iargc and getarg extensions to F90.  If these are not
! available then please uncomment the dummy versions below.

!!$SUBROUTINE getarg(i,c)
!!$  IMPLICIT NONE
!!$  INTEGER,INTENT(in) :: i
!!$  CHARACTER(*),INTENT(out) :: c
!!$  c=''
!!$END SUBROUTINE getarg
!!$ 
!!$INTEGER FUNCTION iargc()
!!$  IMPLICIT NONE
!!$  iargc=0
!!$END FUNCTION iargc


MODULE utils
  ! Miscellaneous utilities.
  USE machine_constants,ONLY : dp ! Double-precision type.
  IMPLICIT NONE


CONTAINS


  SUBROUTINE errstop(sub,message)
    ! Report an error and stop.
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: sub,message
    WRITE(*,*)
    WRITE(*,*)'ERROR in subroutine '//TRIM(ADJUSTL(sub))//'.'
    WRITE(*,*)
    CALL wordwrap(TRIM(ADJUSTL(message)))
    WRITE(*,*)
    STOP
  END SUBROUTINE errstop


  SUBROUTINE wordwrap(text,unit_in,linelength_in)
    ! This subroutine prints out the contents of the character string 'text',
    ! ensuring that line breaks only occur at space characters.  The output
    ! is written to unit unit_in if this parameter is supplied; otherwise the
    ! output is written to unit o.  The maximum length of each line is given
    ! by linelength_in if this is supplied; otherwise the default line length
    ! is 79 characters.
    IMPLICIT NONE
    INTEGER,INTENT(in),OPTIONAL :: unit_in,linelength_in
    CHARACTER(*),INTENT(in) :: text
    CHARACTER(260) :: temp
    INTEGER :: i,unit,lentext,startpoint,stoppoint,lastpos,linelength
    IF(PRESENT(unit_in))THEN
      unit=unit_in
    ELSE
      unit=6
    ENDIF ! unit supplied.
    lentext=LEN_TRIM(text)
    IF(lentext<1)THEN
      WRITE(unit,*)
      RETURN
    ENDIF ! No text
    IF(PRESENT(linelength_in))THEN
      IF(linelength_in>=2)THEN
        linelength=linelength_in
      ELSE
        linelength=2
      ENDIF ! sensible line-length supplied.
    ELSE
      linelength=79
    ENDIF ! linelength present.
    startpoint=1
    DO i=1,HUGE(1)
      stoppoint=startpoint+linelength-1
      IF(stoppoint<=lentext)THEN
        lastpos=INDEX(TRIM(text(startpoint:stoppoint))," ",.TRUE.)
        IF(lastpos>0)stoppoint=startpoint+lastpos-1
      ELSE
        stoppoint=lentext
      ENDIF ! stoppoint <= length of text
      IF(i==1)THEN
        ! Allow the user to indent the first line, if (s)he wishes.
        temp=text(startpoint:stoppoint) ! or else pathscale f90 fails to compile
        WRITE(unit,*)TRIM(temp)
      ELSE
        temp=text(startpoint:stoppoint) ! or else pathscale f90 fails to compile
        WRITE(unit,*)TRIM(ADJUSTL(temp))
      ENDIF ! i=1
      IF(stoppoint==lentext)THEN
        EXIT
      ELSE
        startpoint=stoppoint+1
      ENDIF ! Finished text?
    ENDDO ! Loop over lines.
  END SUBROUTINE wordwrap


  SUBROUTINE swapreal(a,b)
    ! This subroutine swaps two real numbers a and b.
    IMPLICIT NONE
    REAL(dp),INTENT(inout) :: a,b
    REAL(dp) :: temp
    temp=a ; a=b ; b=temp
  END SUBROUTINE swapreal


  LOGICAL FUNCTION isdataline(c)
    ! This function returns T if & only if the first word of c is a number.
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: c
    INTEGER :: ierr
    REAL(dp) :: r
    READ(c,*,iostat=ierr)r
    isdataline=(ierr==0)
  END FUNCTION isdataline


  CHARACTER(12) FUNCTION i2s(n)
    ! Convert integers to left justified strings that can be printed in the
    ! middle of a sentence without introducing large amounts of white space.
    INTEGER,INTENT(in) :: n
    INTEGER i,j
    INTEGER,PARAMETER :: ichar0=ICHAR('0')
    i2s=''
    i=ABS(n)
    DO j=LEN(i2s),1,-1
      i2s(j:j)=ACHAR(ichar0+MOD(i,10))
      i=i/10 ; IF(i==0)EXIT
    ENDDO ! j
    IF(n<0)THEN
      i2s='-'//ADJUSTL(i2s)
    ELSE
      i2s=ADJUSTL(i2s)
    ENDIF ! n<0
  END FUNCTION i2s


  CHARACTER(72) FUNCTION write_mean(av,std_err_in_mean,err_prec_in)
    ! Write out a mean value with the standard error in the mean in the
    ! form av(std_err_in_mean), e.g. 0.123546(7).  err_prec_in specifies
    ! the number of digits of precision to which the error should be
    ! quoted (by default 1).
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: av,std_err_in_mean
    INTEGER,INTENT(in),OPTIONAL :: err_prec_in
    INTEGER :: lowest_digit_to_quote,err_quote,err_prec,int_part,dec_part,i
    INTEGER,PARAMETER :: err_prec_default=1
    REAL(dp) :: av_quote
    CHARACTER(1) sgn
    CHARACTER(72) zero_pad

    IF(std_err_in_mean<=0.d0)THEN
      write_mean='ERROR: NON-POSITIVE ERROR BAR!!!'
      RETURN
    ENDIF ! Error is negative

    IF(PRESENT(err_prec_in))THEN
      IF(err_prec_in>=1)THEN
        err_prec=err_prec_in
      ELSE
        write_mean='ERROR: NON-POSITIVE PRECISION!!!'
        RETURN
      ENDIF ! err_prec_in sensible.
    ELSE
      err_prec=err_prec_default
    ENDIF ! Accuracy of error supplied.

    ! Work out lowest digit of precision that should be retained in the
    ! mean (i.e. the digit in terms of which the error is specified).
    ! Calculate the error in terms of this digit and round.
    lowest_digit_to_quote=FLOOR(LOG(std_err_in_mean)/LOG(10.d0)) &
      &+1-err_prec
    err_quote=NINT(std_err_in_mean*10.d0**DBLE(-lowest_digit_to_quote))
    IF(err_quote==10**err_prec)THEN
      lowest_digit_to_quote=lowest_digit_to_quote+1
      err_quote=err_quote/10
    ENDIF ! err_quote rounds up to next figure.

    IF(err_quote>=10**err_prec.OR.err_quote<10**(err_prec-1))THEN
      write_mean='ERROR: BUG IN WRITE_MEAN!!!'
      RETURN
    ENDIF ! Check error is in range.

    ! Truncate the mean to the relevant precision.  Establish its sign,
    ! then take the absolute value and work out the integer part.
    av_quote=ANINT(av*10.d0**DBLE(-lowest_digit_to_quote)) &
      &*10.d0**DBLE(lowest_digit_to_quote)
    IF(av_quote<0.d0)THEN
      sgn='-'
      av_quote=-av_quote
    ELSE
      sgn=''
    ENDIF ! Sign
    IF(AINT(av_quote)>DBLE(HUGE(1)))THEN
      write_mean='ERROR: NUMBERS ARE TOO LARGE IN WRITE_MEAN!'
      RETURN
    ENDIF ! Vast number
    int_part=FLOOR(av_quote)

    IF(lowest_digit_to_quote<0)THEN
      ! If the error is in a decimal place then construct string using
      ! integer part and decimal part, noting that the latter may need to
      ! be padded with zeros, e.g. if we want "0001" rather than "1".
      IF(ANINT((av_quote-DBLE(int_part)) &
        &*10.d0**DBLE(-lowest_digit_to_quote))>DBLE(HUGE(1)))THEN
        write_mean='ERROR: NUMBERS ARE TOO LARGE IN WRITE_MEAN!'
        RETURN
      ENDIF ! Vast number
      dec_part=NINT((av_quote-DBLE(int_part)) &
        &*10.d0**DBLE(-lowest_digit_to_quote))
      zero_pad=''
      IF(dec_part<0)THEN
        write_mean='ERROR: BUG IN WRITE_MEAN! (2)'
        RETURN
      ENDIF ! dec
      DO i=1,-lowest_digit_to_quote-no_digits_int(dec_part)
        zero_pad(i:i)='0'
      ENDDO ! i
      write_mean=sgn//TRIM(i2s(int_part))//'.'//TRIM(zero_pad) &
        &//TRIM(i2s(dec_part))//'('//TRIM(i2s(err_quote))//')'
    ELSE
      ! If the error is in a figure above the decimal point then, of
      ! course, we don't have to worry about a decimal part.
      write_mean=sgn//TRIM(i2s(int_part))//'(' &
        &//TRIM(i2s(err_quote*10**lowest_digit_to_quote))//')'
    ENDIF ! lowest_digit_to_quote<0

  END FUNCTION write_mean


  INTEGER FUNCTION no_digits_int(i)
    ! Calculate the number of digits in integer i.  For i>0 this
    ! should be floor(log(i)/log(10))+1, but sometimes rounding errors
    ! cause this expression to give the wrong result.  
    IMPLICIT NONE
    INTEGER,INTENT(in) :: i
    INTEGER :: j,k
    j=i
    k=1
    DO
      j=j/10
      IF(j==0)EXIT
      k=k+1
    ENDDO
    no_digits_int=k
  END FUNCTION no_digits_int


  REAL(dp) FUNCTION gammq(a,x)
    ! This function returns Q(a,x)=Gamma_incomplete(a,x)/Gamma(a).
    ! Adapted from Numerical Recipes.
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: a,x
    IF(x<0.d0.OR.a<=0.d0)CALL errstop('GAMMQ','Undefined.')
    IF(x<a+1.d0)THEN
      gammq=1.d0-gser(a,x)
    ELSE
      gammq=gcf(a,x)
    ENDIF
  END FUNCTION gammq


  REAL(dp) FUNCTION gser(a,x)
    ! This function returns the Gamma P(a,x) function, evaluated as a series
    ! expansion.  Adapted from Numerical Recipes.
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: a,x
    REAL(dp),PARAMETER :: eps=6.d-16
    INTEGER,PARAMETER :: itmax=10000
    INTEGER :: n
    REAL(dp) :: ap,del,sum
    IF(x<=0.d0)THEN
      IF(x<0.d0)CALL errstop('GSER','x<0.')
      gser=0.d0
      RETURN
    ENDIF ! x<=0
    ap=a
    sum=1.d0/a
    del=sum
    DO n=1,itmax
      ap=ap+1.d0
      del=del*x/ap
      sum=sum+del
      IF(ABS(del)<ABS(sum)*eps)THEN
        gser=sum*EXP(-x+a*LOG(x)-gammln(a))
        RETURN
      ENDIF ! Converged
    ENDDO ! n
    CALL errstop('GSER','Failed to converge.')
  END FUNCTION gser


  REAL(dp) FUNCTION gcf(a,x)
    ! This function returns the Gamma Q(a,x) function, evaluated as a
    ! continued fraction.  Adapted from Numerical Recipes.
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: a,x
    REAL(dp),PARAMETER :: eps=6.d-16,fpmin=1.d-300
    INTEGER,PARAMETER :: itmax=10000
    INTEGER :: i
    REAL(dp) :: an,b,c,d,del,h
    b=x+1.d0-a
    c=1.d0/fpmin
    d=1.d0/b
    h=d
    DO i=1,itmax
      an=-DBLE(i)*(DBLE(i)-a)
      b=b+2.d0
      d=an*d+b
      IF(ABS(d)<fpmin)d=fpmin
      c=b+an/c
      IF(ABS(c)<fpmin)c=fpmin
      d=1.d0/d
      del=d*c
      h=h*del
      IF(ABS(del-1.d0)<eps)THEN
        gcf=EXP(-x+a*LOG(x)-gammln(a))*h
        RETURN
      ENDIF ! Converged
    ENDDO ! i
    CALL errstop('GCF','Failed to converge.')
  END FUNCTION gcf


  REAL(dp) FUNCTION gammln(xx)
    ! This function returns the logarithm of the Gamma function.  Adapted from
    ! Numerical Recipes.
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: xx
    INTEGER :: j
    REAL(dp) :: ser,tmp,x,y
    REAL(dp),PARAMETER :: cof(6)=(/76.18009172947146d0, &
      &-86.50532032941677d0,24.01409824083091d0,-1.231739572450155d0, &
      &0.1208650973866179d-2,-0.5395239384953d-5/),stp=2.5066282746310005d0
    x=xx
    y=x
    tmp=x+5.5d0
    tmp=(x+0.5d0)*LOG(tmp)-tmp
    ser=1.000000000190015d0
    DO j=1,6
      y=y+1.d0
      ser=ser+cof(j)/y
    ENDDO ! j
    gammln=tmp+LOG(stp*ser/x)
  END FUNCTION gammln


END MODULE utils



MODULE rand_no_gen
  ! Pseudo-random number generator.
  USE utils,ONLY : dp,errstop
  IMPLICIT NONE
  PRIVATE
  PUBLIC ranx,rang
  INTEGER :: iseed=-1 ! Seed.  Supply a negative integer.


CONTAINS


  REAL(dp) FUNCTION ranx()
    ! Random number generator, adapted from ran2 in Numerical Recipes.
    ! (Method of l'Ecuyer with Bays-Durham shuffle.)
    IMPLICIT NONE
    INTEGER,PARAMETER :: im1=2147483563,im2=2147483399,imm1=im1-1,ia1=40014, &
      &ia2=40692,iq1=53668,iq2=52774,ir1=12211,ir2=3791,ntab=32, &
      &ndiv=1+imm1/ntab,ntabp8=ntab+8
    INTEGER :: j,k
    INTEGER,SAVE :: iseed2=123456789,iv(ntab)=0,iy=0
    REAL(dp),PARAMETER :: am=1.d0/im1,rnmx=1.d0-EPSILON(1.d0)
    IF(iseed<=0)THEN
      iseed=MAX(-iseed,1)
      iseed2=iseed
      DO j=ntabp8,1,-1
        k=iseed/iq1
        iseed=ia1*(iseed-k*iq1)-k*ir1
        IF(iseed<0)iseed=iseed+im1
        IF(j<=ntab)iv(j)=iseed
      ENDDO ! j
      iy=iv(1)
    ENDIF ! iseed<=0
    k=iseed/iq1
    iseed=ia1*(iseed-k*iq1)-k*ir1
    IF(iseed<0)iseed=iseed+im1
    k=iseed2/iq2
    iseed2=ia2*(iseed2-k*iq2)-k*ir2
    IF(iseed2<0)iseed2=iseed2+im2
    j=1+iy/ndiv
    iy=iv(j)-iseed2
    iv(j)=iseed
    IF(iy<1)iy=iy+imm1
    ranx=MIN(am*iy,rnmx)
  END FUNCTION ranx


  REAL(dp) FUNCTION rang(variance)
    ! Generate Gaussian-distributed random numbers.
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: variance
    REAL(dp) :: v1,v2,rad
    DO
      v1=2.d0*ranx()-1.d0
      v2=2.d0*ranx()-1.d0
      rad=v1*v1+v2*v2
      IF(rad<=1.d0.AND.rad>0.d0)EXIT
    ENDDO
    rang=v2*SQRT(-2.d0*variance*LOG(rad)/rad)
  END FUNCTION rang


END MODULE rand_no_gen


MODULE bootstrap
  ! Subroutines for Monte Carlo fitting.
  USE utils
  IMPLICIT NONE
  PRIVATE
  PUBLIC get_data,monte_carlo_fit_data,model_y,ndata,nparam,data_y_offset,&
    &rec_errbar_y,data_x
  ! Number of data points and number of parameters in model.
  INTEGER :: ndata,nparam
  ! Number of derived properties
  INTEGER :: nprop
  ! Raw data and error bars.
  REAL(dp),ALLOCATABLE :: data_x(:),data_y(:),rec_errbar_y(:),errbar_y_sq(:)
  ! Randomly offset data.
  REAL(dp),ALLOCATABLE :: data_y_offset(:)
  ! Initial guess at parameter vector.
  REAL(dp),ALLOCATABLE :: params_init(:)
  ! Names of the parameters.
  CHARACTER(2),ALLOCATABLE :: param_name(:)
  ! Names of the properties.
  CHARACTER(2),ALLOCATABLE :: prop_name(:)
  ! Arrays required by nl2sol.
  INTEGER,ALLOCATABLE :: iv(:)
  REAL(dp),ALLOCATABLE :: v(:)


CONTAINS


  SUBROUTINE get_file(in_file)
    ! Find out the name of the data file by one method or another.
    IMPLICIT NONE
    CHARACTER(*),INTENT(out) :: in_file
    INTEGER :: nargs,ierr
!!$    INTERFACE
!!$      SUBROUTINE getarg(i,c)
!!$        INTEGER,INTENT(in) :: i
!!$        CHARACTER(*),INTENT(out) :: c
!!$      END SUBROUTINE getarg
!!$      INTEGER FUNCTION iargc()
!!$      END FUNCTION iargc
!!$    END INTERFACE
    nargs=iargc()
    IF(nargs>0)THEN
      CALL getarg(1,in_file)
    ELSE
      DO
        WRITE(*,*)'Please enter name of data file.'
        READ(*,*,iostat=ierr)in_file
        IF(ierr==0)THEN
          EXIT
        ELSE
          WRITE(*,*)'Please try again.'
        ENDIF ! No error
      ENDDO
      WRITE(*,*)
    ENDIF ! Command-line argument supplied.
    in_file=ADJUSTL(in_file)
  END SUBROUTINE get_file


  SUBROUTINE get_data
    ! Count the lines in the input file, then read in the data.
    IMPLICIT NONE
    INTEGER :: ierr,i,ialloc
    INTEGER,PARAMETER :: io=8
    CHARACTER(200) :: char200
    CHARACTER(80) :: in_file
    REAL(dp) :: tempr(3)

    ! Find out file name.
    CALL get_file(in_file)
    WRITE(*,*)'Reading data from '//TRIM(in_file)//'.'

    OPEN(unit=io,file=TRIM(in_file),status='old',iostat=ierr)
    IF(ierr/=0)CALL errstop('GET_DATA','Problem opening file '//TRIM(in_file)&
      &//'.')

    ! Look for first line of data, to determine whether error bars are
    ! present.
    DO
      READ(io,'(a)',iostat=ierr)char200
      IF(ierr>0)THEN
        CALL errstop('GET_DATA','Error reading '//TRIM(in_file)//'.')
      ELSEIF(ierr<0)THEN
        CALL errstop('GET_DATA','File '//TRIM(in_file)&
          &//' doesn''t contain any data.')
      ENDIF ! ierr>0
      char200=ADJUSTL(char200)
      IF(isdataline(char200))THEN
        READ(char200,*,iostat=ierr)tempr(1:3)
        IF(ierr/=0)THEN
          WRITE(*,*)'File '//TRIM(in_file)//' has the wrong format.'
          WRITE(*,*)'Data should be in the form x y delta_y, where delta_y &
            &is the error in y.'
        ENDIF ! ierr/=0
        EXIT
      ENDIF ! Line begins with a number.
    ENDDO ! Lines
    ndata=1

    ! Count up the remaining lines of data in the file.
    DO
      READ(io,'(a)',iostat=ierr)char200
      IF(ierr>0)THEN
        CALL errstop('GET_DATA','Error reading '//TRIM(in_file)//'.')
      ELSEIF(ierr<0)THEN
        EXIT
      ENDIF ! ierr>0
      IF(isdataline(char200))ndata=ndata+1
    ENDDO ! Lines

    REWIND(io)

    WRITE(*,*)'Number of data lines: '//trim(i2s(ndata))

    ALLOCATE(data_y(ndata),data_y_offset(ndata),data_x(ndata),&
      &rec_errbar_y(ndata),errbar_y_sq(ndata),stat=ialloc)
    IF(ialloc/=0)CALL errstop('GET_DATA','Allocation error: data arrays.')

    i=0
    DO
      READ(io,'(a)',iostat=ierr)char200
      IF(ierr>0)THEN
        CALL errstop('GET_DATA','Error reading '//TRIM(in_file)//'.')
      ELSEIF(ierr<0)THEN
        EXIT
      ENDIF ! ierr>0
      char200=ADJUSTL(char200)
      IF(isdataline(char200))THEN
        i=i+1
        READ(char200,*,iostat=ierr)data_x(i),data_y(i),rec_errbar_y(i)
        IF(ierr/=0)THEN
          WRITE(*,*)'Problem reading '//TRIM(in_file)//'.'
          WRITE(*,*)'Problem line: '//TRIM(char200)
          CALL errstop('GET_DATA','Halting.')
        ENDIF ! ierr/=0
        IF(rec_errbar_y(i)<=0.d0)THEN
          WRITE(*,*)'Found a non-positive error bar.'
          WRITE(*,*)'Problem line: '//TRIM(char200)
          CALL errstop('GET_DATA','Halting.')
        ENDIF ! Error bar non-positive
        errbar_y_sq(i)=rec_errbar_y(i)**2
        rec_errbar_y(i)=1.d0/rec_errbar_y(i)
      ENDIF ! Line starts with a number.
    ENDDO ! Lines
    IF(i/=ndata)CALL errstop('GET_DATA','Bug.')

    CLOSE(io)

    WRITE(*,*)

    ! Make sure data are in ascending order of x.
    CALL sort_data

  END SUBROUTINE get_data


  SUBROUTINE sort_data
    ! Sort the x data into ascending order.  Warn if data not already in
    ! ascending order, or if there are repeated data items.
    IMPLICIT NONE
    INTEGER :: i,j
    LOGICAL :: warning_given_order=.FALSE.,warning_given_equal=.FALSE.
    DO i=1,ndata-1
      DO j=i+1,ndata
        IF(data_x(j)<data_x(i))THEN
          IF(.NOT.warning_given_order)THEN
            WRITE(*,*)'Warning: x data are not in ascending order.'
            WRITE(*,*)
            warning_given_order=.TRUE.
          ENDIF ! Haven't yet warned about order.
          CALL swapreal(data_x(i),data_x(j))
          CALL swapreal(data_y(i),data_y(j))
          CALL swapreal(rec_errbar_y(i),rec_errbar_y(j))
          CALL swapreal(errbar_y_sq(i),errbar_y_sq(j))
        ENDIF ! swap needed
        IF(.NOT.warning_given_equal.AND.data_x(j)==data_x(i))THEN
          WRITE(*,*)'Warning: two x data are equal.'
          WRITE(*,*)'Value: ',data_x(i)
          WRITE(*,*)
          warning_given_equal=.TRUE.
        ENDIF ! Give warning about two equal x data.
      ENDDO ! j
    ENDDO ! i
  END SUBROUTINE sort_data


  SUBROUTINE offset_data(randomise)
    ! Prepare the data actually to be used in the fit (which may be randomly
    ! offset).
    USE rand_no_gen, ONLY : rang
    IMPLICIT NONE
    LOGICAL,INTENT(in) :: randomise
    INTEGER :: i
    IF(randomise)THEN
      DO i=1,ndata
        data_y_offset(i)=data_y(i)+rang(errbar_y_sq(i))
      ENDDO ! i
    ELSE
      data_y_offset=data_y
    ENDIF ! randomise
  END SUBROUTINE offset_data


  SUBROUTINE fit_data(params)
    ! Perform the fit to the (offset) data using nl2sno.  Run it twice.
    USE toms573,ONLY : nl2sno,dfault
    IMPLICIT NONE
    REAL(dp),INTENT(inout) :: params(nparam)
    INTERFACE
      SUBROUTINE madr(n,p,params,nf,r,uiparm,urparm,ufparm)
        USE machine_constants
        INTEGER,INTENT(IN) :: n,p
        INTEGER,INTENT(INOUT) :: nf
        INTEGER,INTENT(INOUT),OPTIONAL :: uiparm(:)
        REAL(dp),INTENT(IN) :: params(:)
        REAL(dp),INTENT(INOUT),OPTIONAL :: urparm(:), ufparm
        REAL(dp),INTENT(OUT) :: r(:)
      END SUBROUTINE madr
    END INTERFACE
    CALL dfault(iv,v)
    iv(19:24)=0 ! Make nl2sol silent.
    CALL nl2sno(ndata,nparam,params,madr,iv,v)
    CALL dfault(iv,v)
    iv(19:24)=0 ! Make nl2sol silent.
    CALL nl2sno(ndata,nparam,params,madr,iv,v)
  END SUBROUTINE fit_data


  SUBROUTINE monte_carlo_fit_data
    ! Repeatedly choose offset data to be original data offset by a Gaussian-
    ! distributed random amount with variance given by the squared error bar
    ! on the data, then fit the parameters, and then compute the average and
    ! statistical error in the fitted parameters.
    IMPLICIT NONE
    INTEGER :: ialloc,i,j,ierr
    INTEGER,PARAMETER :: nrand_points=10000
    REAL(dp) :: rtemp
    REAL(dp),ALLOCATABLE :: params(:),params_opt(:),params2_opt(:), &
      &err_params_opt(:),props(:),props_opt(:),props2_opt(:),err_props_opt(:)
    LOGICAL,PARAMETER :: make_histogram=.TRUE.

    ! Get number of parameters and set initial parameter values.
    CALL initialise_model

    WRITE(*,*)'Initial parameter set:'
    CALL display_params(params_init)
    WRITE(*,*)

    ALLOCATE(iv(nparam+60),v(93+ndata*(nparam+3)+nparam*(3*nparam+33)/2),&
      &stat=ialloc)
    IF(ialloc/=0)CALL errstop('MONTE_CARLO_FIT_DATA',&
      &'Allocation error: nl2sol arrays.')
    ALLOCATE(params(nparam),params_opt(nparam),params2_opt(nparam), &
      &err_params_opt(nparam),props(nprop),props_opt(nprop),&
      &props2_opt(nprop),err_props_opt(nprop),stat=ialloc)
    IF(ialloc/=0)CALL errstop('MONTE_CARLO_FIT_DATA',&
      &'Allocation error: parameter arrays.')

    ! Perform an initial fit with non-offset data, to get actual
    ! parameter values and to generate starting guess for subsequent fits.
    ! Note that the averaged parameter values can give a large chi^2, so
    ! one should use the ones calculated here using the non-offset data.
    params=params_init
    CALL offset_data(.FALSE.) ! DON'T apply random offset to data.
    CALL fit_data(params)     ! Fit parameters to data.
    params_init=params        ! New starting guess.
    CALL derived_properties(params,props)
    WRITE(*,*)'Fitted parameter set:'
    CALL display_params(params_init,props=props)
    WRITE(*,*)

    ! Fit the data a large number of times and average the fitted parameter
    ! values.
    params_opt=0.d0 ; params2_opt=0.d0
    props_opt=0.d0 ; props2_opt=0.d0
    IF(make_histogram)THEN
      DO j=1,nprop
        OPEN(unit=8+j,status='replace',file='histogram_'&
          &//TRIM(ADJUSTL(prop_name(j)))//'.dat',iostat=ierr)
        IF(ierr/=0)CALL errstop('MONTE_CARLO_FIT_DATA',&
          &'Error opening histogram_'//TRIM(ADJUSTL(prop_name(j)))//'.dat.')
      ENDDO ! j
    ENDIF ! make_histogram
    DO i=1,nrand_points
      params=params_init       ! First guess at parameters.
      CALL offset_data(.TRUE.) ! Apply random offset to data.
      CALL fit_data(params)    ! Fit parameters to data.
      CALL derived_properties(params,props)
      params_opt=params_opt+params
      params2_opt=params2_opt+params**2
      props_opt=props_opt+props
      props2_opt=props2_opt+props**2
      IF(make_histogram)THEN
        DO j=1,nprop
          WRITE(8+j,*)props(j)
        ENDDO ! j
      ENDIF ! make_histogram
    ENDDO ! i
    rtemp=1.d0/DBLE(nrand_points)
    params_opt=params_opt*rtemp
    params2_opt=params2_opt*rtemp
    props_opt=props_opt*rtemp
    props2_opt=props2_opt*rtemp

    ! Work out the standard errors in the fitted parameters.
    ! This is the standard deviation in the fitted parameters.
    IF(nrand_points>1)THEN
      rtemp=DBLE(nrand_points)/DBLE(nrand_points-1)
      DO i=1,nparam
        err_params_opt(i)=SQRT(MAX(params2_opt(i)-params_opt(i)**2,0.d0)*rtemp)
      ENDDO ! i
      DO i=1,nprop
        err_props_opt(i)=SQRT(MAX(props2_opt(i)-props_opt(i)**2,0.d0)*rtemp)
      ENDDO ! i
    ELSE
      err_params_opt=0.d0
      err_props_opt=0.d0
    ENDIF ! nrand_points>1

    WRITE(*,*)'Monte Carlo fitted parameter set (using ' &
      &//TRIM(i2s(nrand_points))//' random points):'
    CALL display_params(params_opt,err_params_opt,props_opt,err_props_opt)
    WRITE(*,*)

    IF(make_histogram)THEN
      DO j=1,nprop
        CLOSE(8+j)
        WRITE(*,*)'Histogram of '//TRIM(ADJUSTL(prop_name(j))) &
          &//' written to histogram_'//TRIM(ADJUSTL(prop_name(j)))//'.dat.'
      ENDDO ! j
      WRITE(*,*)
    ENDIF ! make_histogram

    ! Deallocate arrays.
    DEALLOCATE(params,params_opt,params2_opt,err_params_opt,props, &
      &props_opt,props2_opt,err_props_opt)
    CALL finalise

  END SUBROUTINE monte_carlo_fit_data


  SUBROUTINE initialise_model
    ! In this subroutine, the number of parameters is specified, the 
    ! parameter list is allocated, and the initial values are assigned.
    ! Information about the model can also be written out here.
    ! The parameter names are also specified here.
    ! YOU NEED TO CHANGE THIS SUBROUTINE IF YOU HAVE A NEW MODEL.
    IMPLICIT NONE
    INTEGER :: ialloc
    nparam=4 ! Number of parameters in model.
    nprop=3  ! Number of parameter-dependent derived properties to report.
    WRITE(*,*)'Model used is: E(d) = E0 + E1.exp(-d/d0) + E2.log(d)'
    WRITE(*,*)'Derived properties: minimum dm, energy at min Em and frequency w.'
    WRITE(*,*)
    ALLOCATE(params_init(nparam),param_name(nparam),prop_name(nprop),stat=ialloc)
    IF(ialloc/=0)CALL errstop('INITIALISE_MODEL',&
      &'Problem allocating parameter vector.')
    param_name=(/'E0','E1','d0','E2'/) ! Parameter names.
    prop_name=(/'dm','Em','w '/)       ! Derived-property names.
    params_init(1)=-0.436759d0 ! Initial parameter values.
    params_init(2)=2.22176d0
    params_init(3)=0.581409d0
    params_init(4)=0.288348d0
  END SUBROUTINE initialise_model


  REAL(dp) FUNCTION model_y(params,x)
    ! The model function.
    ! YOU NEED TO CHANGE THIS SUBROUTINE IF YOU HAVE A NEW MODEL.
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: params(nparam),x
    model_y=params(1)+params(2)*EXP(-x/params(3))+params(4)*LOG(x)
  END FUNCTION model_y


  SUBROUTINE derived_properties(params,props)
    ! This function returns useful or interesting functions of the parameter
    ! values, which will have their mean and standard error evaluated.
    ! YOU NEED TO CHANGE THIS SUBROUTINE IF YOU HAVE A NEW MODEL.
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: params(nparam)
    REAL(dp),INTENT(out) :: props(nprop)
    REAL(dp) :: dm,dmp
    INTEGER :: i
    INTEGER,PARAMETER :: imin=10,imax=200
    REAL(dp),PARAMETER :: tol=1.d-15
    ! Find minimum of fitted potential using Newton-Raphson.
    dm=1.869d0
    i=0
    DO
      i=i+1
      dmp=dm
      dm=(2.d0*params(4)*params(3)**2*dm*exp(dm/params(3)) &
        &-params(2)*params(3)*dm**2-params(2)*dm**3) &
        &/(params(4)*params(3)**2*exp(dm/params(3))-params(2)*dm**2)
      IF(ABS(dm-dmp)<tol*ABS(dmp).AND.i>imin)EXIT
      IF(i>imax)THEN
        WRITE(*,*)'Warning: error converging Newton-Raphson iteration.'
        EXIT
      ENDIF ! i>imax
    ENDDO
    props(1)=dm ! Minimum of fitted potential.
    ! Value of fitted potential at minimum.
    props(2)=params(1)+params(2)*EXP(-dm/params(3))+params(4)*LOG(dm)
    ! Frequency omega (not including mass) for fitted potential.
    props(3)=sqrt(params(2)/params(3)**2*exp(-dm/params(3))-params(4)/dm**2)
  END SUBROUTINE derived_properties


  SUBROUTINE display_params(params,err_params,props,err_props)
    ! Write out the fitted parameters and the corresponding chi^2 fn.
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: params(nparam)
    REAL(dp),INTENT(in),OPTIONAL :: props(nprop),err_params(nparam),&
      &err_props(nprop)
    INTEGER :: i,j
    REAL(dp) :: chi2
    IF(PRESENT(err_params))THEN
      DO j=1,nparam
        IF(err_params(j)>0.d0)THEN
          WRITE(*,*)'Parameter '//param_name(j)//' = ',params(j),' +/- ', &
            &err_params(j)
          WRITE(*,*)'             = ' &
            &//TRIM(write_mean(params(j),err_params(j)))
        ELSE
          WRITE(*,*)'Parameter '//param_name(j)//' = ',params(j)
        ENDIF ! err_params
      ENDDO ! i
    ELSE
      DO j=1,nparam
        WRITE(*,*)'Parameter '//param_name(j)//' = ',params(j)
      ENDDO ! i
    ENDIF ! err_params present
    IF(nprop>0.AND.PRESENT(props))THEN
      IF(PRESENT(err_props))THEN
        DO j=1,nprop
          IF(err_props(j)>0.d0)THEN
            WRITE(*,*)'Property  '//prop_name(j)//' = ',props(j),' +/- ', &
              &err_props(j)
            WRITE(*,*)'             = ' &
              &//TRIM(write_mean(props(j),err_props(j)))
          ELSE
            WRITE(*,*)'Property  '//prop_name(j)//' = ',props(j)
          ENDIF ! err_props
        ENDDO ! i
      ELSE
        DO j=1,nprop
          WRITE(*,*)'Property  '//prop_name(j)//' = ',props(j)
        ENDDO ! i
      ENDIF ! err_props present
    ENDIF ! nprop>0
    chi2=0.d0
    DO i=1,ndata
      chi2=chi2+((data_y(i)-model_y(params,data_x(i)))*rec_errbar_y(i))**2
    ENDDO ! i
    WRITE(*,*)'chi^2 function per data point = ',chi2/DBLE(ndata)
    if(ndata>nparam)WRITE(*,*)'Probability chi^2 is this bad = ', &
      &gammq(DBLE(ndata-nparam)*0.5d0,chi2*0.5d0)
  END SUBROUTINE display_params


  SUBROUTINE finalise
    ! Deallocate arrays.
    IMPLICIT NONE
    IF(ALLOCATED(data_x))DEALLOCATE(data_x)
    IF(ALLOCATED(data_y))DEALLOCATE(data_y)
    IF(ALLOCATED(data_y_offset))DEALLOCATE(data_y_offset)
    IF(ALLOCATED(rec_errbar_y))DEALLOCATE(rec_errbar_y)
    IF(ALLOCATED(errbar_y_sq))DEALLOCATE(errbar_y_sq)
    IF(ALLOCATED(params_init))DEALLOCATE(params_init)
    IF(ALLOCATED(param_name))DEALLOCATE(param_name)
    IF(ALLOCATED(prop_name))DEALLOCATE(prop_name)
    IF(ALLOCATED(iv))DEALLOCATE(iv,v)
  END SUBROUTINE finalise


END MODULE bootstrap


PROGRAM fitter_mc
  ! Main program starts here.
  USE bootstrap,ONLY : get_data,monte_carlo_fit_data
  IMPLICIT NONE

  WRITE(*,*)
  WRITE(*,*)'FITTER_MC'
  WRITE(*,*)'========='
  WRITE(*,*)

  ! Read in the raw data.
  CALL get_data

  ! Carry out the Monte Carlo fit.
  CALL monte_carlo_fit_data

  WRITE(*,*)'Program finished.'
  WRITE(*,*)

END PROGRAM fitter_mc


SUBROUTINE madr(n,p,params,nf,r,uiparm,urparm,ufparm)
  ! Return the vector of residuals: the difference between the actual data 
  ! values and the model value for the given parameter set at each data point.
  USE utils,ONLY : dp,errstop
  USE bootstrap,ONLY : model_y,data_y_offset,rec_errbar_y,data_x,ndata,nparam
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: n,p
  INTEGER,INTENT(INOUT) :: nf
  INTEGER,INTENT(INOUT),OPTIONAL :: uiparm(:)
  REAL(dp),INTENT(IN) :: params(:)
  REAL(dp),INTENT(INOUT),OPTIONAL :: urparm(:), ufparm
  REAL(dp),INTENT(OUT) :: r(:)
  INTEGER :: i

  ! Check that p and n are number of params and data points, respectively.
  IF(p/=nparam)CALL errstop('MADR','p /= nparam in madr.')
  IF(n/=ndata)CALL errstop('MADR','n /= ndata in madr.')

  ! Evaluate the vector of residuals.
  DO i=1,ndata
    r(i)=(data_y_offset(i)-model_y(params,data_x(i)))*rec_errbar_y(i)
  ENDDO ! i

END SUBROUTINE madr
