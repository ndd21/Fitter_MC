! FITTER_MC, Neil Drummond, 18/4/06.

! This is a general purpose utility for using nl2sno to fit a model to a set
! of data.  The name of the data file can be specified as a command-line
! argument,  Only lines in the data file that start with a number will
! be regarded as data lines; other lines will be ignored.  Error bars should be
! supplied, as the data will be Monte-Carlo simulated in order to estimate
! error bars on the fitted parameters.

! The program uses the iargc and getarg extensions to F90.  If these are not
! available then please uncomment the dummy versions below.

! SUBROUTINE getarg(i,c)
!  IMPLICIT NONE
!  INTEGER,INTENT(in) :: i
!  CHARACTER(*),INTENT(out) :: c
!  c=''
! END SUBROUTINE getarg

! INTEGER FUNCTION iargc()
!  IMPLICIT NONE
!  iargc=0
! END FUNCTION iargc


MODULE rand_no_gen
  ! Pseudo-random number generator.
  IMPLICIT NONE
  PRIVATE
  PUBLIC ranx,rang
  INTEGER :: iseed=-1 ! Seed.  Supply a negative integer.


CONTAINS


  DOUBLE PRECISION FUNCTION ranx()
    ! Random number generator, adapted from ran2 in Numerical Recipes.
    ! (Method of l'Ecuyer with Bays-Durham shuffle.)
    IMPLICIT NONE
    INTEGER,PARAMETER :: im1=2147483563,im2=2147483399,imm1=im1-1,ia1=40014, &
      &ia2=40692,iq1=53668,iq2=52774,ir1=12211,ir2=3791,ntab=32, &
      &ndiv=1+imm1/ntab,ntabp8=ntab+8
    INTEGER :: j,k
    INTEGER,SAVE :: iseed2=123456789,iv(ntab)=0,iy=0
    DOUBLE PRECISION,PARAMETER :: am=1.d0/im1,rnmx=1.d0-EPSILON(1.d0)
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


  DOUBLE PRECISION FUNCTION rang(variance)
    ! Generate Gaussian-distributed random numbers.
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(in) :: variance
    DOUBLE PRECISION :: v1,v2,rad
    DO
      v1=2.d0*ranx()-1.d0
      v2=2.d0*ranx()-1.d0
      rad=v1*v1+v2*v2
      IF(rad<=1.d0.AND.rad>0.d0)EXIT
    ENDDO
    rang=v2*SQRT(-2.d0*variance*LOG(rad)/rad)
  END FUNCTION rang


END MODULE rand_no_gen


MODULE utils
  IMPLICIT NONE
  PRIVATE
  PUBLIC get_data,monte_carlo_fit_data,model_y,no_data,no_params, &
    &data_y_offset,rec_errbar_y,data_x
  ! Number of data points and number of parameters in model.
  INTEGER :: no_data,no_params
  ! Raw data and error bars.
  DOUBLE PRECISION,ALLOCATABLE :: data_x(:),data_y(:),rec_errbar_y(:), &
    &errbar_y_sq(:)
  ! Randomly offset data.
  DOUBLE PRECISION,ALLOCATABLE :: data_y_offset(:)
  ! Initial guess at parameter vector.
  DOUBLE PRECISION,ALLOCATABLE :: params_init(:)
  ! Names of the parameters.
  CHARACTER(2),ALLOCATABLE :: param_name(:)
  ! Arrays required by nl2sol.
  INTEGER,ALLOCATABLE :: iv(:)
  DOUBLE PRECISION,ALLOCATABLE :: v(:)


CONTAINS


  SUBROUTINE get_file(in_file)
    ! Find out the name of the data file by one method or another.
    IMPLICIT NONE
    CHARACTER(*),INTENT(out) :: in_file
    INTEGER :: nargs,ierr
!    INTERFACE
!      SUBROUTINE getarg(i,c)
!        INTEGER,INTENT(in) :: i
!        CHARACTER(*),INTENT(out) :: c
!      END SUBROUTINE getarg
!      INTEGER FUNCTION iargc()
!      END FUNCTION iargc
!    END INTERFACE
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
    DOUBLE PRECISION :: tempr(3)

    ! Find out file name.
    CALL get_file(in_file)
    WRITE(*,*)'Reading data from '//TRIM(in_file)//'.'

    OPEN(unit=io,file=TRIM(in_file),status='old',iostat=ierr)
    IF(ierr/=0)THEN
      WRITE(*,*)'Problem opening file '//TRIM(in_file)//'.'
      STOP
    ENDIF ! ierr/=0

    ! Look for first line of data, to determine whether error bars are
    ! present.
    DO
      READ(io,'(a)',iostat=ierr)char200
      IF(ierr>0)THEN
        WRITE(*,*)'Error reading '//TRIM(in_file)//'.'
        STOP
      ELSEIF(ierr<0)THEN
        WRITE(*,*)'File '//TRIM(in_file)//' doesn''t contain any data.'
        STOP
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
    no_data=1

    ! Count up the remaining lines of data in the file.
    DO
      READ(io,'(a)',iostat=ierr)char200
      IF(ierr>0)THEN
        WRITE(*,*)'Error reading '//TRIM(in_file)//'.'
        STOP
      ELSEIF(ierr<0)THEN
        EXIT
      ENDIF ! ierr>0
      IF(isdataline(char200))no_data=no_data+1
    ENDDO ! Lines

    REWIND(io)

    WRITE(*,*)'Number of data lines: ',no_data

    ALLOCATE(data_y(no_data),data_y_offset(no_data),data_x(no_data), &
      &rec_errbar_y(no_data),errbar_y_sq(no_data),stat=ialloc)
    IF(ialloc/=0)THEN
      WRITE(*,*)'Allocation error: data arrays.'
      STOP
    ENDIF ! ialloc/=0

    i=0
    DO
      READ(io,'(a)',iostat=ierr)char200
      IF(ierr>0)THEN
        WRITE(*,*)'Error reading '//TRIM(in_file)//'.'
        STOP
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
          STOP
        ENDIF ! ierr/=0
        IF(rec_errbar_y(i)<=0.d0)THEN
          WRITE(*,*)'Found a non-positive error bar.'
          WRITE(*,*)'Problem line: '//TRIM(char200)
          STOP
        ENDIF ! Error bar non-positive
        errbar_y_sq(i)=rec_errbar_y(i)**2
        rec_errbar_y(i)=1.d0/rec_errbar_y(i)
      ENDIF ! Line starts with a number.
    ENDDO ! Lines
    IF(i/=no_data)THEN
      WRITE(*,*)'Bug.'
      STOP
    ENDIF ! i/=no_data

    CLOSE(io)

    WRITE(*,*)

    ! Make sure data are in ascending order of x.
    CALL sort_data

  END SUBROUTINE get_data


  LOGICAL FUNCTION isdataline(c)
    ! This function returns T if & only if the first word of c is a number.
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: c
    INTEGER :: ierr
    DOUBLE PRECISION :: r
    READ(c,*,iostat=ierr)r
    isdataline=(ierr==0)
  END FUNCTION isdataline


  SUBROUTINE sort_data
    ! Sort the x data into ascending order.  Warn if data not already in
    ! ascending order, or if there are repeated data items.
    IMPLICIT NONE
    INTEGER :: i,j
    LOGICAL :: warning_given_order=.FALSE.,warning_given_equal=.FALSE.
    DO i=1,no_data-1
      DO j=i+1,no_data
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


  SUBROUTINE swapreal(a,b)
    ! This subroutine swaps two real numbers a and b.
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(inout) :: a,b
    DOUBLE PRECISION :: temp
    temp=a ; a=b ; b=temp
  END SUBROUTINE swapreal


  SUBROUTINE offset_data(randomise)
    ! Prepare the data that is actually to be used in the fit (which may
    ! be randomly offset).
    USE rand_no_gen, ONLY : rang
    IMPLICIT NONE
    LOGICAL,INTENT(in) :: randomise
    INTEGER :: i
    IF(randomise)THEN
      DO i=1,no_data
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
    DOUBLE PRECISION,INTENT(inout) :: params(no_params)
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
    CALL nl2sno(no_data,no_params,params,madr,iv,v)
    CALL dfault(iv,v)
    iv(19:24)=0 ! Make nl2sol silent.
    CALL nl2sno(no_data,no_params,params,madr,iv,v)
  END SUBROUTINE fit_data


  SUBROUTINE monte_carlo_fit_data
    ! Repeatedly choose offset data to be original data offset by a Gaussian-
    ! distributed random amount with variance given by the squared error bar
    ! on the data, then fit the parameters, and then compute the average and
    ! statistical error in the fitted parameters.
    IMPLICIT NONE
    INTEGER :: ialloc,i
    INTEGER,PARAMETER :: no_rand_points=10000
    DOUBLE PRECISION :: rtemp
    DOUBLE PRECISION,ALLOCATABLE :: params(:),params_opt(:),params2_opt(:), &
      &err_params_opt(:)

    ! Get number of parameters and set initial parameter values.
    CALL initialise_model

    WRITE(*,*)'Initial parameter set:'
    CALL display_params(params_init)
    WRITE(*,*)

    ALLOCATE(iv(no_params+60),v(93+no_data*(no_params+3)+no_params &
      &*(3*no_params+33)/2),stat=ialloc)
    IF(ialloc/=0)THEN
      WRITE(*,*)'Allocation error: nl2sol arrays.'
      STOP
    ENDIF ! ialloc/=0
    ALLOCATE(params(no_params),params_opt(no_params),params2_opt(no_params), &
      &err_params_opt(no_params),stat=ialloc)
    IF(ialloc/=0)THEN
      WRITE(*,*)'Allocation error: parameter arrays.'
      STOP
    ENDIF ! ialloc/=0

    ! Perform an initial fit with non-offset data, to get actual
    ! parameter values and to generate starting guess for subsequent fits.
    ! Note that the averaged parameter values can give a large chi^2, so
    ! one should use the ones calculated here using the non-offset data.
    params=params_init
    CALL offset_data(.FALSE.) ! DON'T apply random offset to data.
    CALL fit_data(params)     ! Fit parameters to data.
    params_init=params        ! New starting guess.
    WRITE(*,*)'Fitted parameter set:'
    CALL display_params(params_init)
    WRITE(*,*)

    ! Fit the data a large number of times and average the fitted parameter
    ! values.
    params_opt=0.d0 ; params2_opt=0.d0
    DO i=1,no_rand_points
      params=params_init       ! First guess at parameters.
      CALL offset_data(.TRUE.) ! Apply random offset to data.
      CALL fit_data(params)    ! Fit parameters to data.
      params_opt=params_opt+params
      params2_opt=params2_opt+params**2
    ENDDO ! i
    rtemp=1.d0/DBLE(no_rand_points)
    params_opt=params_opt*rtemp
    params2_opt=params2_opt*rtemp

    ! Work out the standard errors in the fitted parameters.
    ! This is the standard deviation in the fitted parameters.
    IF(no_rand_points>1)THEN
      rtemp=DBLE(no_rand_points)/DBLE(no_rand_points-1)
      DO i=1,no_params
        err_params_opt(i)=SQRT(MAX(params2_opt(i)-params_opt(i)**2,0.d0) &
          &*rtemp)
      ENDDO ! i
    ELSE
      err_params_opt=0.d0
    ENDIF ! no_rand_points>1

    WRITE(*,*)'Monte Carlo fitted parameter set (using ' &
      &//TRIM(i2s(no_rand_points))//' random points):'
    CALL display_params(params_opt,err_params_opt)
    WRITE(*,*)

    ! Deallocate arrays.
    DEALLOCATE(params,params_opt,params2_opt,err_params_opt)
    CALL finalise

  END SUBROUTINE monte_carlo_fit_data


  SUBROUTINE initialise_model
    ! In this subroutine, the number of parameters is specified, the 
    ! parameter list is allocated, and the initial values are assigned.
    ! Information about the model can also be written out here.
    ! The parameter names are also specified here.
    IMPLICIT NONE
    INTEGER :: ialloc
    no_params=2
    !    WRITE(*,*)'Model used is: rho(p) = A + B.p^2 + C.p^4'
    !    WRITE(*,*)
    ALLOCATE(params_init(no_params),param_name(no_params),stat=ialloc)
    IF(ialloc/=0)THEN
      WRITE(*,*)'Problem allocating parameter vector.'
      STOP
    ENDIF ! ialloc/=0
    param_name=(/'A ','B '/)
    params_init(1)=1.d0
    params_init(2)=2.4d0
  END SUBROUTINE initialise_model


  DOUBLE PRECISION FUNCTION model_y(params,x)
    ! The model function.
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(in) :: params(no_params),x
    model_y=params(1)*x**params(2)
  END FUNCTION model_y


  SUBROUTINE display_params(params,err_params)
    ! Write out the fitted parameters and the corresponding chi^2 fn.
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(in) :: params(no_params)
    DOUBLE PRECISION,INTENT(in),OPTIONAL :: err_params(no_params)
    INTEGER :: i,j
    DOUBLE PRECISION :: chi2
    IF(PRESENT(err_params))THEN
      DO j=1,no_params
        IF(err_params(j)>0.d0)THEN
          WRITE(*,*)'Parameter '//param_name(j)//' = ',params(j),' +/- ', &
            &err_params(j)
          WRITE(*,*)'            = ' &
            &//TRIM(write_mean(params(j),err_params(j)))
        ELSE
          WRITE(*,*)'Parameter '//param_name(j)//' = ',params(j)
        ENDIF ! err_params
      ENDDO ! i
    ELSE
      DO j=1,no_params
        WRITE(*,*)'Parameter '//param_name(j)//' = ',params(j)
      ENDDO ! i
    ENDIF ! err_params present
    chi2=0.d0
    DO i=1,no_data
      chi2=chi2+((data_y(i)-model_y(params,data_x(i)))*rec_errbar_y(i))**2
    ENDDO ! i
    WRITE(*,*)'chi^2 function per data point = ',chi2/DBLE(no_data)
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
    IF(ALLOCATED(iv))DEALLOCATE(iv,v)
  END SUBROUTINE finalise


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
    DOUBLE PRECISION,INTENT(in) :: av,std_err_in_mean
    INTEGER,INTENT(in),OPTIONAL :: err_prec_in
    INTEGER :: lowest_digit_to_quote,err_quote,err_prec,int_part,dec_part,i
    INTEGER,PARAMETER :: err_prec_default=1
    DOUBLE PRECISION :: av_quote
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


END MODULE utils


PROGRAM fitter_mc
  ! Main program starts here.
  USE utils,ONLY : get_data,monte_carlo_fit_data
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
  USE machine_constants, ONLY : dp
  USE utils,ONLY : model_y,data_y_offset,rec_errbar_y,data_x,no_data,no_params
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: n,p
  INTEGER,INTENT(INOUT) :: nf
  INTEGER,INTENT(INOUT),OPTIONAL :: uiparm(:)
  REAL(dp),INTENT(IN) :: params(:)
  REAL(dp),INTENT(INOUT),OPTIONAL :: urparm(:), ufparm
  REAL(dp),INTENT(OUT) :: r(:)
  INTEGER :: i

  ! Check that p and n are number of params and data points, respectively.
  IF(p/=no_params)THEN
    WRITE(*,*)'p /= no_params in madr.'
    STOP
  ENDIF ! Parameter vector wrong length.
  IF(n/=no_data)THEN
    WRITE(*,*)'n /= no_data in madr.'
    STOP
  ENDIF ! Residual vector wrong length.

  ! Evaluate the vector of residuals.
  DO i=1,no_data
    r(i)=(data_y_offset(i)-model_y(params,data_x(i)))*rec_errbar_y(i)
  ENDDO ! i

END SUBROUTINE madr
