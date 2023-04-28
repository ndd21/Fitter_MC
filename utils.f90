MODULE utils
  ! Miscellaneous utilities.
  USE iso_fortran_env,ONLY : uo=>output_unit,ue=>error_unit,ui=>input_unit
  USE machine_constants,ONLY : dp ! Double-precision type.
  IMPLICIT NONE
  INTEGER,PARAMETER :: luxlevel=3 ! Random "luxury" level.
  INTEGER,PARAMETER :: seed=31415 ! Random seed.
  INTERFACE
    SUBROUTINE ranlux(rvec,len)
      USE machine_constants,ONLY : dp
      REAL(dp),INTENT(out) :: rvec(*)
      INTEGER,INTENT(in) :: len
    END SUBROUTINE ranlux
    SUBROUTINE rluxgo(lux,int,k1,k2)
      INTEGER,INTENT(in) :: lux,int,k1,k2
    END SUBROUTINE rluxgo
  END INTERFACE
  PRIVATE
  PUBLIC uo,ui,errstop,wordwrap,i2s,write_mean,isdataline,gammq,initialise_rng,&
   &rang,sort_ranked


CONTAINS


  SUBROUTINE errstop(sub,message)
    ! Report an error and stop.
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: sub,message
    WRITE(ue,*)
    WRITE(ue,*)'ERROR in subroutine '//TRIM(ADJUSTL(sub))//'.'
    WRITE(ue,*)
    CALL wordwrap(TRIM(ADJUSTL(message)),ue)
    WRITE(ue,*)
    FLUSH(ue)
    STOP 1
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
      unit=uo
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
    DO i=1,HUGE(1)-1
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


  SUBROUTINE initialise_rng()
    ! Initialise the ranlux random-number generator.
    IMPLICIT NONE
    CALL rluxgo(luxlevel,seed,0,0)
  END SUBROUTINE initialise_rng


  SUBROUTINE rang(n,r)
    ! Generate a set of Gaussian-distributed random numbers using the
    ! Box-Muller algorithm.
    IMPLICIT NONE
    INTEGER,INTENT(in) :: n
    REAL(dp),INTENT(out) :: r(n)
    INTEGER :: i
    REAL(dp) :: v1,v2,f,rad,ranx(2)
    DO i=1,n,2
      DO
        CALL ranlux(ranx,2)
        v1=2.d0*ranx(1)-1.d0 ; v2=2.d0*ranx(2)-1.d0
        rad=v1*v1+v2*v2
        IF(rad<=1.d0.AND.rad>0.d0)EXIT
      ENDDO
      f=SQRT(-2.d0*LOG(rad)/rad)
      r(i)=v1*f
      IF(i<n)r(i+1)=v2*f
    ENDDO ! i
  END SUBROUTINE rang


  SUBROUTINE sort_ranked(ira,datara)
    ! Sort data in datara using a precomputed ranking array ira.
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ira(:)
    REAL(dp),INTENT(inout) :: datara(SIZE(ira))
    REAL(dp) :: tempdata(SIZE(ira))
    INTEGER :: i
    DO i=1,SIZE(ira)
      tempdata(i)=datara(ira(i))
    ENDDO ! i
    datara=tempdata
  END SUBROUTINE sort_ranked


END MODULE utils
