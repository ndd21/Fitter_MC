MODULE bootstrap
  ! Subroutines for Monte Carlo fitting.  See comments at top of fitter_mc.f90.
  USE machine_constants,ONLY : dp
  USE utils,ONLY : errstop,wordwrap,i2s,write_mean,isdataline,gammq,&
    &rang,sort_ranked
  IMPLICIT NONE
  PRIVATE
  PUBLIC get_data,monte_carlo_fit_data
  INTEGER,PARAMETER :: nrand_points=100000 ! No. of samples for bootstrap MC.
  ! Number of data points and number of parameters in model.
  INTEGER :: ndata,nparam
  INTEGER :: nprop ! Number of derived properties
  ! Raw data and error bars.
  REAL(dp),ALLOCATABLE :: data_x(:),data_y(:),rec_errbar_y(:),errbar_y(:)
  REAL(dp),ALLOCATABLE :: data_y_offset(:) ! Randomly offset data.
  !$OMP THREADPRIVATE(data_y_offset)
  REAL(dp),ALLOCATABLE :: params_init(:) ! Initial guess at parameter vector.
  CHARACTER(2),ALLOCATABLE :: param_name(:) ! Names of the parameters.
  CHARACTER(2),ALLOCATABLE :: prop_name(:) ! Names of the properties.
  ! Arrays required by nl2sol.
  INTEGER,ALLOCATABLE :: iv(:)
  !$OMP THREADPRIVATE(iv)
  REAL(dp),ALLOCATABLE :: v(:)
  !$OMP THREADPRIVATE(v)


CONTAINS


  SUBROUTINE get_file(in_file)
    ! Find out the name of the data file by one method or another.
    IMPLICIT NONE
    CHARACTER(*),INTENT(out) :: in_file
    INTEGER :: nargs,ierr
    nargs=command_argument_COUNT()
    IF(nargs>0)THEN
      CALL get_command_ARGUMENT(1,in_file)
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
    CHARACTER(200) :: in_file
    REAL(dp) :: tempr(3)

    ! Find out file name.
    CALL get_file(in_file)
    CALL wordwrap('Reading data from '//TRIM(in_file)//'.')

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
          CALL wordwrap('File '//TRIM(in_file)//' has the wrong format.')
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

    WRITE(*,*)'Number of data lines: '//TRIM(i2s(ndata))

    ALLOCATE(data_y(ndata),data_x(ndata),rec_errbar_y(ndata),&
      &errbar_y(ndata),stat=ialloc)
    IF(ialloc/=0)CALL errstop('GET_DATA','Allocation error: data arrays.')
    !$OMP PARALLEL DEFAULT(none) PRIVATE(ialloc) SHARED(ndata)
    ALLOCATE(data_y_offset(ndata),stat=ialloc)
    IF(ialloc/=0)CALL errstop('GET_DATA','Allocation error: data_y_offset.')
    !$OMP END PARALLEL

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
        errbar_y(i)=rec_errbar_y(i)
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
    USE m_mrgrnk,ONLY : mrgrnk
    IMPLICIT NONE
    INTEGER :: i,ira(ndata)
    LOGICAL :: out_of_order,warning_given_order=.FALSE.,&
      &warning_given_equal=.FALSE.
    REAL(dp),PARAMETER :: tol=1.d-13
    CALL mrgrnk(data_x,ira)
    out_of_order=.FALSE.
    DO i=1,ndata
      IF(ira(i)/=i)THEN
        out_of_order=.TRUE.
        EXIT
      ENDIF ! ira/=i
    ENDDO ! i
    IF(out_of_order)THEN
      IF(.NOT.warning_given_order)THEN
        WRITE(*,*)'Warning: raw x data are not in ascending order.'
        WRITE(*,*)
        warning_given_order=.TRUE.
      ENDIF ! not warning_given_order
      CALL sort_ranked(ira,data_x)
      CALL sort_ranked(ira,data_y)
      CALL sort_ranked(ira,rec_errbar_y)
      CALL sort_ranked(ira,errbar_y)
    ENDIF ! out_of_order
    IF(.NOT.warning_given_equal)THEN
      DO i=2,ndata
        IF(ABS(data_x(i)-data_x(i-1))<tol)THEN
          WRITE(*,*)'Warning: two x data are equal.'
          WRITE(*,*)'Value: ',data_x(i)
          WRITE(*,*)
          warning_given_equal=.TRUE.
          EXIT
        ENDIF ! Give warning about two equal x data.
      ENDDO ! i
    ENDIF ! not warning_equal
  END SUBROUTINE sort_data


  SUBROUTINE offset_data(randomise)
    ! Prepare the data actually to be used in the fit (which may be randomly
    ! offset).
    USE utils,ONLY : rang
    IMPLICIT NONE
    LOGICAL,INTENT(in) :: randomise
    REAL(dp) :: offset(ndata)
    IF(randomise)THEN
      !$OMP CRITICAL
      call rang(ndata,offset)
      !$OMP END CRITICAL
      data_y_offset=data_y+offset*errbar_y
    ELSE
      data_y_offset=data_y
    ENDIF ! randomise
  END SUBROUTINE offset_data


  SUBROUTINE fit_data(params)
    ! Perform the fit to the (offset) data using nl2sno.  Run it twice.
    USE toms573,ONLY : nl2sno,dfault
    IMPLICIT NONE
    REAL(dp),INTENT(inout) :: params(nparam)
    CALL dfault(iv,v)
    iv(19:24)=0 ! Make nl2sol silent.
    CALL nl2sno(ndata,nparam,params,madr,iv,v)
    CALL dfault(iv,v)
    iv(19:24)=0 ! Make nl2sol silent.
    CALL nl2sno(ndata,nparam,params,madr,iv,v)
  END SUBROUTINE fit_data


  SUBROUTINE madr(n,p,params,nf,r,uiparm,urparm,ufparm)
    ! Return the vector of residuals: the difference between the actual data 
    ! values and the model value for the given parameter set at each data point.
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


  SUBROUTINE monte_carlo_fit_data
    ! Repeatedly choose offset data to be original data offset by a Gaussian-
    ! distributed random amount with variance given by the squared error bar
    ! on the data, then fit the parameters, and then compute the average and
    ! statistical error in the fitted parameters.
    IMPLICIT NONE
    INTEGER :: ialloc,i,j,ierr
    REAL(dp) :: rtemp
    REAL(dp),ALLOCATABLE :: params(:),params_opt(:),params2_opt(:), &
      &err_params_opt(:),props(:),props_opt(:),props2_opt(:),err_props_opt(:)
    LOGICAL,PARAMETER :: make_histogram=.TRUE.

    ! Get number of parameters and set initial parameter values.
    CALL initialise_model

    WRITE(*,*)'Initial parameter set:'
    CALL display_params(params_init)
    WRITE(*,*)

    !$OMP PARALLEL DEFAULT(none) PRIVATE(ialloc) SHARED(nparam,ndata)
    ALLOCATE(iv(nparam+60),v(93+ndata*(nparam+3)+nparam*(3*nparam+33)/2),&
      &stat=ialloc)
    IF(ialloc/=0)CALL errstop('MONTE_CARLO_FIT_DATA',&
      &'Allocation error: nl2sol arrays.')
    !$OMP END PARALLEL
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
    !$OMP PARALLEL DO DEFAULT(none) PRIVATE(i,params,props,j) &
    !$OMP   &SHARED(params_init,nprop) REDUCTION(+:params_opt,params2_opt,&
    !$OMP   &props_opt,props2_opt)
    DO i=1,nrand_points
      params=params_init       ! First guess at parameters.
      CALL offset_data(.TRUE.) ! Apply random offset to data.
      CALL fit_data(params)    ! Fit parameters to data.
      CALL derived_properties(params,props)
      params_opt=params_opt+params ; params2_opt=params2_opt+params**2
      props_opt=props_opt+props ; props2_opt=props2_opt+props**2
      IF(make_histogram)THEN
        !$OMP CRITICAL
        DO j=1,nprop
          WRITE(8+j,*)props(j)
        ENDDO ! j
        !$OMP END CRITICAL
      ENDIF ! make_histogram
    ENDDO ! i
    !$OMP END PARALLEL DO
    rtemp=1.d0/DBLE(nrand_points)
    params_opt=params_opt*rtemp ; params2_opt=params2_opt*rtemp
    props_opt=props_opt*rtemp ; props2_opt=props2_opt*rtemp

    ! Work out the standard errors in the fitted parameters.
    ! This is the standard deviation in the fitted parameters.
    IF(nrand_points>1)THEN
      rtemp=DBLE(nrand_points)/DBLE(nrand_points-1)
      err_params_opt=SQRT(MAX(params2_opt-params_opt**2*rtemp,0.d0))
      err_props_opt=SQRT(MAX(props2_opt-props_opt**2*rtemp,0.d0))
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
        WRITE(*,*)'  Histogram of '//TRIM(ADJUSTL(prop_name(j))) &
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
    nparam=3 ! Number of parameters in model.
    nprop=1  ! Number of parameter-dependent derived properties to report.
    WRITE(*,*)'Model used is: test'
    WRITE(*,'(1x,a,i0)')'Number of parameters: ',nparam
    WRITE(*,*)
    ALLOCATE(params_init(nparam),param_name(nparam),prop_name(nprop),&
      &stat=ialloc)
    IF(ialloc/=0)CALL errstop('INITIALISE_MODEL',&
      &'Problem allocating parameter vector.')
    param_name=(/'a ','b ','al'/) ! Parameter names.
    prop_name=(/'ab'/)       ! Derived-property names.
    params_init(1)=0.d0 ! Initial parameter values.
    params_init(2)=2.d0
    params_init(3)=1.d0
  END SUBROUTINE initialise_model


  REAL(dp) FUNCTION model_y(params,x)
    ! The model function.
    ! YOU NEED TO CHANGE THIS SUBROUTINE IF YOU HAVE A NEW MODEL.
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: params(nparam),x
    model_y=params(1)+params(2)*EXP(-params(3)*x)
  END FUNCTION model_y


  SUBROUTINE derived_properties(params,props)
    ! This function returns useful or interesting functions of the parameter
    ! values, which will have their mean and standard error evaluated.
    ! YOU NEED TO CHANGE THIS SUBROUTINE IF YOU HAVE A NEW MODEL.
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: params(nparam)
    REAL(dp),INTENT(out) :: props(nprop)
    props(1)=params(1)+params(2)
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
          WRITE(*,*)'  Parameter '//param_name(j)//' = ',params(j),' +/- ', &
            &err_params(j)
          WRITE(*,*)'               = ' &
            &//TRIM(write_mean(params(j),err_params(j)))
        ELSE
          WRITE(*,*)'Parameter '//param_name(j)//' = ',params(j)
        ENDIF ! err_params
      ENDDO ! j
    ELSE
      DO j=1,nparam
        WRITE(*,*)'  Parameter '//param_name(j)//' = ',params(j)
      ENDDO ! j
    ENDIF ! err_params present
    IF(nprop>0.AND.PRESENT(props))THEN
      IF(PRESENT(err_props))THEN
        DO j=1,nprop
          IF(err_props(j)>0.d0)THEN
            WRITE(*,*)'  Property  '//prop_name(j)//' = ',props(j),' +/- ', &
              &err_props(j)
            WRITE(*,*)'               = ' &
              &//TRIM(write_mean(props(j),err_props(j)))
          ELSE
            WRITE(*,*)'Property  '//prop_name(j)//' = ',props(j)
          ENDIF ! err_props
        ENDDO ! j
      ELSE
        DO j=1,nprop
          WRITE(*,*)'  Property  '//prop_name(j)//' = ',props(j)
        ENDDO ! j
      ENDIF ! err_props present
    ENDIF ! nprop>0
    chi2=0.d0
    !$OMP PARALLEL DO DEFAULT(none) PRIVATE(i) SHARED(ndata,data_x,data_y,&
    !$OMP   &params,rec_errbar_y) REDUCTION(+:chi2)
    DO i=1,ndata
      chi2=chi2+((data_y(i)-model_y(params,data_x(i)))*rec_errbar_y(i))**2
    ENDDO ! i
    !$OMP END PARALLEL DO
    WRITE(*,*)'  Reduced chi^2 function = ',chi2/DBLE(ndata-nparam)
    IF(ndata>nparam)WRITE(*,*)'  Probability chi^2 is this bad = ', &
      &gammq(DBLE(ndata-nparam)*0.5d0,chi2*0.5d0)
  END SUBROUTINE display_params


  SUBROUTINE finalise
    ! Deallocate arrays.
    IMPLICIT NONE
    IF(ALLOCATED(data_x))DEALLOCATE(data_x)
    IF(ALLOCATED(data_y))DEALLOCATE(data_y)
    !$OMP PARALLEL DEFAULT(none)
    IF(ALLOCATED(data_y_offset))DEALLOCATE(data_y_offset)
    IF(ALLOCATED(iv))DEALLOCATE(iv,v)
    !$OMP END PARALLEL
    IF(ALLOCATED(rec_errbar_y))DEALLOCATE(rec_errbar_y)
    IF(ALLOCATED(errbar_y))DEALLOCATE(errbar_y)
    IF(ALLOCATED(params_init))DEALLOCATE(params_init)
    IF(ALLOCATED(param_name))DEALLOCATE(param_name)
    IF(ALLOCATED(prop_name))DEALLOCATE(prop_name)
  END SUBROUTINE finalise


END MODULE bootstrap
