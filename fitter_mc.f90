PROGRAM fitter_mc
  ! FITTER_MC, Neil Drummond, 18/4/06.

  ! This is a general-purpose utility for using nl2sno to fit a model to a set
  ! of data.  The name of the data file can be specified as a command-line
  ! argument.  Error bars should be supplied, as the data will be Monte Carlo
  ! simulated in order to estimate error bars on the fitted parameters and
  ! properties expressed in terms of the parameters.

  ! If you wish to change the model, you only need to alter the three
  ! subroutines in bootstrap.f90 containing the comment

  ! YOU NEED TO CHANGE THIS SUBROUTINE IF YOU HAVE A NEW MODEL.

  USE bootstrap,ONLY : get_data,monte_carlo_fit_data
  !$ USE omp_lib,ONLY : omp_set_dynamic,omp_set_nested,omp_get_num_threads
  USE utils,ONLY : initialise_rng
  IMPLICIT NONE
  !$ INTEGER :: nthreads

  WRITE(*,*)
  WRITE(*,*)'FITTER_MC'
  WRITE(*,*)'========='
  !$ CALL omp_set_dynamic(.false.)
  !$ CALL omp_set_nested(.false.)
  !$OMP PARALLEL DEFAULT(none) SHARED(nthreads)
  !$ nthreads=omp_get_num_threads()
  !$OMP END PARALLEL
  !$ WRITE(*,'(1x,a,i0,a)')'Running in parallel using ',nthreads,' threads.'
  WRITE(*,*)

  CALL initialise_rng

  CALL get_data

  CALL monte_carlo_fit_data

  WRITE(*,*)'Program finished.'
  WRITE(*,*)

END PROGRAM fitter_mc
