MODULE toms573
!-------------------------------------------------------------------!
! NL2SOL -- an adaptive nonlinear least-squares algorithm           !
! ------                                                            !
!                                                                   !
! Authors : John E. Dennis, Jr., David M. Gay, and Roy E. Welsch    !
!                                                                   !
! ACM Transactions on Mathematical Software, September, 1981.       !
!                                                                   !
! Translation from Fortran 66/77 format to compatibility with ELF90 !
! by Alan Miller, 5 June 1997 (this revision 10th July 1997)        !
!                                                                   !
! This version shoehorned into CASINO 10.2002 (MDT/NDD)             !
!                                                                   !
! NOTE : apparently pgf90 refused to compile this (see the note     !
!        at the end of the file.). Modifications were needed.       !
!                                                                   !
! Alterations:                                                      !
! 10.2002 NDD  Set dp to give default double precision at runtime.  !
!  4.2003 NDD  Pass iteration number to least-squares function.     !
!              Change made to both nl2sno and nl2sol.               !
!  7.2004 PLR  Report mode (evaluation/derivative) and status of    !
!              previous correction (success/failure) to evalr.      !
!  7.2004 PLR  Array indices made explicit in calls to subroutines  !
!              so that pgf compiles the code properly.              !
! 12.2004 PLR  Beautification of output.                            !
!-------------------------------------------------------------------!
 USE machine_constants
 IMPLICIT NONE
 LOGICAL :: stop_nl2sol=.false.

CONTAINS


 SUBROUTINE nl2sol(n,p,x,calcr,calcj,iv,v,uiparm,urparm,ufparm)
!-----------------------------------------------------------------------!
! Minimize nonlinear sum of squares using analytic Jacobian             !
! (NL2SOL version 2.2)                                                  !
!-----------------------------------------------------------------------!
 INTEGER,INTENT(in) :: n, p
 INTEGER,INTENT(inout) :: iv(:)
 INTEGER,INTENT(inout),OPTIONAL :: uiparm(:)
 REAL(dp),INTENT(inout) :: x(:),v(:)
 REAL(dp),INTENT(inout),OPTIONAL :: urparm(:),ufparm

 INTERFACE
  SUBROUTINE calcr(n,p,x,nf,r,uiparm,urparm,ufparm)
   USE machine_constants
   IMPLICIT NONE
   INTEGER,INTENT(in) :: n,p
   INTEGER,INTENT(inout) :: nf
   INTEGER,INTENT(inout),OPTIONAL :: uiparm(:)
   REAL(dp),INTENT(in) :: x(:)
   REAL(dp),INTENT(inout),OPTIONAL :: urparm(:),ufparm
   REAL(dp),INTENT(out) :: r(:)
  END SUBROUTINE calcr
  SUBROUTINE calcj(n,p,x,nf,j,uiparm,urparm,ufparm)
   USE machine_constants
   IMPLICIT NONE
   INTEGER,INTENT(in) :: n, p
   INTEGER,INTENT(inout) :: nf
   INTEGER,INTENT(inout),OPTIONAL :: uiparm(:)
   REAL(dp),INTENT(in) :: x(:)
   REAL(dp),INTENT(inout),OPTIONAL :: urparm(:),ufparm
   REAL(dp),INTENT(out) :: j(:)
  END SUBROUTINE calcj
 END INTERFACE

!------------------------------------------------------------------------!
! Given a p-vector x of parameters, calcr computes an n-vector           !
! r = r(x) of residuals corresponding to x.  (r(x) probably arises       !
! from a nonlinear model involving p parameters and n observations.)     !
! this routine interacts with nl2itr to seek a parameter vector x        !
! that minimizes the sum of the squares of (the components of) r(x),     !
! i.e., that minimizes the sum-of-squares function                       !
! f(x) = (r(x)**t) * r(x) / 2.  r(x) is assumed to be a twice            !
! continuously differentiable function of x.                             !
!------------------------------------------------------------------------!
!
!--------------------------  parameter usage  --------------------------
!
! n........ (input) the number of observations, i.e., the number of
!                  components in r(x).  n must be >= p.
! p........ (input) the number of parameters (components in x).  p must
!                  be positive.
! x........ (input/output).  on input, x is an initial guess at the
!                  desired parameter estimate.  on output, x contains
!                  the best parameter estimate found.
! calcr.... (input) a subroutine which, given x, computes r(x).  calcr
!                  must be declared external in the calling program.
!                  it is invoked by
!                       call calcr(n,p,x,nf,r,uiparm,urparm,ufparm)
!                  when calcr is called, nf is the invocation count
!                  for calcr.  it is included for possible use with
!                  calcj.  if x is out of bounds (e.g. if it would
!                  cause overflow in computing r(x)), then calcr should
!                  set nf to 0.  this will cause a shorter step to be
!                  attempted.  the other parameters are as described
!                  above and below.  calcr should not change n, p, or x.
! calcj.... (input) a subroutine which, given x, computes the jacobian
!                  matrix j of r at x, i.e., the n by p matrix whose
!                  (i,k) entry is the partial derivative of the i-th
!                  component of r with respect to x(k).  calcj must be
!                  declared external in the calling program.  it is
!                  invoked by
!                       call calcj(n,p,x,nf,j,uiparm,urparm,ufparm)
!                  nf is the invocation count for calcr at the time
!                  r(x) was evaluated.  the x passed to calcj is
!                  usually the one passed to calcr on either its most
!                  recent invocation or the one prior to it.  if calcr
!                  saves intermediate results for use by calcj, then it
!                  is possible to tell from nf whether they are valid
!                  for the current x (or which copy is valid if two
!                  copies are kept).  if j cannot be computed at x,
!                  then calcj should set nf to 0.  in this case, nl2sol
!                  will return with iv(1) = 15.  the other parameters
!                  to calcj are as described above and below.  calcj
!                  should not change n, p, or x.
! iv....... (input/output) an integer value array of length at least
!                  60 + p that helps control the nl2sol algorithm and
!                  that is used to store various intermediate quanti-
!                  ties.  of particular interest are the initialization/
!                  return code iv(1) and the entries in iv that control
!                  printing and limit the number of iterations and func-
!                  tion evaluations.  see the section on iv input
!                  values below.
! v........ (input/output) a floating-point value array of length at
!                  least 93 + n*p + 3*n + p*(3*p+33)/2 that helps con-
!                  trol the nl2sol algorithm and that is used to store
!                  various intermediate quantities.  of particular in-
!                  terest are the entries in v that limit the length of
!                  the first step attempted (lmax0), specify conver-
!                  gence tolerances (afctol, rfctol, xctol, xftol),
!                  and help choose the step size used in computing the
!                  covariance matrix (delta0).  see the section on
!                  (selected) v input values below.
! uiparm... (input) user integer parameter array passed without change
!                  to calcr and calcj.
! urparm... (input) user floating-point parameter array passed without
!                  change to calcr and calcj.
! ufparm... (input) user external subroutine or function passed without
!                  change to calcr and calcj.
!
!  ***  iv input values (from subroutine dfault)  ***
!
! iv(1)...  on input, iv(1) should have a value between 0 and 12......
!             0 and 12 mean this is a fresh start.  0 means that
!             dfault(iv, v) is to be called to provide all default
!             values to iv and v.  12 (the value that dfault assigns to
!             iv(1)) means the caller has already called dfault(iv, v)
!             and has possibly changed some iv and/or v entries to non-
!             default values.  default = 12.
! iv(covprt)... iv(14) = 1 means print a covariance matrix at the solu-
!             tion.  (this matrix is computed just before a return with
!             iv(1) = 3, 4, 5, 6.)
!             iv(covprt) = 0 means skip this printing.  default = 1.
! iv(covreq)... iv(15) = nonzero means compute a covariance matrix
!             just before a return with iv(1) = 3, 4, 5, 6.  in
!             this case, an approximate covariance matrix is obtained
!             in one of several ways.  let k = abs(iv(covreq)) and let
!             scale = 2*f(x)/max(1,n-p),  where 2*f(x) is the residual
!             sum of squares.  if k = 1 or 2, then a finite-difference
!             Hessian approximation h is obtained.  if h is positive
!             definite (or, for k = 3, if the jacobian matrix j at x
!             is nonsingular), then one of the following is computed...
!                  k = 1....  scale * h**-1 * (j**t * j) * h**-1.
!                  k = 2....  scale * h**-1.
!                  k = 3....  scale * (j**t * j)**-1.
!             (j**t is the transpose of j, while **-1 means inverse.)
!             if iv(covreq) is positive, then both function and grad-
!             ient values (calls on calcr and calcj) are used in com-
!             puting h (with step sizes determined using v(delta0) --
!             see below), while if iv(covreq) is negative, then only
!             function values (calls on calcr) are used (with step
!             sizes determined using v(dltfdc) -- see below).  if
!             iv(covreq) = 0, then no attempt is made to compute a co-
!             variance matrix (unless iv(covprt) = 1, in which case
!             iv(covreq) = 1 is assumed).  see iv(covmat) below.
!             default = 1.
! iv(dtype).... iv(16) tells how the scale vector d (see ref. 1) should
!             be chosen.  iv(dtype)>=1 means choose d as described
!             below with v(dfac).  iv(dtype)<=0 means the caller
!             has chosen d and has stored it in v starting at
!             v(94 + 2*n + p*(3*p + 31)/2).  default = 1.
! iv(inits).... iv(25) tells how the s matrix (see ref. 1) should be
!             initialized.  0 means initialize s to 0 (and start with
!             the gauss-newton model).  1 and 2 mean that the caller
!             has stored the lower triangle of the initial s rowwise in
!             v starting at v(87+2*p).  iv(inits) = 1 means start with
!             the gauss-newton model, while iv(inits) = 2 means start
!             with the augmented model (see ref. 1).  default = 0.
! iv(mxfcal)... iv(17) gives the maximum number of function evaluations
!             (calls on calcr, excluding those used to compute the co-
!             variance matrix) allowed.  if this number does not suf-
!             fice, then nl2sol returns with iv(1) = 9.  default = 200.
! iv(mxiter)... iv(18) gives the maximum number of iterations allowed.
!             it also indirectly limits the number of gradient evalua-
!             tions (calls on calcj, excluding those used to compute
!             the covariance matrix) to iv(mxiter) + 1.  if iv(mxiter)
!             iterations do not suffice, then nl2sol returns with
!             iv(1) = 10.  default = 150.
! iv(outlev)... iv(19) controls the number and length of iteration sum-
!             mary lines printed (by itsmry).  iv(outlev) = 0 means do
!             not print any summary lines.  otherwise, print a summary
!             line after each abs(iv(outlev)) iterations.  if iv(outlev)
!             is positive, then summary lines of length 117 (plus carri-
!             age control) are printed, including the following...  the
!             iteration and function evaluation counts, current func-
!             tion value (v(f) = half the sum of squares), relative
!             difference in function values achieved by the latest step
!             (i.e., reldf = (f0-v(f))/f0, where f0 is the function
!             value from the previous iteration), the relative function
!             reduction predicted for the step just taken (i.e.,
!             preldf = v(preduc) / f0, where v(preduc) is described
!             below), the scaled relative change in x (see v(reldx)
!             below), the models used in the current iteration (g =
!             gauss-newton, s=augmented), the marquardt parameter
!             stppar used in computing the last step, the sizing factor
!             used in updating s, the 2-norm of the scale vector d
!             times the step just taken (see ref. 1), and npreldf, i.e.,
!             v(nreduc)/f0, where v(nreduc) is described below -- if
!             npreldf is positive, then it is the relative function
!             reduction predicted for a newton step (one with
!             stppar = 0).  if npreldf is zero, either the gradient
!             vanishes (as does preldf) or else the augmented model
!             is being used and its Hessian is indefinite (with preldf
!             positive).  if npreldf is negative, then it is the nega-
!             of the relative function reduction predicted for a step
!             computed with step bound v(lmax0) for use in testing for
!             singular convergence.
!                  if iv(outlev) is negative, then lines of maximum
!             length 79 (or 55 is iv(covprt) = 0) are printed, includ-
!             ing only the first 6 items listed above (through reldx).
!             default = 1.
! iv(parprt)... iv(20) = 1 means print any nondefault v values on a
!             fresh start or any changed v values on a restart.
!             iv(parprt) = 0 means skip this printing.  default = 1.
! iv(prunit)... iv(21) is the output unit number on which all printing
!             is done.  iv(prunit) = -1 means suppress all printing.
!             (setting iv(prunit) to -1 is the only way to suppress the
!             one-line termination reason message printed by itsmry.)
!             default = standard output unit (unit 6 on most systems).
! iv(solprt)... iv(22) = 1 means print out the value of x returned (as
!             well as the corresponding gradient and scale vector d).
!             iv(solprt) = 0 means skip this printing.  default = 1.
! iv(statpr)... iv(23) = 1 means print summary statistics upon return-
!             ing.  these consist of the function value (half the sum
!             of squares) at x, v(reldx) (see below), the number of
!             function and gradient evaluations (calls on calcr and
!             calcj respectively, excluding any calls used to compute
!             the covariance), the relative function reductions predict-
!             ed for the last step taken and for a newton step (or per-
!             haps a step bounded by v(lmax0) -- see the descriptions
!             of preldf and npreldf under iv(outlev) above), and (if an
!             attempt was made to compute the covariance) the number of
!             calls on calcr and calcj used in trying to compute the
!             covariance.  iv(statpr) = 0 means skip this printing.
!             default = 1.
! iv(x0prt).... iv(24) = 1 means print the initial x and scale vector d
!             (on a fresh start only).  iv(x0prt) = 0 means skip this
!             printing.  default = 1.
!
!  ***  (selected) iv output values  ***
!
! iv(1)........ on output, iv(1) is a return code....
!             3 = x-convergence.  the scaled relative difference be-
!                  tween the current parameter vector x and a locally
!                  optimal parameter vector is very likely at most
!                  v(xctol).
!             4 = relative function convergence.  the relative differ-
!                  ence between the current function value and its lo-
!                  cally optimal value is very likely at most v(rfctol).
!             5 = both x- and relative function convergence (i.e., the
!                  conditions for iv(1) = 3 and iv(1) = 4 both hold).
!             6 = absolute function convergence.  the current function
!                  value is at most v(afctol) in absolute value.
!             7 = singular convergence.  the Hessian near the current
!                  iterate appears to be singular or nearly so, and a
!                  step of length at most v(lmax0) is unlikely to yield
!                  a relative function decrease of more than v(rfctol).
!             8 = false convergence.  the iterates appear to be converg-
!                  ing to a noncritical point.  this may mean that the
!                  convergence tolerances (v(afctol), v(rfctol),
!                  v(xctol)) are too small for the accuracy to which
!                  the function and gradient are being computed, that
!                  there is an error in computing the gradient, or that
!                  the function or gradient is discontinuous near x.
!             9 = function evaluation limit reached without other con-
!                  vergence (see iv(mxfcal)).
!            10 = iteration limit reached without other convergence
!                  (see iv(mxiter)).
!            11 = stopx returned .true. (external interrupt).  see the
!                  usage notes below.
!            13 = f(x) cannot be computed at the initial x.
!            14 = bad parameters passed to assess (which should not
!                  occur).
!            15 = the jacobian could not be computed at x (see calcj
!                  above).
!            16 = n or p (or parameter nn to nl2itr) out of range --
!                  p<=0 or n<p or nn<n.
!            17 = restart attempted with n or p (or par. nn to nl2itr)
!                  changed.
!            18 = iv(inits) is out of range.
!            19...45 = v(iv(1)) is out of range.
!            50 = iv(1) was out of range.
!            87...(86+p) = jtol(iv(1)-86) (i.e., v(iv(1)) is not
!                  positive (see v(dfac) below).
! iv(covmat)... iv(26) tells whether a covariance matrix was computed.
!             if (iv(covmat) is positive, then the lower triangle of
!             the covariance matrix is stored rowwise in v starting at
!             v(iv(covmat)).  if iv(covmat) = 0, then no attempt was
!             made to compute the covariance.  if iv(covmat) = -1,
!             then the finite-difference Hessian was indefinite.  and
!             and if iv(covmat) = -2, then a successful finite-differ-
!             encing step could not be found for some component of x
!             (i.e., calcr set nf to 0 for each of two trial steps).
!             note that iv(covmat) is reset to 0 after each successful
!             step, so if such a step is taken after a restart, then
!             the covariance matrix will be recomputed.
! iv(d)........ iv(27) is the starting subscript in v of the current
!             scale vector d.
! iv(g)........ iv(28) is the starting subscript in v of the current
!             least-squares gradient vector (j**t)*r.
! iv(nfcall)... iv(6) is the number of calls so far made on calcr (i.e.,
!             function evaluations, including those used in computing
!             the covariance).
! iv(nfcov).... iv(40) is the number of calls made on calcr when
!             trying to compute covariance matrices.
! iv(ngcall)... iv(30) is the number of gradient evaluations (calls on
!             calcj) so far done (including those used for computing
!             the covariance).
! iv(ngcov).... iv(41) is the number of calls made on calcj when
!             trying to compute covariance matrices.
! iv(niter).... iv(31) is the number of iterations performed.
! iv(r)........ iv(50) is the starting subscript in v of the residual
!             vector r corresponding to x.
!
!  ***  (selected) v input values (from subroutine dfault)  ***
!
! v(afctol)... v(31) is the absolute function convergence tolerance.
!             if nl2sol finds a point where the function value (half
!             the sum of squares) is less than v(afctol), and if nl2sol
!             does not return with iv(1) = 3, 4, or 5, then it returns
!             with iv(1) = 6.  default = max(10**-20, machep**2), where
!             machep is the unit roundoff.
! v(delta0)... v(44) is a factor used in choosing the finite-difference
!             step size used in computing the covariance matrix when
!             iv(covreq) = 1 or 2.  for component i, step size
!                  v(delta0) * max(abs(x(i)), 1/d(i)) * sign(x(i))
!             is used, where d is the current scale vector (see ref. 1).
!             (if this step results in calcr setting nf to 0, then -0.5
!             times this step is also tried.)  default = machep**0.5,
!             where machep is the unit roundoff.
! v(dfac)..... v(41) and the d0 and jtol arrays (see v(d0init) and
!             v(jtinit)) are used in updating the scale vector d when
!             iv(dtype)>0.  (d is initialized according to
!             v(dinit).)  let d1(i) =
!               max(sqrt(jcnorm(i)**2 + max(s(i,i),0)), v(dfac)*d(i)),
!             where jcnorm(i) is the 2-norm of the i-th column of the
!             current jacobian matrix and s is the s matrix of ref. 1.
!             if iv(dtype) = 1, then d(i) is set to d1(i) unless
!             d1(i)<jtol(i), in which case d(i) is set to
!                                max(d0(i), jtol(i)).
!             if iv(dtype)>=2, then d is updated during the first
!             iteration as for iv(dtype) = 1 (after any initialization
!             due to v(dinit)) and is left unchanged thereafter.
!             default = 0.6.
! v(dinit).... v(38), if nonnegative, is the value to which the scale
!             vector d is initialized.  default = 0.
! v(dltfdc)... v(40) helps choose the step size used when computing the
!             covariance matrix when iv(covreq) = -1 or -2.  for
!             differences involving x(i), the step size first tried is
!                       v(dltfdc) * max(abs(x(i)), 1/d(i)),
!             where d is the current scale vector (see ref. 1).  (if
!             this step is too big the first time it is tried, i.e., if
!             calcr sets nf to 0, then -0.5 times this step is also
!             tried.)  default = machep**(1/3), where machep is the
!             unit roundoff.
! v(d0init)... v(37), if positive, is the value to which all components
!             of the d0 vector (see v(dfac)) are initialized.  if
!             v(dfac) = 0, then it is assumed that the caller has
!             stored d0 in v starting at v(p+87).  default = 1.0.
! v(jtinit)... v(39), if positive, is the value to which all components
!             of the jtol array (see v(dfac)) are initialized.  if
!             v(jtinit) = 0, then it is assumed that the caller has
!             stored jtol in v starting at v(87).  default = 10**-6.
! v(lmax0).... v(35) gives the maximum 2-norm allowed for d times the
!             very first step that nl2sol attempts.  it is also used
!             in testing for singular convergence -- if the function
!             reduction predicted for a step of length bounded by
!             v(lmax0) is at most v(rfctol) * abs(f0), where  f0 is
!             the function value at the start of the current iteration,
!             and if nl2sol does not return with iv(1) = 3, 4, 5, or 6,
!             then it returns with iv(1) = 7.    default = 100.
! v(rfctol)... v(32) is the relative function convergence tolerance.
!             if the current model predicts a maximum possible function
!             reduction (see v(nreduc)) of at most v(rfctol)*abs(f0) at
!             the start of the current iteration, where  f0 is the
!             then current function value, and if the last step attempt-
!             ed achieved no more than twice the predicted function
!             decrease, then nl2sol returns with iv(1) = 4 (or 5).
!             default = max(10**-10, machep**(2/3)), where machep is
!             the unit roundoff.
! v(tuner1)... v(26) helps decide when to check for false convergence
!             and to consider switching models.  this is done if the
!             actual function decrease from the current step is no more
!             than v(tuner1) times its predicted value.  default = 0.1.
! v(xctol).... v(33) is the x-convergence tolerance.  if a newton step
!             (see v(nreduc)) is tried that has v(reldx)<=v(xctol)
!             and if this step yields at most twice the predicted func-
!             tion decrease, then nl2sol returns with iv(1) = 3 (or 5).
!             (see the description of v(reldx) below.)
!             default = machep**0.5, where machep is the unit roundoff.
! v(xftol).... v(34) is the false convergence tolerance.  if a step is
!             tried that gives no more than v(tuner1) times the predict-
!             ed function decrease and that has v(reldx)<=v(xftol),
!             and if nl2sol does not return with iv(1) = 3, 4, 5, 6, or
!             7, then it returns with iv(1) = 8.  (see the description
!             of v(reldx) below.)  default = 100*machep, where
!             machep is the unit roundoff.
! v(*)........ dfault supplies to v a number of tuning constants, with
!             which it should ordinarily be unnecessary to tinker.  see
!             version 2.2 of the nl2sol usage summary (which is an
!             appendix to ref. 1).
!
!  ***  (selected) v output values  ***
!
! v(dgnorm)... v(1) is the 2-norm of (d**-1)*g, where g is the most re-
!             cently computed gradient and d is the corresponding scale
!             vector.
! v(dstnrm)... v(2) is the 2-norm of d*step, where step is the most re-
!             cently computed step and d is the current scale vector.
! v(f)........ v(10) is the current function value (half the sum of
!             squares).
! v(f0)....... v(13) is the function value at the start of the current
!             iteration.
! v(nreduc)... v(6), if positive, is the maximum function reduction
!             possible according to the current model, i.e., the func-
!             tion reduction predicted for a newton step (i.e.,
!             step = -h**-1 * g,  where  g = (j**t) * r  is the current
!             gradient and h is the current Hessian approximation --
!             h = (j**t)*j  for the gauss-newton model and
!             h = (j**t)*j + s  for the augmented model).
!                  v(nreduc) = zero means h is not positive definite.
!                  if v(nreduc) is negative, then it is the negative of
!             the function reduction predicted for a step computed with
!             a step bound of v(lmax0) for use in testing for singular
!             convergence.
! v(preduc)... v(7) is the function reduction predicted (by the current
!             quadratic model) for the current step.  this (divided by
!             v(f0)) is used in testing for relative function
!             convergence.
! v(reldx).... v(17) is the scaled relative change in x caused by the
!             current step, computed as
!                  max(abs(d(i)*(x(i)-x0(i)), 1<=i<=p) /
!                     max(d(i)*(abs(x(i))+abs(x0(i))), 1<=i<=p),
!             where x = x0 + step.
!
!-------------------------------  notes  -------------------------------
!
!  ***  algorithm notes  ***
!
!        see ref. 1 for a description of the algorithm used.
!        on problems which are naturally well scaled, better perform-
!     ance may be obtained by setting v(d0init) = 1.0 and iv(dtype) = 0,
!     which will cause the scale vector d to be set to all ones.
!
!  ***  usage notes  ***
!
!        after a return with iv(1)<=11, it is possible to restart,
!     i.e., to change some of the iv and v input values described above
!     and continue the algorithm from the point where it was interrupt-
!     ed.  iv(1) should not be changed, nor should any entries of iv
!     and v other than the input values (those supplied by dfault).
!        those who do not wish to write a calcj which computes the ja-
!     cobian matrix analytically should call nl2sno rather than nl2sol.
!     nl2sno uses finite differences to compute an approximate jacobian.
!        those who would prefer to provide r and j (the residual and
!     jacobian) by reverse communication rather than by writing subrou-
!     tines calcr and calcj may call on nl2itr directly.  see the com-
!     ments at the beginning of nl2itr.
!        those who use nl2sol interactively may wish to supply their
!     own stopx function, which should return .true. if the break key
!     has been pressed since stopx was last invoked.  this makes it pos-
!     sible to externally interrupt nl2sol (which will return with
!     iv(1) = 11 if stopx returns .true.).
!        storage for j is allocated at the end of v.  thus the caller
!     may make v longer than specified above and may allow calcj to use
!     elements of j beyond the first n*p as scratch storage.
!
!  ***  portability notes  ***
!
!        the nl2sol distribution tape contains both single- and double-
!     precision versions of the nl2sol source code, so it should be un-
!     necessary to change precisions.
!        only the functions imdcon and rmdcon contain machine-dependent
!     constants.  to change from one machine to another, it should
!     suffice to change the (few) relevant lines in these functions.
!        intrinsic functions are explicitly declared.  on certain com-
!     puters (e.g. univac), it may be necessary to comment out these
!     declarations.  so that this may be done automatically by a simple
!     program, such declarations are preceded by a comment having c/+
!     in columns 1-3 and blanks in columns 4-72 and are followed by
!     a comment having c/ in columns 1 and 2 and blanks in columns 3-72.
!        the nl2sol source code is expressed in 1966 ansi standard
!     fortran.  it may be converted to fortran 77 by
!     commenting out all lines that fall between a line having c/6 in
!     columns 1-3 and a line having c/7 in columns 1-3 and by removing
!     (i.e., replacing by a blank) the c in column 1 of the lines that
!     follow the c/7 line and preceed a line having c/ in columns 1-2
!     and blanks in columns 3-72.  these changes convert some data
!     statements into parameter statements, convert some variables from
!     real to character*4, and make the data statements that initialize
!     these variables use character strings delimited by primes instead
!     of hollerith constants.  (such variables and data statements
!     appear only in modules itsmry and parchk.  parameter statements
!     appear nearly everywhere.)
!
!  ***  references  ***
!
! 1.  dennis, j.e., gay, d.m., and welsch, r.e. (1981), an adaptive
!             nonlinear least-squares algorithm, acm trans. math.
!             software, vol. 7, no. 3.
!
!  ***  general  ***
!
!     coded by david m. gay (winter 1979 - winter 1980).
!     this subroutine was written in connection with research
!     supported by the national science foundation under grants
!     mcs-7600324, dcr75-10143, 76-14311dss, mcs76-11989, and
!     mcs-7906671.
!
!----------------------------  declarations  ---------------------------

! EXTERNAL itsmry,nl2itr
! itsmry... prints iteration summary and info about initial and final x.
! nl2itr... reverse-communication routine that carries out nl2sol algorithm.

 LOGICAL strted
 INTEGER d1,j1,nf,r1
! iv subscript values
 INTEGER,PARAMETER :: nfcall=6,nfgcal=7,toobig=2
! v subscript values
 INTEGER,PARAMETER :: d=27,j=33,r=50

 d1=94+2*n+p*(3*p+31)/2
 iv(d)=d1
 r1=d1+p
 iv(r)=r1
 j1=r1+n
 iv(j)=j1
 strted=.true.
 if(iv(1)/=0.and.iv(1)/=12)goto 40
 strted=.false.
 iv(nfcall)=1
 iv(nfgcal)=1

10 nf=iv(nfcall)
 call calcr(n,p,x,nf,v(r1:r1+n-1),uiparm,urparm,ufparm)
 if(strted)goto 20
 if(nf>0)goto 30
 iv(1)=13
 goto 60


20 if(nf<=0)iv(toobig)=1
 goto 40

30 call calcj(n,p,x,iv(nfgcal),v(j1:j1+n*p-1),uiparm,urparm,ufparm)
 if(iv(nfgcal)==0)goto 50
 strted=.true.

40 call nl2itr(v(d1:d1+p-1),iv,v(j1:j1+n*p-1),n,n,p,v(r1:r1+n-1),v(1:d1-1),x)
 if(iv(1)-2<0)then
  goto 10
 elseif(iv(1)-2==0)then
  goto 30
 else
  goto 70
 endif

50 iv(1)=15
60 call itsmry(v(d1:),iv,p,v,x)

70 return

END SUBROUTINE nl2sol


SUBROUTINE nl2sno(n,p,x,calcr,iv,v,uiparm,urparm,ufparm)
!-----------------------------------------------------------------------!
! Like NL2SOL, but without calcj -- minimize nonlinear sum of           !
! squares using finite-difference jacobian approximations               !
! (NL2SOL version 2.2)                                                  !
!-----------------------------------------------------------------------!

 INTEGER,INTENT(in) :: n,p
 INTEGER,INTENT(inout) :: iv(:)
 INTEGER,INTENT(inout),OPTIONAL :: uiparm(:)
 REAL(dp),INTENT(inout) :: x(:),v(:)
 REAL(dp),INTENT(inout),OPTIONAL :: urparm(:),ufparm
 REAL(dp) x_h(p),h_x(p)

 INTERFACE
  SUBROUTINE calcr(n,p,x,nf,r,uiparm,urparm,ufparm)
   USE machine_constants
   IMPLICIT NONE
   INTEGER,INTENT(in) :: n,p
   INTEGER,INTENT(inout) :: nf
   INTEGER,INTENT(inout),OPTIONAL :: uiparm(:)
   REAL(dp),INTENT(in) :: x(:)
   REAL(dp),INTENT(inout),OPTIONAL :: urparm(:),ufparm
   REAL(dp),INTENT(out) :: r(:)
  END SUBROUTINE calcr
 END INTERFACE

!-----------------------------  discussion  ----------------------------
!
!        The parameters for nl2sno are the same as those for nl2sol
!     (which see), except that calcj is omitted.  instead of calling
!     calcj to obtain the jacobian matrix of r at x, nl2sno computes
!     an approximation to it by finite (forward) differences -- see
!     v(dltfdj) below.  nl2sno uses function values only when comput-
!     the covariance matrix (rather than the functions and gradients
!     that nl2sol may use).  to do so, nl2sno sets iv(covreq) to -1 if
!     iv(covprt) = 1 with iv(covreq) = 0 and to minus its absolute
!     value otherwise.  thus v(delta0) is never referenced and only
!     v(dltfdc) matters -- see nl2sol for a description of v(dltfdc).
!        The number of extra calls on calcr used in computing the jaco-
!     bian approximation are not included in the function evaluation
!     count iv(nfcall) and are not otherwise reported.
!
! v(dltfdj)... v(36) helps choose the step size used when computing the
!             finite-difference jacobian matrix.  for differences in-
!             volving x(i), the step size first tried is
!                       v(dltfdj) * max(abs(x(i)), 1/d(i)),
!             where d is the current scale vector (see ref. 1).  (if
!             this step is too big, i.e., if calcr sets nf to 0, then
!             smaller steps are tried until the step size is shrunk be-
!             low 1000 * machep, where machep is the unit roundoff.
!             default = machep**0.5.
!
!  ***  references  ***
!
! 1.  dennis, j.e., gay, d.m., and welsch, r.e. (1981), an adaptive
!             nonlinear least-squares algorithm, acm trans. math.
!             software, vol. 7, no. 3.
!
!  ***  general  ***
!
!     coded by David M. Gay.
!     This subroutine was written in connection with research
!     supported by the national science foundation under grants
!     mcs-7600324, dcr75-10143, 76-14311dss, mcs76-11989, and
!     mcs-7906671.
!
!+++++++++++++++++++++++++++  declarations  ++++++++++++++++++++++++++++
!
! dfault... supplies default parameter values.
! itsmry... prints iteration summary and info about initial and final x.
! nl2itr... reverse-communication routine that carries out nl2sol algo-
!             rithm.
! rmdcon... returns machine-dependent constants.
! vscopy... sets all elements of a vector to a scalar.

 LOGICAL strted,bad_point
 INTEGER dk,d1,i,j1,j1k,k,nf,rn,r1
 REAL(dp) h,xk
 REAL(dp),PARAMETER :: hfac=1.d3,negpt5=-0.5d0
 REAL(dp),PARAMETER :: one=1.d0,zero=0.d0
 INTEGER,PARAMETER :: covprt=14,covreq=15,d=27,dtype=16,j=33,nfcall=6,nfgcal=7,&
  &r=50,toobig=2
 INTEGER,PARAMETER :: dltfdj=36,dinit=38
 REAL(dp),SAVE :: hlim=0.d0

! Initialize the iteration counter to zero here. (NDD).
 iv(31)=0

 d1=94+2*n+p*(3*p+31)/2
 iv(d)=d1
 r1=d1+p
 iv(r)=r1
 j1=r1+n
 iv(j)=j1
 rn=j1-1
 if(iv(1)==0)call dfault(iv,v)
 iv(covreq)=-abs(iv(covreq))
 if(iv(covprt)/=0.and.iv(covreq)==0)iv(covreq)=-1
 strted=.true.
 if(iv(1)/=12)goto 80
 strted=.false.
 iv(nfcall)=1
 iv(nfgcal)=1
! Initialize scale vector d to ones for computing initial jacobian.
 if(iv(dtype)>0)call vscopy(p,v(d1:d1+p-1),one)
 if(v(dinit)>zero)call vscopy(p,v(d1:d1+p-1),v(dinit))

10 nf=iv(nfcall)
 bad_point=.false.
 call calcr(n,p,x,nf,v(r1:r1+n-1),uiparm,urparm,ufparm)
 if(bad_point)iv(2)=1
 if(strted)goto 20
 if(nf>0)goto 30
 iv(1)=13
 goto 90

20 if(nf<=0)iv(toobig)=1
 goto 80

! Compute finite-difference Jacobian.

30 j1k=j1
 dk=d1
 do k=1,p
  xk=x(k)
  h=v(dltfdj)*max(abs(xk),one/v(dk))
  dk=dk+1
40 x(k)=xk+h
  nf=iv(nfgcal)
  call calcr(n,p,x,nf,v(j1k:j1k+n-1),uiparm,urparm,ufparm)
  if(nf>0)goto 50
  if(hlim==zero)hlim=hfac*rmdcon(3)
! hlim = hfac times the unit roundoff
  h=negpt5*h
  if(abs(h)>=hlim)goto 40
  iv(1)=15
  goto 90
50 x(k)=xk
  do i=r1,rn
   v(j1k)=(v(j1k)-v(i))/h
   j1k=j1k+1
  enddo
 enddo
!!$30 j1k=j1
!!$ dk=d1
!!$ do k=1,p
!!$  xk=x(k)
!!$  h=v(dltfdj)*max(abs(xk),one/v(dk))
!!$  dk=dk+1
!!$  h_x(k)=h
!!$  x_h(k)=xk+h
!!$ enddo
!!$ nf=iv(nfgcal)
!!$ call calcr(n,p,x,nf,v(j1:j1+n*p-1),uiparm,urparm,ufparm)
!!$ j1k=j1
!!$ do k=1,p
!!$  h=h_x(k)
!!$  do i=r1,rn
!!$   v(j1k)=(v(j1k)-v(i))/h
!!$   j1k=j1k+1
!!$  enddo
!!$ enddo

 strted=.true.

80 call nl2itr(v(d1:d1+p-1),iv,v(j1:j1+n*p-1),n,n,p,v(r1:r1+n-1),v(1:d1-1),x)
 if(iv(1)-2<0)then
  goto 10
 elseif(iv(1)-2==0)then
  goto 30
 else
  goto 100
 endif

90 call itsmry(v(d1:),iv,p,v,x)

100 return

END SUBROUTINE nl2sno


SUBROUTINE nl2itr(d,iv,jac,n,nn,p,r,v,x)
!--------------------------------------------------------------!
! Carry out NL2SOL (nonlinear least-squares) iterations.       !
! (NL2SOL version 2.2)                                         !
!--------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: n, nn, p
 INTEGER,INTENT(inout) :: iv(:)
 REAL(dp),INTENT(inout) :: d(:),jac(:),r(:),v(:),x(:)
!--------------------------  parameter usage  --------------------------
!
! d.... scale vector.
! iv... integer value array.
! j.... n by p jacobian matrix (lead dimension nn).
! n.... number of observations (components in r).
! nn... lead dimension of j.
! p.... number of parameters (components in x).
! r.... residual vector.
! v.... floating-point value array.
! x.... parameter vector.
!
!  ***  discussion  ***
!
!        parameters iv, n, p, v, and x are the same as the correspond-
!     ing ones to nl2sol (which see), except that v can be shorter
!     (since the part of v that nl2sol uses for storing d, j, and r is
!     not needed).  moreover, compared with nl2sol, iv(1) may have the
!     two additional output values 1 and 2, which are explained below,
!     as is the use of iv(toobig) and iv(nfgcal).  the values iv(d),
!     iv(j), and iv(r), which are output values from nl2sol (and
!     nl2sno), are not referenced by nl2itr or the subroutines it calls.
!        on a fresh start, i.e., a call on nl2itr with iv(1) = 0 or 12,
!     nl2itr assumes that r = r(x), the residual at x, and j = j(x),
!     the corresponding jacobian matrix of r at x.
!
! iv(1) = 1 means the caller should set r to r(x), the residual at x,
!             and call nl2itr again, having changed none of the other
!             parameters.  an exception occurs if r cannot be evaluated
!             at x (e.g. if r would overflow), which may happen because
!             of an oversized step.  in this case the caller should set
!             iv(toobig) = iv(2) to 1, which will cause nl2itr to ig-
!             nore r and try a smaller step.  the parameter nf that
!             nl2sol passes to calcr (for possible use by calcj) is a
!             copy of iv(nfcall) = iv(6).
! iv(1) = 2 means the caller should set j to j(x), the jacobian matrix
!             of r at x, and call nl2itr again.  the caller may change
!             d at this time, but should not change any of the other
!             parameters.  the parameter nf that nl2sol passes to
!             calcj is iv(nfgcal) = iv(7).  if j cannot be evaluated
!             at x, then the caller may set iv(nfgcal) to 0, in which
!             case nl2itr will return with iv(1) = 15.
!
!  ***  general  ***
!
!     coded by david m. gay.
!     this subroutine was written in connection with research
!     supported by the national science foundation under grants
!
!     mcs-7600324, dcr75-10143, 76-14311dss, mcs76-11989, and
!     mcs-7906671.
!        (see nl2sol for references.)
!
 INTEGER dig1,g1,g01,h0,h1,i,im1,ipivi,ipivk,ipiv1,ipk,k,km1,&
  &l,lky1,lmat1,lstgst,m,pp1o2,qtr1,rdk,rd0,rd1,rsave1,smh,sstep,&
  &step1,stpmod,s1,temp1,temp2,w1,x01
 REAL(dp) e,rdof1,sttsst,t,t1
 REAL(dp),ALLOCATABLE :: j(:,:)

! assess... assesses candidate step.
! covclc... computes covariance matrix.
! dotprd... returns inner product of two vectors.
! dupdat... updates scale vector d.
! gqtstp... computes goldfeld-quandt-trotter step (augmented model).
! itsmry... prints iteration summary and info about initial and final x.
! lmstep... computes levenberg-marquardt step (gauss-newton model).
! parchk... checks validity of input iv and v values.
! qapply... applies orthogonal matrix q from qrfact to a vector.
! qrfact... computes qr decomposition of a matrix via householder trans.
! rptmul... multiplies vector by the r matrix (and/or its transpose)
!             stored by qrfact.
! slupdt... performs quasi-newton update on compactly stored lower tri-
!             angle of a symmetric matrix.
! stopx.... returns .true. if the break key has been pressed.
! vaxpy.... computes scalar times one vector plus another.
! vcopy.... copies one vector to another.
! vscopy... sets all elements of a vector to a scalar.
! v2norm... returns the 2-norm of a vector.

!  ***  iv subscript values  ***
 INTEGER,PARAMETER :: cnvcod=34,covmat=26,covprt=14,covreq=15,dig=43,&
  &dtype=16,g=28,h=44,ierr=32,inits=25,ipivot=61,ipiv0=60,irc=3,kagqt=35,&
  &kalm=36,lky=37,lmat=58,mode=38,model=5,mxfcal=17,mxiter=18,nfcall=6,&
  &nfgcal=7,nfcov=40,ngcov=41,ngcall=30,niter=31,qtr=49,radinc=8,rd=51,&
  &restor=9,rsave=52,s=53,step=55,stglim=11,stlstg=56,sused=57,switch=12,&
  &toobig=2,w=59,xirc=13,x0=60

 INTEGER,PARAMETER :: cosmin=43,dgnorm=1,dinit=38,dstnrm=2,d0init=37,&
  &f=10,fdif=11,fuzz=45,f0=13,gtstep=4,incfac=23,&
  &jtinit=39,jtol1=87,lmax0=35,nvsave=9,phmxfc=21,&
  &preduc=7,radfac=16,radius=8,rad0=9,rlimit=42,&
  &size=47,stppar=5,tuner4=29,tuner5=30,vsave1=78,wscale=48

 REAL(dp),PARAMETER :: half=0.5d0,negone=-1.d0,one=1.d0,zero=0.d0

 allocate(j(nn,p))
! Trouble with the RESHAPE intrinsic in ifort 8.1/9.0 (PLR 10.2005)
! j=reshape(jac,(/nn,p/))
 l=0
 do k=1,p
  do i=1,nn
   l=l+1 ; j(i,k)=jac(l)
  enddo
 enddo
 i=iv(1)
 if(i==1)goto 20
 if(i==2)goto 50

! Check validity of iv and v input values.

! Note -- if iv(1) = 0, then parchk calls dfault(iv, v)
 call parchk(iv,n,nn,p,v)
 i=iv(1)-2

 select case(i)
 case (1:6)   ; goto 300
 case (7,9)   ; goto 150
 case (8)     ; goto 100
 case (10)    ; goto 10
 case default ; goto 590
 endselect

! Initialization and storage allocation
 10 iv(niter)=0
 iv(nfcall)=1
 iv(ngcall)=1
 iv(nfgcal)=1
 iv(mode)=-1
 iv(stglim)=2
 iv(toobig)=0
 iv(cnvcod)=0
 iv(covmat)=0
 iv(nfcov)=0
 iv(ngcov)=0
 iv(kalm)=-1
 iv(radinc)=0
 iv(s)=jtol1+2*p
 pp1o2=p*(p+1)/2
 iv(x0)=iv(s)+pp1o2
 iv(step)=iv(x0)+p
 iv(stlstg)=iv(step)+p
 iv(dig)=iv(stlstg)+p
 iv(g)=iv(dig)+p
 iv(lky)=iv(g)+p
 iv(rd)=iv(lky)+p
 iv(rsave)=iv(rd)+p
 iv(qtr)=iv(rsave)+n
 iv(h)=iv(qtr)+n
 iv(w)=iv(h)+pp1o2
 iv(lmat)=iv(w)+4*p+7
! Length of w = p*(p+9)/2 + 7.  lmat is contained in w.
 if(v(dinit)>=zero)call vscopy(p,d(1:p),v(dinit))
 if(v(jtinit)>zero)call vscopy(p,v(jtol1:jtol1+p-1),v(jtinit))
 i=jtol1+p
 if(v(d0init)>zero)call vscopy(p,v(i:i+p-1),v(d0init))
 v(rad0)=zero
 v(stppar)=zero
 v(radius)=v(lmax0)/(one+v(phmxfc))

! Set initial model and s matrix.

 iv(model)=1
 if(iv(inits)==2)iv(model)=2
 s1=iv(s)
 if(iv(inits)==0)call vscopy(pp1o2,v(s1:s1+pp1o2-1),zero)

! Compute function value (half the sum of squares)

20 t=v2norm(n,r)
 if(t>v(rlimit))iv(toobig)=1
 if(iv(toobig)/=0)goto 30
 v(f)=half*t**2
30 if(iv(mode)<0)then
  goto 40
 elseif(iv(mode)==0)then
  goto 300
 else
  goto 540
 endif

40 if(iv(toobig)==0)goto 60
 iv(1)=13
 goto 580

! Make sure Jacobian could be computed.

50 if(iv(nfgcal)/=0)goto 60
 iv(1)=15
 goto 580

! Compute gradient.

60 iv(kalm)=-1
 g1=iv(g)
 do i=1,p
  v(g1)=dotprd(n,r,j(:,i))
  g1=g1+1
 enddo
 if(iv(mode)>0) goto 520

! Update d and make copies of r for possible use later.

 if(iv(dtype)>0)call dupdat(d(1:p),iv,j,n,p,v)
 rsave1=iv(rsave)
 call vcopy (n,v(rsave1:rsave1+n-1),r)
 qtr1=iv(qtr)
 call vcopy(n,v(qtr1:qtr1+n-1),r)

! Compute d**-1 * gradient.

 g1=iv(g)
 dig1=iv(dig)
 k=dig1
 do i=1,p
  v(k)=v(g1)/d(i)
  k=k+1
  g1=g1+1
 enddo
 v(dgnorm)=v2norm(p,v(dig1:))

 if(iv(cnvcod)/=0)goto 510
 if(iv(mode)==0)goto 460
 iv(mode)=0

!-----------------------------  main loop  -----------------------------

! Print iteration summary, check iteration limit.

 90 call itsmry(d,iv,p,v,x)
 100 k=iv(niter)
 if(k<iv(mxiter))goto 110
 iv(1)=10
 goto 580
 110 iv(niter)=k+1

! Update radius.

 if(k==0) goto 130
 step1=iv(step)
 do i=1,p
  v(step1)=d(i)*v(step1)
  step1=step1+1
 enddo
 step1=iv(step)
 v(radius)=v(radfac)*v2norm(p,v(step1:))

! Initialize for start of next iteration.

130 x01=iv(x0)
 v(f0)=v(f)
 iv(kagqt)=-1
 iv(irc)=4
 iv(h)=-abs(iv(h))
 iv(sused)=iv(model)

! Copy x to x0.

 call vcopy(p,v(x01:x01+p-1),x)

! Check stopx and function evaluation limit.

 140 if(.not.stopx())goto 160
 iv(1)=11
 goto 170

! Come here when restarting after func. eval. limit or stopx.

150 if(v(f)>=v(f0))goto 160
 v(radfac)=one
 k=iv(niter)
 goto 110

 160 if(iv(nfcall)<iv(mxfcal)+iv(nfcov))goto 180
 iv(1)=9
 170 if(v(f)>=v(f0))goto 580

! In case of stopx or function evaluation limit with
! improved v(f), evaluate the gradient at x.

 iv(cnvcod)=iv(1)
 goto 450

!. . . . . . . . . . . . .  compute candidate step  . . . . . . . . . .

 180 step1=iv(step)
 w1=iv(w)
 if(iv(model)==2)goto 220

! Compute levenberg-marquardt step.

 qtr1=iv(qtr)
 if(iv(kalm)>=0)goto 190
 rd1=iv(rd)
 if(-1==iv(kalm))call qrfact(n,p,j,v(rd1:rd1+p-1),iv(ipivot:ipivot+p-1),&
  &iv(ierr),0,v(w1:w1+p-1))
 call qapply(n,p,j,v(qtr1:),iv(ierr))
 190 h1=iv(h)
 if(h1>0)goto 210

! Copy r matrix to h.

 h1=-h1
 iv(h)=h1
 k=h1
 rd1=iv(rd)
 v(k)=v(rd1)
 if(p==1)goto 210
 do i=2,p
  call vcopy(i-1,v(k+1:k+1+i-1-1),j(:,i))
  k=k+i
  rd1=rd1+1
  v(k)=v(rd1)
 enddo

210 g1=iv(g)
 call lmstep(d,v(g1:g1+p-1),iv(ierr),iv(ipivot:ipivot+p-1),iv(kalm),p,&
  &v(qtr1:qtr1+n-1),v(h1:h1+(p*(p+1))/2-1),v(step1:step1+p-1),v(1:86),v(w1:))
 goto 290

! Compute goldfeld-quandt-trotter step (augmented model).

220 if(iv(h)>0)goto 280

! Set h to  d**-1 * ( (j**t)*j + s) ) * d**-1.

 h1=-iv(h)
 iv(h)=h1
 s1=iv(s)
 if(-1/=iv(kalm))goto 250

! j is in its original form.

 do i=1,p
  t=one/d(i)
  do k=1,i
   v(h1)=t*(dotprd(n,j(:,i),j(:,k))+v(s1))/d(k)
   h1=h1+1
   s1=s1+1
  enddo
 enddo
 goto 280

! lmstep has applied qrfact to j

250 smh=s1-h1
 h0=h1-1
 ipiv1=iv(ipivot)
 t1=one/d(ipiv1)
 rd0=iv(rd)-1
 rdof1=v(rd0+1)
 do i=1,p
  l=ipiv0+i
  ipivi=iv(l)
  h1=h0+ipivi*(ipivi-1)/2
  l=h1+ipivi
  m=l+smh
! ***  v(l) = h(ipivot(i), ipivot(i))  ***
! ***  v(m) = s(ipivot(i), ipivot(i))  ***
  t=one/d(ipivi)
  rdk=rd0+i
  e=v(rdk)**2
  if(i>1)e=e+dotprd(i-1,j(:,i),j(:,i))
  v(l)=(e+v(m))*t**2
  if(i==1)cycle
  l=h1+ipiv1
  if(ipivi<ipiv1)l=l+((ipiv1-ipivi)*(ipiv1+ipivi-3))/2
  m=l+smh
! ***  v(l) = h(ipivot(i), ipivot(1))  ***
! ***  v(m) = s(ipivot(i), ipivot(1))  ***
  v(l)=t*(rdof1*j(1,i)+v(m))*t1
  if(i==2)cycle
  im1=i-1
  do k=2,im1
   ipk=ipiv0+k
   ipivk=iv(ipk)
   l=h1+ipivk
   if(ipivi<ipivk)l=l+((ipivk-ipivi)*(ipivk+ipivi-3))/2
   m=l+smh
! ***  v(l) = h(ipivot(i), ipivot(k))  ***
! ***  v(m) = s(ipivot(i), ipivot(k))  ***
   km1=k-1
   rdk=rd0+k
   v(l)=t*(dotprd(km1,j(:,i),j(:,k))+v(rdk)*j(k,i)+v(m))/d(ipivk)
  enddo
 enddo

! Compute actual Goldfeld-Quandt-Trotter step.

280 h1=iv(h)
 dig1=iv(dig)
 lmat1=iv(lmat)
 call gqtstp(d,v(dig1:dig1+p-1),v(h1:h1+(p*(p+1))/2-1),iv(kagqt),&
  &v(lmat1:lmat1+(p*(p+1))/2-1),p,v(step1:step1+p-1),v(1:86),v(w1:w1+4*p+7-1))

! Compute r(x0 + step).

290 if(iv(irc)==6)goto 300
 x01=iv(x0)
 step1=iv(step)
 call vaxpy(p,x(1:p),one,v(step1:step1+p-1),v(x01:x01+p-1))
 iv(nfcall)=iv(nfcall)+1
 iv(1)=1
 iv(toobig)=0
 goto 590

! Assess candidate step.

300 step1=iv(step)
 lstgst=iv(stlstg)
 x01=iv(x0)
 call assess(d,iv,p,v(step1:step1+p-1),v(lstgst:lstgst+p-1),v(1:86),x,&
  &v(x01:x01+p-1))

! If necessary, switch models and/or restore r.

 if(iv(switch)==0)goto 310
 iv(h)=-abs(iv(h))
 iv(sused)=iv(sused)+2
 call vcopy(nvsave,v(1:nvsave),v(vsave1:))
310 if(iv(restor)==0)goto 320
 rsave1=iv(rsave)
 call vcopy(n,r(1:n),v(rsave1:))
320 l=iv(irc)-4
 stpmod=iv(model)
 if(l>0)then
  select case(l)
  case (1)   ; goto 340
  case (2)   ; goto 360
  case (3:8) ; goto 370
  case (9)   ; goto 500
  case (10)  ; goto 460
  endselect
 endif

! Decide whether to change models.

 e=v(preduc)-v(fdif)
 sstep=iv(lky)
 s1=iv(s)
 call slvmul(p,v(sstep:sstep+p*(p+1)/2-1),v(s1:s1+(p*(p+1))/2-1),&
  &v(step1:step1+p-1))
 sttsst=half*dotprd(p,v(step1:),v(sstep:))
 if(iv(model)==1)sttsst=-sttsst
 if(abs(e+sttsst)*v(fuzz)>=abs(e))goto 330

! Switch models.

 iv(model)=3-iv(model)
 if(iv(model)==1)iv(kagqt)=-1
 if(iv(model)==2.and.iv(kalm)>0)iv(kalm)=0
 if(-2<l) goto 380
 iv(h)=-abs(iv(h))
 iv(sused)=iv(sused)+2
 call vcopy(nvsave,v(vsave1:vsave1+nvsave-1),v)
 goto 350

330 if(-3<l)goto 380

! Recompute step with decreased radius.

 v(radius)=v(radfac)*v(dstnrm)
 goto 140

! Recompute step, saving v values and r if necessary.

340 v(radius)=v(radfac)*v(dstnrm)
350 if(v(f)>=v(f0))goto 140
 rsave1=iv(rsave)
 call vcopy(n,v(rsave1:rsave1+n-1),r)
 goto 140

! Compute step of length v(lmax0) for singular convergence test.

360 v(radius)=v(lmax0)
 goto 180

! Convergence or false convergence.

370 iv(cnvcod)=l
 if(v(f)>=v(f0))goto 510
 if(iv(xirc)==14)goto 510
 iv(xirc)=14

! Process acceptable step.

380 iv(covmat)=0

! Set lky = (j(x0)**t) * r(x).

 lky1=iv(lky)
 if(iv(kalm)>=0)goto 400

! Jacobian has not been modified.

 do i=1,p
  v(lky1)=dotprd(n,j(:,i),r)
  lky1=lky1+1
 enddo
 goto 410

! qrfact has been applied to j.  store copy of r in qtr and apply q to it.

400 qtr1=iv(qtr)
 call vcopy(n,v(qtr1:qtr1+n-1),r)
 call qapply(n,p,j,v(qtr1:qtr1+n-1),iv(ierr))

! Multiply top p-vector in qtr by permuted upper triangle stored by qrfact
! in j and rd.

 rd1=iv(rd)
 temp1=iv(stlstg)
 call rptmul(3,iv(ipivot:ipivot+p-1),j,p,v(rd1:rd1+p-1),v(qtr1:qtr1+p-1),&
  &v(lky1:lky1+p-1),v(temp1:temp1+p-1))

! See whether to set v(radfac) by gradient tests.

410 if(iv(irc)/=3)goto 450
 step1=iv(step)
 temp1=iv(stlstg)
 temp2=iv(x0)

! Set temp1 = hessian * step  for use in gradient tests.

 if(stpmod==2)goto 420

! Step computed using gauss-newton model -- qrfact has been applied to j

 rd1=iv(rd)
 call rptmul(2,iv(ipivot:ipivot+p-1),j,p,v(rd1:rd1+p-1),v(step1:step1+p-1),&
  &v(temp1:temp1+p-1),v(temp2:temp2+p-1))
 goto 450

! Step computed using augmented model.

420 h1=iv(h)
 k=temp2
 do i=1,p
  v(k)=d(i)*v(step1)
  k=k+1
  step1=step1+1
 enddo
 call slvmul(p,v(temp1:temp1+p*(p+1)/2-1),v(h1:),v(temp2:))
 do i=1,p
  v(temp1)=d(i)*v(temp1)
  temp1=temp1+1
 enddo

! Save old gradient and compute new one.

450 iv(ngcall)=iv(ngcall)+1
 g1=iv(g)
 g01=iv(w)
 call vcopy(p,v(g01:g01+p-1),v(g1:))
 iv(1)=2
 goto 590

! Initializations -- g0 = g - g0, etc.

460 g01=iv(w)
 g1=iv(g)
 call vaxpy(p,v(g01:g01+p-1),negone,v(g01:),v(g1:))
 step1=iv(step)
 temp1=iv(stlstg)
 temp2=iv(x0)
 if(iv(irc)/=3)goto 490

! Set v(radfac) by gradient tests.

! Set  temp1 = d**-1 * (hessian * step  +  (g(x0) - g(x))).

 k=temp1
 l=g01
 do i=1,p
  v(k)=(v(k)-v(l))/d(i)
  k=k+1
  l=l+1
 enddo

! Do gradient tests.

 if(v2norm(p,v(temp1:))<=v(dgnorm)*v(tuner4))goto 480
 if(dotprd(p,v(g1:),v(step1:))>=v(gtstep)*v(tuner5))goto 490
480 v(radfac)=v(incfac)

! Finish computing lky = ((j(x) - j(x0))**t) * r.

! Currently lky = (j(x0)**t) * r.

490 lky1=iv(lky)
 call vaxpy (p,v(lky1:lky1+p-1),negone,v(lky1:),v(g1:))

! Determine sizing factor v(size).

! Set temp1 = s * step.
 s1=iv(s)
 call slvmul(p,v(temp1:temp1+(p*(p+1)/2)-1),v(s1:),v(step1:))

 t1=abs(dotprd(p,v(step1:),v(temp1:)))
 t=abs(dotprd(p,v(step1:),v(lky1:)))
 v(size)=one
 if(t<t1)v(size)=t/t1

! Update s.

 call slupdt(v(s1:s1+(p*(p+1))/2-1),v(cosmin),p,v(size),v(step1:step1+p-1),&
  &v(temp1:temp1+p-1),v(temp2:temp2+p-1),v(g01:g01+p-1),v(wscale),&
  &v(lky1:lky1+p-1))
 iv(1)=2
 goto 90

! Misc. details.

! Bad parameters to assess.

500 iv(1)=14
 goto 580

! Convergence obtained -- compute covariance matrix if desired.

510 if(iv(covreq)==0.and.iv(covprt)==0)goto 570
 if(iv(covmat)/=0)goto 570
 if(iv(cnvcod)>=7)goto 570
 iv(mode)=0
520 call covclc(i,d,iv,j,n,p,r,v,x)

 select case(i)
 case(1:2)
  iv(nfcov)=iv(nfcov)+1
  iv(nfcall)=iv(nfcall)+1
  iv(restor)=i
  iv(1)=1
  goto 590
 case(3) ; goto 550
 case(4)
  iv(mode)=0
  if(iv(niter)==0)iv(mode)=-1
  goto 570
 endselect

540 if(iv(restor)==1.or.iv(toobig)/=0)goto 520
 iv(nfgcal)=iv(nfcall)
550 iv(ngcov)=iv(ngcov)+1
 iv(ngcall)=iv(ngcall)+1
 iv(1)=2
 goto 590

570 iv(1)=iv(cnvcod)
 iv(cnvcod)=0

! Print summary of final iteration and other requested items.

580 call itsmry(d,iv,p,v,x)

590 l=0
! Trouble with the RESHAPE intrinsic in ifort 8.1/9.0 (PLR 10.2005)
! jac(:nn*p)=reshape(j,(/nn*p/))
 do k=1,p
  do i=1,nn
   l=l+1 ; jac(l)=j(i,k)
  enddo
 enddo
 deallocate(j)

END SUBROUTINE nl2itr


SUBROUTINE assess(d,iv,p,step,stlstg,v,x,x0)
!-----------------------------------------------!
! Assess candidate step (nl2sol version 2.2).   !
!-----------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: p
 INTEGER,INTENT(inout) :: iv(:)
 REAL(dp),INTENT(in) :: d(:),x0(:)
 REAL(dp),INTENT(inout) :: step(:),stlstg(:),v(:),x(:)

! Purpose
!        this subroutine is called by an unconstrained minimization
!     routine to assess the next candidate step.  it may recommend one
!     of several courses of action, such as accepting the step, recom-
!     puting it using the same or a new quadratic model, or halting due
!     to convergence or false convergence.  see the return code listing
!     below.

!--------------------------  parameter usage  --------------------------

!     iv (i/o) integer parameter and scratch vector -- see description
!             below of iv values referenced.
!      d (in)  scale vector used in computing v(reldx) -- see below.
!      p (in)  number of parameters being optimized.
!   step (i/o) on input, step is the step to be assessed.  it is un-
!             changed on output unless a previous step achieved a
!             better objective function reduction, in which case stlstg
!             will have been copied to step.
! stlstg (i/o) when assess recommends recomputing step even though the
!             current (or a previous) step yields an objective func-
!             tion decrease, it saves in stlstg the step that gave the
!             best function reduction seen so far (in the current itera-
!             tion).  if the recomputed step yields a larger function
!             value, then step is restored from stlstg and
!             x = x0 + step is recomputed.
!      v (i/o) real parameter and scratch vector -- see description
!             below of v values referenced.
!      x (i/o) on input, x = x0 + step is the point at which the objec-
!             tive function has just been evaluated.  if an earlier
!             step yielded a bigger function decrease, then x is
!             restored to the corresponding earlier value.  otherwise,
!             if the current step does not give any function decrease,
!             then x is restored to x0.
!     x0 (in)  initial objective function parameter vector (at the
!             start of the current iteration).

!  iv values referenced

!    iv(irc) (i/o) on input for the first step tried in a new iteration,
!             iv(irc) should be set to 3 or 4 (the value to which it is
!             set when step is definitely to be accepted).  on input
!             after step has been recomputed, iv(irc) should be
!             unchanged since the previous return of assess.
!                on output, iv(irc) is a return code having one of the
!             following values...
!                  1 = switch models or try smaller step.
!                  2 = switch models or accept step.
!                  3 = accept step and determine v(radfac) by gradient
!                       tests.
!                  4 = accept step, v(radfac) has been determined.
!                  5 = recompute step (using the same model).
!                  6 = recompute step with radius = v(lmax0) but do not
!                       evaluate the objective function.
!                  7 = x-convergence (see v(xctol)).
!                  8 = relative function convergence (see v(rfctol)).
!                  9 = both x- and relative function convergence.
!                 10 = absolute function convergence (see v(afctol)).
!                 11 = singular convergence (see v(lmax0)).
!                 12 = false convergence (see v(xftol)).
!                 13 = iv(irc) was out of range on input.
!             return code i has precedence over i+1 for i = 9, 10, 11.
! iv(mlstgd) (i/o) saved value of iv(model).
!  iv(model) (i/o) on input, iv(model) should be an integer identifying
!             the current quadratic model of the objective function.
!             if a previous step yielded a better function reduction,
!             then iv(model) will be set to iv(mlstgd) on output.
! iv(nfcall) (in)  invocation count for the objective function.
! iv(nfgcal) (i/o) value of iv(nfcall) at step that gave the biggest
!             function reduction this iteration.  iv(nfgcal) remains
!             unchanged until a function reduction is obtained.
! iv(radinc) (i/o) the number of radius increases (or minus the number
!             of decreases) so far this iteration.
! iv(restor) (out) set to 0 unless x and v(f) have been restored, in
!             which case assess sets iv(restor) = 1.
!  iv(stage) (i/o) count of the number of models tried so far in the
!             current iteration.
! iv(stglim) (in)  maximum number of models to consider.
! iv(switch) (out) set to 0 unless a new model is being tried and it
!             gives a smaller function value than the previous model,
!             in which case assess sets iv(switch) = 1.
! iv(toobig) (in)  is nonzero if step was too big (e.g. if it caused
!             overflow).
!   iv(xirc) (i/o) value that iv(irc) would have in the absence of
!             convergence, false convergence, and oversized steps.

!  ***  v values referenced  ***

! v(afctol) (in)  absolute function convergence tolerance.  if the
!             absolute value of the current function value v(f) is less
!             than v(afctol), then assess returns with iv(irc) = 10.
! v(decfac) (in)  factor by which to decrease radius when iv(toobig) is
!             nonzero.
! v(dstnrm) (in)  the 2-norm of d*step.
! v(dstsav) (i/o) value of v(dstnrm) on saved step.
!   v(dst0) (in)  the 2-norm of d times the newton step (when defined,
!             i.e., for v(nreduc)>=0).
!      v(f) (i/o) on both input and output, v(f) is the objective func-
!             tion value at x.  if x is restored to a previous value,
!             then v(f) is restored to the corresponding value.
!   v(fdif) (out) the function reduction v(f0) - v(f) (for the output
!             value of v(f) if an earlier step gave a bigger function
!             decrease, and for the input value of v(f) otherwise).
! v(flstgd) (i/o) saved value of v(f).
!     v(f0) (in)  objective function value at start of iteration.
! v(gtslst) (i/o) value of v(gtstep) on saved step.
! v(gtstep) (in)  inner product between step and gradient.
! v(incfac) (in)  minimum factor by which to increase radius.
!  v(lmax0) (in)  maximum reasonable step size (and initial step bound).
!             if the actual function decrease is no more than twice
!             what was predicted, if a return with iv(irc) = 7, 8, 9,
!             or 10 does not occur, if v(dstnrm)>v(lmax0), and if
!             v(preduc)<=v(rfctol) * abs(v(f0)), then assess re-
!             turns with iv(irc) = 11.  if so doing appears worthwhile,
!             then assess repeats this test with v(preduc) computed for
!             a step of length v(lmax0) (by a return with iv(irc) = 6).
! v(nreduc) (i/o)  function reduction predicted by quadratic model for
!             newton step.  if assess is called with iv(irc) = 6, i.e.,
!             if v(preduc) has been computed with radius = v(lmax0) for
!             use in the singular convergence test, then v(nreduc) is
!             set to -v(preduc) before the latter is restored.
! v(plstgd) (i/o) value of v(preduc) on saved step.
! v(preduc) (i/o) function reduction predicted by quadratic model for
!             current step.
! v(radfac) (out) factor to be used in determining the new radius,
!             which should be v(radfac)*dst, where  dst  is either the
!             output value of v(dstnrm) or the 2-norm of
!             diag(newd)*step  for the output value of step and the
!             updated version, newd, of the scale vector d.  for
!             iv(irc) = 3, v(radfac) = 1.0 is returned.
! v(rdfcmn) (in)  minimum value for v(radfac) in terms of the input
!             value of v(dstnrm) -- suggested value = 0.1.
! v(rdfcmx) (in)  maximum value for v(radfac) -- suggested value = 4.0.
!  v(reldx) (out) scaled relative change in x caused by step, computed
!             by function  reldst  as
!                 max (d(i)*abs(x(i)-x0(i)), 1<=i<=p) /
!                    max (d(i)*(abs(x(i))+abs(x0(i))), 1<=i<=p).
!             if an acceptable step is returned, then v(reldx) is com-
!             puted using the output (possibly restored) values of x
!             and step.  otherwise it is computed using the input
!             values.
! v(rfctol) (in)  relative function convergence tolerance.  if the
!             actual function reduction is at most twice what was pre-
!             dicted and  v(nreduc)<=v(rfctol)*abs(v(f0)),  then
!             assess returns with iv(irc) = 8 or 9.  see also v(lmax0).
! v(stppar) (in)  marquardt parameter -- 0 means full newton step.
! v(tuner1) (in)  tuning constant used to decide if the function
!             reduction was much less than expected.  suggested
!             value = 0.1.
! v(tuner2) (in)  tuning constant used to decide if the function
!             reduction was large enough to accept step.  suggested
!             value = 10**-4.
! v(tuner3) (in)  tuning constant used to decide if the radius
!             should be increased.  suggested value = 0.75.
!  v(xctol) (in)  x-convergence criterion.  if step is a newton step
!             (v(stppar) = 0) having v(reldx)<=v(xctol) and giving
!             at most twice the predicted function decrease, then
!             assess returns iv(irc) = 7 or 9.
!  v(xftol) (in)  false convergence tolerance.  if step gave no or only
!             a small function decrease and v(reldx)<=v(xftol),
!             then assess returns with iv(irc) = 12.

!-------------------------------  notes  -------------------------------

!  ***  application and usage restrictions  ***

!        this routine is called as part of the nl2sol (nonlinear
!     least-squares) package.  it may be used in any unconstrained
!     minimization solver that uses dogleg, goldfeld-quandt-trotter,
!     or levenberg-marquardt steps.

!  ***  algorithm notes  ***

!        see (1) for further discussion of the assessing and model
!     switching strategies.  while nl2sol considers only two models,
!     assess is designed to handle any number of models.

!  ***  usage notes  ***

!        on the first call of an iteration, only the i/o variables
!     step, x, iv(irc), iv(model), v(f), v(dstnrm), v(gtstep), and
!     v(preduc) need have been initialized.  between calls, no i/o
!     values except step, x, iv(model), v(f) and the stopping toler-
!     ances should be changed.
!        after a return for convergence or false convergence, one can
!     change the stopping tolerances and call assess again, in which
!     case the stopping tests will be repeated.

!  ***  references  ***

!     (1) dennis, j.e., jr., gay, d.m., and welsch, r.e. (1981),
!        an adaptive nonlinear least-squares algorithm,
!        acm trans. math. software, vol. 7, no. 3.

!     (2) powell, m.j.d. (1970)  a fortran subroutine for solving
!        systems of nonlinear algebraic equations, in numerical
!        methods for nonlinear algebraic equations, edited by
!        p. rabinowitz, gordon and breach, london.

!  ***  history  ***

!        john dennis designed much of this routine, starting with
!     ideas in (2). roy welsch suggested the model switching strategy.
!        david gay and stephen peters cast this subroutine into a more
!     portable form (winter 1977), and david gay cast it into its
!     present form (fall 1978).

!  ***  general  ***

!     this subroutine was written in connection with research
!     supported by the national science foundation under grants
!     mcs-7600324, dcr75-10143, 76-14311dss, mcs76-11989, and
!     mcs-7906671.

!------------------------  external quantities  ------------------------

!  ***  no common blocks  ***

!--------------------------  local variables  --------------------------

 LOGICAL goodx
 INTEGER i,nfc
 REAL(dp) emax,gts,reldx1,rfac1,xmax

! Data initializations

 REAL(dp),PARAMETER :: half=0.5d0,one=1.d0,two=2.d0,zero=0.d0

 INTEGER,PARAMETER :: irc=3,mlstgd=4,model=5,nfcall=6,nfgcal=7,radinc=8, &
  &restor=9,stage=10,stglim=11,switch=12,toobig=2,xirc=13

 INTEGER,PARAMETER :: afctol=31,decfac=22,dstnrm=2,dst0=3,dstsav=18,f=10,&
  &fdif=11,flstgd=12,f0=13,gtslst=14,gtstep=4,incfac=23,lmax0=35,nreduc=6, &
  &plstgd=15,preduc=7,radfac=16,rdfcmn=24,rdfcmx=25,reldx=17,rfctol=32,&
  &stppar=5,tuner1=26,tuner2=27,tuner3=28,xctol=33,xftol=34

 nfc=iv(nfcall)
 iv(switch)=0
 iv(restor)=0
 rfac1=one
 goodx=.true.
 i=iv(irc)
 select case (i)
 case (1)     ; goto 20
 case (2)     ; goto 30
 case (3:4)   ; goto 10
 case (5)     ; goto 40
 case (6)     ; goto 280
 case (7:11)  ; goto 220
 case (12)    ; goto 170
 case default ; iv(irc)=13 ; goto 300
 endselect

! Initialize for new iteration.

10 iv(stage)=1
 iv(radinc)=0
 v(flstgd)=v(f0)
 if(iv(toobig)==0)goto 90
 iv(stage)=-1
 iv(xirc)=i
 goto 60

! Step was recomputed with new model or smaller radius, first decide which

20 if(iv(model)/=iv(mlstgd))goto 30
! Old model retained, smaller radius tried
! Do not consider any more new models this iteration
 iv(stage)=iv(stglim)
 iv(radinc)=-1
 goto 90

! A new model is being tried.  decide whether to keep it.

30 iv(stage)=iv(stage)+1

! Now we add the possibility that step was recomputed with
! the same model, perhaps because of an oversized step.

40 if(iv(stage)>0)goto 50

! Step was recomputed because it was too big.

 if(iv(toobig)/=0)goto 60

! Restore iv(stage) and pick up where we left off.

 iv(stage)=-iv(stage)
 i=iv(xirc)
 select case(i)
 case(1)   ; goto 20
 case(2)   ; goto 30
 case(3:4) ; goto 90
 case(5)   ; goto 70
 endselect

50 if(iv(toobig)==0)goto 70

! Handle oversize step.

 if(iv(radinc)>0)goto 80
 iv(stage)=-iv(stage)
 iv(xirc)=iv(irc)

60 v(radfac)=v(decfac)
 iv(radinc)=iv(radinc)-1
 iv(irc)=5
 goto 300

70 if(v(f)<v(flstgd))goto 90

! The new step is a loser. Restore old model.

 if(iv(model)==iv(mlstgd))goto 80
 iv(model)=iv(mlstgd)
 iv(switch)=1

! Restore step, etc. only if a previous step decreased v(f).

80 if(v(flstgd)>=v(f0))goto 90
 iv(restor)=1
 v(f)=v(flstgd)
 v(preduc)=v(plstgd)
 v(gtstep)=v(gtslst)
 if(iv(switch)==0)rfac1=v(dstnrm)/v(dstsav)
 v(dstnrm)=v(dstsav)
 nfc=iv(nfgcal)
 goodx=.false.

! Compute relative change in x by current step.

90 reldx1=reldst(p,d,x,x0)

! Restore x and step if necessary.

 if(goodx)goto 110
 do i=1,p
  step(i)=stlstg(i)
  x(i)=x0(i)+stlstg(i)
 enddo

110 v(fdif)=v(f0)-v(f)
 if(v(fdif)>v(tuner2)*v(preduc))goto 140

! No (or only a trivial) function decrease
! -- so try new model or smaller radius.

 v(reldx)=reldx1
 if(v(f)<v(f0))goto 120
 iv(mlstgd)=iv(model)
 v(flstgd)=v(f)
 v(f)=v(f0)
 call vcopy(p,x(1:p),x0)
 iv(restor)=1
 goto 130
120 iv(nfgcal)=nfc
130 iv(irc)=1
 if(iv(stage)<iv(stglim))goto 160
 iv(irc)=5
 iv(radinc)=iv(radinc)-1
 goto 160

! Nontrivial function decrease achieved.

140 iv(nfgcal)=nfc
 rfac1=one
 if(goodx)v(reldx)=reldx1
 v(dstsav)=v(dstnrm)
 if(v(fdif)>v(preduc)*v(tuner1))goto 190

! Decrease was much less than predicted -- either change models
! or accept step with decreased radius.

 if(iv(stage)>=iv(stglim))goto 150
! Consider switching models.
 iv(irc)=2
 goto 160

! Accept step with decreased radius.

150 iv(irc)=4

! Set v(radfac) to Fletcher*s decrease factor.

160 iv(xirc)=iv(irc)
 emax=v(gtstep)+v(fdif)
 v(radfac)=half*rfac1
 if(emax<v(gtstep))v(radfac)=rfac1*max(v(rdfcmn),half*v(gtstep)/emax)

! Do false convergence test.

170 if(v(reldx)<=v(xftol))goto 180
 iv(irc)=iv(xirc)
 if(v(f)<v(f0))goto 200
 goto 230

180 iv(irc)=12
 goto 240

! Handle good function decrease.

190 if(v(fdif)<(-v(tuner3)*v(gtstep)))goto 210

! Increasing radius looks worthwhile.  See if we just
! recomputed step with a decreased radius or restored step
! after recomputing it with a larger radius.

 if(iv(radinc)<0)goto 210
 if(iv(restor)==1)goto 210

! We did not. Try a longer step unless this was a newton step.

 v(radfac)=v(rdfcmx)
 gts=v(gtstep)
 if(v(fdif)<(half/v(radfac)-one)*gts) v(radfac)=MAX(v(incfac),half*gts/ &
  &(gts+v(fdif)))
 iv(irc)=4
 if(v(stppar)==zero)goto 230
! Step was not a newton step. Recompute it with a larger radius.
 iv(irc)=5
 iv(radinc)=iv(radinc)+1

! Save values corresponding to good step.

200 v(flstgd)=v(f)
 iv(mlstgd)=iv(model)
 call vcopy(p,stlstg(1:p),step)
 v(dstsav)=v(dstnrm)
 iv(nfgcal)=nfc
 v(plstgd)=v(preduc)
 v(gtslst)=v(gtstep)
 goto 230

! Accept step with radius unchanged.

210 v(radfac)=one
 iv(irc)=3
 goto 230

! Come here for a restart after convergence.

220 iv(irc)=iv(xirc)
 if(v(dstsav)>=zero)goto 240
 iv(irc)=12
 goto 240

! Perform convergence tests.

230 iv(xirc)=iv(irc)
240 if(abs(v(f))<v(afctol))iv(irc)=10
 if(half*v(fdif)>v(preduc))goto 300
 emax=v(rfctol)*abs(v(f0))
 if(v(dstnrm)>v(lmax0).and.v(preduc)<=emax)iv(irc)=11
 if(v(dst0)<zero)goto 250
 i=0
 if((v(nreduc)>zero.and.v(nreduc)<=emax).or. &
  &(v(nreduc)==zero.and.v(preduc)==zero))i=2
 if(v(stppar)==zero.and.v(reldx)<=v(xctol).and.goodx)i=i+1
 if(i>0)iv(irc)=i+6

! Consider recomputing step of length v(lmax0) for singular convergence test.

250 if(abs(iv(irc)-3)>2.and.iv(irc)/=12)goto 300
 if(v(dstnrm)>v(lmax0))goto 260
 if(v(preduc)>=emax)goto 300
 if(v(dst0)<=zero)goto 270
 if(half*v(dst0)<=v(lmax0)) goto 300
 goto 270
260 if(half*v(dstnrm)<=v(lmax0))goto 300
 xmax=v(lmax0)/v(dstnrm)
 if(xmax*(two-xmax)*v(preduc)>=emax)goto 300
270 if(v(nreduc)<zero)goto 290

! Recompute v(preduc) for use in singular convergence test.

 v(gtslst)=v(gtstep)
 v(dstsav)=v(dstnrm)
 if(iv(irc)==12)v(dstsav)=-v(dstsav)
 v(plstgd)=v(preduc)
 iv(irc)=6
 call vcopy(p,stlstg(1:p),step)
 goto 300

! Perform singular convergence test with recomputed v(preduc).

280 v(gtstep)=v(gtslst)
 v(dstnrm)=abs(v(dstsav))
 call vcopy(p,step(1:p),stlstg)
 iv(irc)=iv(xirc)
 if(v(dstsav)<=zero)iv(irc)=12
 v(nreduc)=-v(preduc)
 v(preduc)=v(plstgd)
290 if(-v(nreduc)<=v(rfctol)*abs(v(f0)))iv(irc)=11

300 return

END SUBROUTINE assess


SUBROUTINE covclc(covirc,d,iv,j,n,p,r,v,x)
!----------------------------------------------------------------------!
! Ccompute covariance matrix for nl2itr (NL2SOL version 2.2).          !
!                                                                      !
!   let k = ABS(iv(covreq).  for k <= 2, a finite-difference           !
!   Hessian h is computed (using func. and grad. values if             !
!   iv(covreq) is nonnegative, and using only func. values if          !
!   iv(covreq) is negative).  for scale = 2*f(x) / max(1, n-p),        !
!   where 2*f(x) is the residual sum of squares, covclc computes...    !
!              k = 0 or 1...  scale * h**-1 * (j**t * j) * h**-1.      !
!              k = 2...  scale * h**-1.                                !
!              k >= 3...  scale * (j**t * j)**-1.                      !
!----------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: n, p
 INTEGER,INTENT(inout) :: covirc,iv(:)
 REAL(dp),INTENT(in) :: d(:)
 REAL(dp),INTENT(inout) :: v(:),x(:),j(:,:)
 REAL(dp),INTENT(inout) :: r(:)
 INTEGER gp,gsave1,g1,hc,hmi,hpi,hpm,i,ipivi,ipivk,ip1,irc,k,&
  &kind,kl,l,m,mm1,mm1o2,pp1o2,qtr1,rd1,stpi,stpm,stp0,wl,w0,w1
 INTEGER,SAVE :: cov
 REAL(dp) del,t,wk
 LOGICAL havej

! linvrt... invert lower triangular matrix.
! litvmu... apply inverse-transpose of compact lower triang. matrix.
! livmul... apply inverse of compact lower triang. matrix.
! lsqrt.... compute cholesky factor of (lower trinag. of) a sym. matrix.
! ltsqar... given lower triang. matrix l, compute (l**t)*l.
! qrfact... compute qr decomposition of a matrix.
! vcopy.... copy one vector to another.
! vscopy... set all elements of a vector to a scalar.

! Subscripts for iv and v.

 REAL(dp),PARAMETER :: half=0.5d0,negpt5=-0.5d0,one=1.d0,two=2.d0,zero=0.d0

 INTEGER,PARAMETER :: covmat=26,covreq=15,delta=50,delta0=44,dltfdc=40,f=10,&
  &fx=46,g=28,h=44,ierr=32,ipivot=61,ipiv0=60,kagqt=35,kalm=36,lmat=58,mode=38,&
  &nfgcal=7,qtr=49,rd=51,rsave=52,savei=54,switch=12,toobig=2,w=59,xmsave=49

 covirc=4
 kind=iv(covreq)
 m=iv(mode)
 if(m>0)goto 10
 iv(kagqt)=-1
 if(iv(kalm)>0)iv(kalm)=0
 if(abs(kind)>=3)goto 310
 v(fx)=v(f)
 k=iv(rsave)
 call vcopy (n,v(k:k+n-1),r)
10 if(m>p)goto 220
 if(kind<0)goto 110

! Compute finite-difference Hessian using both function and  gradient values.

 gsave1=iv(w)+p
 g1=iv(g)
 if(m>0)goto 20
! First call on covclc.  set gsave = g, take first step  ***
 call vcopy(p,v(gsave1:gsave1+p-1),v(g1:))
 iv(switch)=iv(nfgcal)
 goto 90

20 del=v(delta)
 x(m)=v(xmsave)
 if(iv(toobig)==0)goto 40

! Handle oversize v(delta).

 if(del*x(m)>zero)goto 30
! We already tried shrinking v(delta), so quit.
 iv(covmat)=-2
 goto 210

! Try shrinking v(delta).
30 del=negpt5*del
 goto 100

40 cov=iv(lmat)
 gp=g1+p-1

! Set g = (g - gsave)/del.

 do i=g1,gp
  v(i)=(v(i)-v(gsave1))/del
  gsave1=gsave1+1
 enddo

! Add g as new col. to finite-diff. Hessian matrix.

 k=cov+m*(m-1)/2
 l=k+m-2
 if(m==1)goto 70

! Set h(i,m) = 0.5 * (h(i,m) + g(i))  for i = 1 to m-1

 do i=k,l
  v(i)=half*(v(i)+v(g1))
  g1=g1+1
 enddo

! Add h(i,m) = g(i)  for i = m to p

70 l=l+1
 do i=m,p
  v(l)=v(g1)
  l=l+i
  g1=g1+1
 enddo

90 m=m+1
 iv(mode)=m
 if(m>p)goto 210

! Choose next finite-difference step, return to get g there.

 del=v(delta0)*max(one/d(m),abs(x(m)))
 if(x(m)<zero)del=-del
 v(xmsave)=x(m)
100 x(m)=x(m)+del
 v(delta)=del
 covirc=2
 goto 390

! Compute finite-difference Hessian using function values only.

110 stp0=iv(w)+p-1
 mm1=m-1
 mm1o2=m*mm1/2
 if(m>0)goto 120
! First call on covclc.
 iv(savei)=0
 goto 200

120 i=iv(savei)
 if(i>0)goto 180
 if(iv(toobig)==0)goto 140

! Handle oversize step.

 stpm=stp0+m
 del=v(stpm)
 if(del*x(xmsave)>zero)goto 130
! We already tried shrinking the step, so quit.
 iv(covmat)=-2
 goto 390

! Try shrinking the step.
130 del=negpt5*del
 x(m)=x(xmsave)+del
 v(stpm)=del
 covirc=1
 goto 390

! Save f(x + stp(m)*e(m)) in h(p,m).

140 pp1o2=p*(p-1)/2
 cov=iv(lmat)
 hpm=cov+pp1o2+mm1
 v(hpm)=v(f)

! Start computing row m of the finite-difference Hessian h.

 hmi=cov+mm1o2
 if(mm1==0)goto 160
 hpi=cov+pp1o2
 do i=1,mm1
  v(hmi)=v(fx)-(v(f)+v(hpi))
  hmi=hmi+1
  hpi=hpi+1
 enddo
160 v(hmi)=v(f)-two*v(fx)

! Compute function values needed to complete row m of h.

 i=1

170 iv(savei)=i
 stpi=stp0+i
 v(delta)=x(i)
 x(i)=x(i)+v(stpi)
 if(i==m)x(i)=v(xmsave)-v(stpi)
 covirc=1
 goto 390

180 x(i)=v(delta)
 if(iv(toobig)==0)goto 190
! Punt in the event of an oversize step.
 iv(covmat)=-2
 goto 390

! Finish computing h(m,i).

190 stpi=stp0+i
 hmi=cov+mm1o2+i-1
 stpm=stp0+m
 v(hmi)=(v(hmi)+v(f))/(v(stpi)*v(stpm))
 i=i+1
 if(i<=m)goto 170
 iv(savei)=0
 x(m)=v(xmsave)

200 m=m+1
 iv(mode)=m
 if(m>p)goto 210

! Prepare to compute row m of the finite-difference Hessian h.
! Compute m-th step size stp(m), then return to obtain
! f(x + stp(m)*e(m)), where e(m) = m-th std. unit vector.

 del=v(dltfdc)*MAX(one/d(m),ABS(x(m)))
 if(x(m)<zero)del=-del
 v(xmsave)=x(m)
 x(m)=x(m)+del
 stpm=stp0+m
 v(stpm)=del
 covirc=1
 goto 390

! Restore r, v(f), etc.

210 k=iv(rsave)
 call vcopy(n,r(1:n),v(k:))
 v(f)=v(fx)
 if(kind<0)goto 220
 iv(nfgcal)=iv(switch)
 qtr1=iv(qtr)
 call vcopy(n,v(qtr1:qtr1+n-1),r)
 if(iv(covmat)<0)goto 390
 covirc=3
 goto 390

220 cov=iv(lmat)

! The complete finite-diff. Hessian is now stored at v(cov).
! Use it to compute the requested covariance matrix.

! Compute Cholesky factor c of h = c*(c**t) and store it at v(hc).

 hc=cov
 if(abs(kind)==2)goto 230
 hc=abs(iv(h))
 iv(h)=-hc
230 call lsqrt(1,p,v(hc:hc+(p*(p+1)/2)-1),v(cov:),irc)
 iv(covmat)=-1
 if(irc/=0)goto 390

 w1=iv(w)+p
 if(abs(kind)>1)goto 340

! Covariance = scale * h**-1 * (j**t * j) * h**-1

 call vscopy(p*(p+1)/2,v(cov:cov+p*(p+1)/2-1),zero)
 havej=iv(kalm)==(-1)
! havej = .true. means j is in its original form, while
! havej = .false. means qrfact has been applied to j.

 m=p
 if(havej)m=n
 w0=w1-1
 rd1=iv(rd)
 do i=1,m
  if(havej) goto 250

! Set w = ipivot * (row i of r matrix from qrfact).

  call vscopy(p,v(w1:w1+p-1),zero)
  ipivi=ipiv0+i
  l=w0+iv(ipivi)
  v(l)=v(rd1)
  rd1=rd1+1
  if(i==p)goto 270
  ip1=i+1
  do k=ip1,p
   ipivk=ipiv0+k
   l=w0+iv(ipivk)
   v(l)=j(i,k)
  enddo
  goto 270

! Set w = (row i of j).

250 l=w0
  do k=1,p
   l=l+1
   v(l)=j(i,k)
  enddo

! Set w = h**-1 * w.

270 call livmul(p,v(w1:w1+p-1),v(hc:),v(w1:))
  call litvmu(p,v(w1:w1+p-1),v(hc:),v(w1:))

! Add w * w**t to covariance matrix.

  kl=cov
  do k=1,p
   l=w0+k
   wk=v(l)
   do l=1,k
    wl=w0+l
    v(kl)=v(kl)+wk*v(wl)
    kl=kl+1
   enddo
  enddo
 enddo
 goto 370

! Covariance = scale * (j**t * j)**-1.

310 rd1=iv(rd)
 if(iv(kalm)/=(-1))goto 320

! Apply qrfact to j.

 qtr1=iv(qtr)
 call vcopy(n,v(qtr1:qtr1+n-1),r)
 w1=iv(w)+p
 call qrfact(n,p,j,v(rd1:rd1+p-1),iv(ipivot:ipivot+p-1),iv(ierr),0,v(w1:w1+p-1))
 iv(kalm)=-2
320 iv(covmat)=-1
 if(iv(ierr)/=0)goto 390
 cov=iv(lmat)
 hc=abs(iv(h))
 iv(h)=-hc

! Set hc = (r matrix from qrfact).

 l=hc
 do i=1,p
  if(i>1)call vcopy(i-1,v(l:l+i-1-1),j(:,i))
  l=l+i-1
  v(l)=v(rd1)
  l=l+1
  rd1=rd1+1
 enddo

! The Cholesky factor c of the unscaled inverse covariance matrix
! (or permutation thereof) is stored at v(hc).

! set c = c**-1.

340 call linvrt(p,v(hc:hc+p*(p+1)/2-1),v(hc:))

! set c = c**t * c.

 call ltsqar(p,v(hc:hc+p*(p+1)/2-1),v(hc:))

 if(hc==cov)goto 370

! c = permuted, unscaled covariance.
! set cov = ipivot * c * ipivot**t.

 do i=1,p
  m=ipiv0+i
  ipivi=iv(m)
  kl=cov-1+ipivi*(ipivi-1)/2
  do k=1,i
   m=ipiv0+k
   ipivk=iv(m)
   l=kl+ipivk
   if(ipivk>ipivi)l=l+(ipivk-ipivi)*(ipivk+ipivi-3)/2
   v(l)=v(hc)
   hc=hc+1
  enddo
 enddo

370 iv(covmat)=cov

! Apply scale factor = (resid. sum of squares) / max(1,n-p).

 t=v(f)/(half*real(max(1,n-p)))
 k=cov-1+p*(p+1)/2
 v(cov:k)=t*v(cov:k)

390 return

END SUBROUTINE covclc


SUBROUTINE dfault(iv,v)
!-----------------------------------------------------------!
! Supply NL2SOL (version 2.2) default values to iv and v.   !
!-----------------------------------------------------------!
INTEGER,INTENT(inout) :: iv(:)
REAL(dp),INTENT(inout) :: v(:)
REAL(dp) machep,mepcrt,sqteps
REAL(dp),PARAMETER :: one=1.d0,three=3.d0
! iv subscript values
INTEGER,PARAMETER :: covprt=14,covreq=15,dtype=16,inits=25,mxfcal=17,&
 &mxiter=18,outlev=19,parprt=20,prunit=21,solprt=22,statpr=23,x0prt=24
! v subscript values
INTEGER,PARAMETER :: afctol=31,cosmin=43,decfac=22,delta0=44,dfac=41,&
 &dinit=38,dltfdc=40,dltfdj=36,d0init=37,epslon=19,fuzz=45,incfac=23,jtinit=39,&
 &lmax0=35,phmnfc=20,phmxfc=21,rdfcmn=24,rdfcmx=25,rfctol=32,rlimit=42,&
 &tuner1=26,tuner2=27,tuner3=28,tuner4=29,tuner5=30,xctol=33,xftol=34

iv(1)=12
iv(covprt)=1
iv(covreq)=1
iv(dtype)=1
iv(inits)=0
iv(mxfcal)=200
iv(mxiter)=150
iv(outlev)=1
iv(parprt)=1
iv(prunit)=imdcon(1)
iv(solprt)=1
iv(statpr)=1
iv(x0prt)=1

machep=rmdcon(3)
v(afctol)=1.d-20
if(machep>1.d-10)v(afctol)=machep**2
v(cosmin)=MAX(1.d-6,1.d2*machep)
v(decfac)=0.5d0
sqteps=rmdcon(4)
v(delta0)=sqteps
v(dfac)=0.6d0
v(dinit)=0.d0
mepcrt=machep**(one/three)
v(dltfdc)=mepcrt
v(dltfdj)=sqteps
v(d0init)=1.d0
v(epslon)=0.1d0
v(fuzz)=1.5d0
v(incfac)=2.d0
v(jtinit)=1.d-6
v(lmax0)=100.d0
v(phmnfc)=-0.1d0
v(phmxfc)=0.1d0
v(rdfcmn)=0.1d0
v(rdfcmx)=4.d0
v(rfctol)=MAX(1.d-10,mepcrt**2)
v(rlimit)=rmdcon(5)
v(tuner1)=0.1d0
v(tuner2)=1.d-4
v(tuner3)=0.75d0
v(tuner4)=0.5d0
v(tuner5)=0.75d0
v(xctol)=sqteps
v(xftol)=1.d2*machep

END SUBROUTINE dfault


FUNCTION dotprd(p,x,y) RESULT(fn_val)
!--------------------------------------------------------!
! Return the inner product of the p-vectors x and y.     !
!--------------------------------------------------------!
INTEGER,INTENT(in) :: p
REAL(dp),INTENT(in) :: x(:),y(:)
REAL(dp) fn_val
INTEGER i
REAL(dp) t
! rmdcon(2) returns a machine-dependent constant, sqteta, which
! is slightly larger than the smallest positive number that
! can be squared without underflowing.
REAL(dp),PARAMETER :: one=1.d0,zero=0.d0
REAL(dp),SAVE :: sqteta=0.d0

fn_val=zero
if(p<=0)goto 30
if(sqteta==zero)sqteta=rmdcon(2)
do i=1,p
 t=max(abs(x(i)),abs(y(i)))
 if(t>one)goto 10
 if(t<sqteta)cycle
 t=(x(i)/sqteta)*y(i)
 if(ABS(t)<sqteta)cycle
 10 fn_val=fn_val+x(i)*y(i)
enddo

30 return
END FUNCTION dotprd


SUBROUTINE dupdat(d,iv,j,n,p,v)
!---------------------------------------------------------!
! Update scale vector d for nl2itr (nl2sol version 2.2)   !
!---------------------------------------------------------!
INTEGER,INTENT(in)   :: iv(:),n,p
REAL(dp),INTENT(in)  :: j(:,:),v(:)
REAL(dp),INTENT(inout) :: d(:)
! dimension iv(*), v(*), j(nn,p)
INTEGER d0,i,jtoli,s1
REAL(dp) sii,t,vdfac
REAL(dp) :: zero=0.d0
! Subscripts for iv and v
INTEGER,PARAMETER :: dfac=41,dtype=16,jtol0=86,niter=31,s=53

i=iv(dtype)
if(i==1)goto 10
if(iv(niter)>0)goto 30

10 vdfac=v(dfac)
d0=jtol0+p
s1=iv(s)-1
do i=1,p
 s1=s1+i
 sii=v(s1)
 t=v2norm(n,j(:,i))
 if(sii>zero)t=sqrt(t*t+sii)
 jtoli=jtol0+i
 d0=d0+1
 if(t<v(jtoli))t=MAX(v(d0),v(jtoli))
 d(i)=max(vdfac*d(i),t)
enddo

30 return

END SUBROUTINE dupdat


SUBROUTINE gqtstp(d,dig,dihdi,ka,l,p,step,v,w)
!-----------------------------------------------------------------------------!
! Compute Goldfeld-Quandt-Trotter step by More-Hebden technique.
! (NL2SOL version 2.2)
!
!  ***  purpose  ***
!
!        Given the (compactly stored) lower triangle of a scaled
!     Hessian (approximation) and a nonzero scaled gradient vector,
!     this subroutine computes a Goldfeld-Quandt-Trotter step of
!     approximate length v(radius) by the More-Hebden technique.  In
!     other words, step is computed to (approximately) minimize
!     psi(step) = (g**t)*step + 0.5*(step**t)*h*step  such that the
!     2-norm of d*step is at most (approximately) v(radius), where
!     g  is the gradient,  h  is the Hessian, and  d  is a diagonal
!     scale matrix whose diagonal is stored in the parameter d.
!     (gqtstp assumes  dig = d**-1 * g  and  dihdi = d**-1 * h * d**-1.)
!     if g = 0, however, step = 0 is returned (even at a saddle point).
!
!  ***  parameter description  ***
!
!     d (in)  = the scale vector, i.e. the diagonal of the scale
!              matrix  d  mentioned above under purpose.
!   dig (in)  = the scaled gradient vector, d**-1 * g.  if g = 0, then
!              step = 0 and  v(stppar) = 0 are returned.
! dihdi (in)  = lower triangle of the scaled Hessian (approximation),
!              i.e., d**-1 * h * d**-1, stored compactly by rows., i.e.,
!              in the order (1,1), (2,1), (2,2), (3,1), (3,2), etc.
!    ka (i/o) = the number of hebden iterations (so far) taken to deter-
!              mine step.  ka<0 on input means this is the first
!              attempt to determine step (for the present dig and dihdi)
!              -- ka is initialized to 0 in this case.  output with
!              ka = 0 (or v(stppar) = 0)  means  step = -(h**-1)*g.
!     l (i/o) = workspace of length p*(p+1)/2 for cholesky factors.
!     p (in)  = number of parameters -- the Hessian is a  p x p  matrix.
!  step (i/o) = the step computed.
!     v (i/o) contains various constants and variables described below.
!     w (i/o) = workspace of length 4*p + 6.
!
!  ***  entries in v  ***
!
! v(dgnorm) (i/o) = 2-norm of (d**-1)*g.
! v(dstnrm) (output) = 2-norm of d*step.
! v(dst0)   (i/o) = 2-norm of d*(h**-1)*g (for pos. def. h only), or
!             overestimate of smallest eigenvalue of (d**-1)*h*(d**-1).
! v(epslon) (in)  = max. rel. error allowed for psi(step).  for the
!             step returned, psi(step) will exceed its optimal value
!             by less than -v(epslon)*psi(step).  suggested value = 0.1.
! v(gtstep) (out) = inner product between g and step.
! v(nreduc) (out) = psi(-(h**-1)*g) = psi(newton step)  (for pos. def.
!             h only -- v(nreduc) is set to zero otherwise).
! v(phmnfc) (in)  = tol. (together with v(phmxfc)) for accepting step
!             (more*s sigma).  the error v(dstnrm) - v(radius) must lie
!             between v(phmnfc)*v(radius) and v(phmxfc)*v(radius).
! v(phmxfc) (in)  (see v(phmnfc).)
!             suggested values -- v(phmnfc) = -0.25, v(phmxfc) = 0.5.
! v(preduc) (out) = psi(step) = predicted obj. func. reduction for step.
! v(radius) (in)  = radius of current (scaled) trust region.
! v(rad0)   (i/o) = value of v(radius) from previous call.
! v(stppar) (i/o) is normally the marquardt parameter, i.e. the alpha
!             described below under algorithm notes.  if h + alpha*d**2
!             (see algorithm notes) is (nearly) singular, however,
!             then v(stppar) = -alpha.
!
!  ***  usage notes  ***
!
!     if it is desired to recompute step using a different value of
!     v(radius), then this routine may be restarted by calling it
!     with all parameters unchanged except v(radius).  (this explains
!     why step and w are listed as i/o).  on an initial call (one with
!     ka<0), step and w need not be initialized and only compo-
!     nents v(epslon), v(stppar), v(phmnfc), v(phmxfc), v(radius), and
!     v(rad0) of v must be initialized.  to compute step from a saddle
!     point (where the true gradient vanishes and h has a negative
!     eigenvalue), a nonzero g with small components should be passed.
!
!  ***  application and usage restrictions  ***
!
!     this routine is called as part of the nl2sol (nonlinear least-
!     squares) package (ref. 1), but it could be used in solving any
!     unconstrained minimization problem.
!
!  ***  algorithm notes  ***
!
!        the desired g-q-t step (ref. 2, 3, 4) satisfies
!     (h + alpha*d**2)*step = -g  for some nonnegative alpha such that
!     h + alpha*d**2 is positive semidefinite.  alpha and step are
!     computed by a scheme analogous to the one described in ref. 5.
!     estimates of the smallest and largest eigenvalues of the Hessian
!     are obtained from the gerschgorin circle theorem enhanced by a
!     simple form of the scaling described in ref. 6.  cases in which
!     h + alpha*d**2 is nearly (or exactly) singular are handled by
!     the technique discussed in ref. 2.  in these cases, a step of
!     (exact) length v(radius) is returned for which psi(step) exceeds
!     its optimal value by less than -v(epslon)*psi(step).
!
!  ***  functions and subroutines called  ***
!
! dotprd - returns inner product of two vectors.
! litvmu - applies inverse-transpose of compact lower triang. matrix.
! livmul - applies inverse of compact lower triang. matrix.
! lsqrt  - finds cholesky factor (of compactly stored lower triang.).
! lsvmin - returns approx. to min. sing. value of lower triang. matrix.
! rmdcon - returns machine-dependent constants.
! v2norm - returns 2-norm of a vector.
!
!  ***  references  ***
!
! 1.  dennis, j.e., gay, d.m., and welsch, r.e. (1981), an adaptive
!             nonlinear least-squares algorithm, acm trans. math.
!             software, vol. 7, no. 3.
! 2.  gay, d.m. (1981), computing optimal locally constrained steps,
!             siam j. sci. statist. computing, vol. 2, no. 2, pp.
!             186-197.
! 3.  goldfeld, s.m., quandt, r.e., and trotter, h.f. (1966),
!             maximization by quadratic hill-climbing, econometrica 34,
!             pp. 541-551.
! 4.  hebden, m.d. (1973), an algorithm for minimization using exact
!             second derivatives, report t.p. 515, theoretical physics
!             div., a.e.r.e. harwell, oxon., england.
! 5.  more, j.j. (1978), the levenberg-marquardt algorithm, implemen-
!             tation and theory, pp.105-116 of springer lecture notes
!             in mathematics no. 630, edited by g.a. watson, springer-
!             verlag, berlin and new york.
! 6.  varga, r.s. (1965), minimal gerschgorin sets, pacific j. math. 15,
!             pp. 719-729.
!
!  ***  general  ***
!
!     coded by david m. gay.
!     this subroutine was written in connection with research
!     supported by the national science foundation under grants
!     mcs-7600324, dcr75-10143, 76-14311dss, mcs76-11989, and
!     mcs-7906671.
!
!-----------------------------------------------------------------------------!
INTEGER,INTENT(in) :: p
INTEGER,INTENT(inout) :: ka
REAL(dp),INTENT(in) :: d(:),dig(:)
REAL(dp),INTENT(inout) :: l(:),v(:),step(:),w(:),dihdi(:)
INTEGER dggdmx,diag,diag0,dstsav,emax,emin,i,im1,inc,irc,j,k,kalim,k1,lk0, &
 &phipin,q,q0,uk0,x,x0
REAL(dp) alphak,aki,akk,delta,dst,epso6,lk,oldphi,phi,phimax,phimin,psifac,&
 &rad,root,si,sk,sw,t,twopsi,t1,uk,wi
LOGICAL restrt

! Subscripts for v
INTEGER,PARAMETER :: dgnorm=1,dstnrm=2,dst0=3,epslon=19,gtstep=4,   &
 &nreduc=6,phmnfc=20,phmxfc=21,preduc=7,radius=8,rad0=9,stppar=5

REAL(dp),PARAMETER :: epsfac=50.d0,four=4.d0,half=0.5d0,kappa=2.d0, &
 &negone=-1.d0,one=1.d0,p001=1.d-3,six=6.d0,three=3.d0,two=2.d0,zero=0.d0
! Save dgxfac
REAL(dp),SAVE :: dgxfac=0.d0

dst=0.d0 ; phi=0.d0 ; alphak=0.d0 ; uk=0.d0

! Store largest abs. entry in (d**-1)*h*(d**-1) at w(dggdmx).
dggdmx=p+1
! Store gerschgorin over- and underestimates of the largest
! and smallest eigenvalues of (d**-1)*h*(d**-1) at w(emax)
! and w(emin) respectively.
emax=dggdmx+1
emin=emax+1
! For use in recomputing step, the final values of lk, uk, dst,
! and the inverse derivative of more*s phi at 0 (for pos. def.
! h) are stored in w(lk0), w(uk0), w(dstsav), and w(phipin) respectively.
lk0=emin+1
phipin=lk0+1
uk0=phipin+1
dstsav=uk0+1
! Store diag of (d**-1)*h*(d**-1) in w(diag),...,w(diag0+p).
diag0=dstsav
diag=diag0+1
! Store -d*step in w(q),...,w(q0+p).
q0=diag0+p
q=q0+1
rad=v(radius)
! phitol = max. error allowed in dst = v(dstnrm) = 2-norm of d*step.
phimax=v(phmxfc)*rad
phimin=v(phmnfc)*rad
! epso6 and psifac are used in checking for the special case
! of (nearly) singular h + alpha*d**2 (see ref. 2).
psifac=two*v(epslon)/(three*(four*(v(phmnfc)+one)*(kappa+one)+&
 &kappa+two)*rad**2)
! oldphi is used to detect limits of numerical accuracy.
! If we recompute step and it does not change, then we accept it.
oldphi=zero
epso6=v(epslon)/six
irc=0
restrt=.false.
kalim=ka+50

! Start or restart, depending on ka

if(ka>=0)goto 290

! Fresh start

k=0
uk=negone
ka=0
kalim=50

! Store diag(dihdi) in w(diag0+1),...,w(diag0+p)

j=0
do i=1,p
 j=j+i
 k1=diag0+i
 w(k1)=dihdi(j)
enddo

! Determine w(dggdmx), the largest element of dihdi

t1=zero
j=p*(p+1)/2
do i=1,j
 t=abs(dihdi(i))
 if(t1<t)t1=t
enddo
w(dggdmx)=t1

! Try alpha = 0

30 call lsqrt (1, p, l(1:(p*(p+1)/2)), dihdi, irc)
if(irc==0)goto 50
! indef. h -- underestimate smallest eigenvalue, use this
! estimate to initialize lower bound lk on alpha.
j=irc*(irc+1)/2
t=l(j)
l(j)=one
w(1:irc) = zero
w(irc)=one
call litvmu(irc,w(1:irc),l,w)
t1=v2norm(irc,w)
lk=-t/t1/t1
v(dst0)=-lk
if(restrt)goto 200
v(nreduc)=zero
goto 60

! Positive definite h -- compute unmodified newton step.
50 lk=zero
call livmul(p,w(q:q+p-1),l,dig)
v(nreduc)=half*dotprd(p,w(q:),w(q:))
call litvmu(p,w(q:q+p-1),l,w(q:))
dst=v2norm(p,w(q:))
v(dst0)=dst
phi=dst-rad
if(phi<=phimax)goto 260
if(restrt)goto 200

! Prepare to compute gerschgorin estimates of largest (and
! smallest) eigenvalues.

60 v(dgnorm)=v2norm(p,dig)
if(v(dgnorm)==zero)goto 430
k=0
do i=1,p
 wi=zero
 if(i==1)goto 80
 im1=i-1
 do j=1,im1
  k=k+1
  t=abs(dihdi(k))
  wi=wi+t
  w(j)=w(j)+t
 enddo
 80 w(i)=wi
 k=k+1
enddo

! (under-)estimate smallest eigenvalue of (d**-1)*h*(d**-1)

k=1
t1=w(diag)-w(1)
if(p<=1)goto 110
do i=2,p
 j=diag0+i
 t=w(j)-w(i)
 if(t>=t1)cycle
 t1=t
 k=i
enddo

110 sk=w(k)
j=diag0+k
akk=w(j)
k1=k*(k-1)/2+1
inc=1
t=zero
do i=1,p
 if(i==k)goto 120
 aki=abs(dihdi(k1))
 si=w(i)
 j=diag0+i
 t1=half*(akk-w(j)+si-aki)
 t1=t1+sqrt(t1*t1+sk*aki)
 if(t<t1)t=t1
 if(i<k)goto 130
 120 inc=i
 130 k1=k1+inc
enddo

w(emin)=akk-t
uk=v(dgnorm)/rad-w(emin)

! Compute gerschgorin (over-)estimate of largest eigenvalue

k=1
t1=w(diag)+w(1)
if(p<=1)goto 160
do i=2,p
 j=diag0+i
 t=w(j)+w(i)
 if(t<=t1)cycle
 t1=t
 k=i
enddo

160 sk=w(k)
j=diag0+k
akk=w(j)
k1=k*(k-1)/2+1
inc=1
t=zero
do i=1,p
 if(i==k)goto 170
 aki=abs(dihdi(k1))
 si=w(i)
 j=diag0+i
 t1=half*(w(j)+si-aki-akk)
 t1=t1+sqrt(t1*t1+sk*aki)
 if(t<t1)t=t1
 if(i<k)goto 180
 170 inc=i
 180 k1=k1+inc
enddo

w(emax)=akk+t
lk=max(lk,v(dgnorm)/rad-w(emax))

! alphak = current value of alpha (see alg. notes above).  we
! use more*s scheme for initializing it.
alphak=abs(v(stppar))*v(rad0)/rad

if(irc/=0)goto 200

! Compute l0 for positive definite h

call livmul(p,w(1:p),l,w(q:))
t=v2norm(p,w)
w(phipin)=dst/t/t
lk=max(lk,phi*w(phipin))

! Safeguard alphak and add alphak*i to (d**-1)*h*(d**-1)

200 ka=ka+1
if(-v(dst0)>=alphak.or.alphak<lk.or.alphak>=uk)alphak=uk*max(p001,sqrt(lk/uk))
k=0
do i=1,p
 k=k+i
 j=diag0+i
 dihdi(k)=w(j)+alphak
enddo

! Try computing cholesky decomposition

call lsqrt(1,p,l(1:p*(p+1)/2),dihdi,irc)
if(irc==0)goto 230

! (d**-1)*h*(d**-1) + alphak*i  is indefinite -- overestimate
! smallest eigenvalue for use in updating lk

j=(irc*(irc+1))/2
t=l(j)
l(j)=one
do i=1,irc
 w(i)=zero
enddo
w(irc)=one
call litvmu(irc,w(1:irc),l,w)
t1=v2norm(irc,w)
lk=alphak-t/t1/t1
v(dst0)=-lk
goto 200

! alphak makes (d**-1)*h*(d**-1) positive definite.
! compute q = -d*step, check for convergence.

230 call livmul(p,w(q:q+p-1),l,dig)
call litvmu(p,w(q:q+p-1),l,w(q:))
dst=v2norm(p,w(q:))
phi=dst-rad
if(phi<=phimax.and.phi>=phimin)goto 270
if(phi==oldphi)goto 270
oldphi=phi
if(phi>zero)goto 240
! Check for the special case of  h + alpha*d**2  (nearly)
! singular.  delta is >= the smallest eigenvalue of
! (d**-1)*h*(d**-1) + alphak*i.
if(v(dst0)>zero)goto 240
delta=alphak+v(dst0)
twopsi=alphak*dst*dst+dotprd(p,dig,w(q:))
if(delta<psifac*twopsi)goto 250

! Unacceptable alphak -- update lk, uk, alphak

240 if(ka>=kalim)goto 270
call livmul(p,w(1:p),l,w(q:))
t1=v2norm(p,w)
! The following MIN is necessary because of restarts
if(phi<zero)uk=min(uk,alphak)
alphak=alphak+(phi/t1)*(dst/t1)*(dst/rad)
lk=max(lk,alphak)
goto 200

! Decide how to handle (nearly) singular h + alpha*d**2

! If not yet available, obtain machine dependent value dgxfac.
250 if(dgxfac==zero)dgxfac=epsfac*rmdcon(3)

! Now decide.
if(delta>dgxfac*w(dggdmx))goto 330
! delta is so small we cannot handle the special case in
! the available arithmetic.  accept step as it is.
goto 270

! Acceptable step on first try

260 alphak=zero

! Successful step in general.  compute step = -(d**-1)*q

270 step(1:p)=-w(q0+1:q0+p)/d(1:p)
v(gtstep)=-dotprd(p,dig,w(q:))
v(preduc)=half*(abs(alphak)*dst*dst-v(gtstep))
goto 410

! Restart with new radius

290 if(v(dst0)<=zero.or.v(dst0)-rad>phimax)goto 310

! Prepare to return Newton step

restrt=.true.
ka=ka+1
k=0
do i=1,p
 k=k+i
 j=diag0+i
 dihdi(k)=w(j)
enddo
uk=negone
goto 30

310 if(ka==0)goto 50

dst=w(dstsav)
alphak=abs(v(stppar))
phi=dst-rad
t=v(dgnorm)/rad
if(rad>v(rad0))goto 320

! Smaller radius
uk=t-w(emin)
lk=zero
if(alphak>zero)lk=w(lk0)
lk=max(lk,t-w(emax))
if(v(dst0)>zero)lk=max(lk,(v(dst0)-rad)*w(phipin))
goto 240

! Bigger radius
320 uk=t-w(emin)
if(alphak>zero)uk=min(uk,w(uk0))
lk=max(zero,-v(dst0),t-w(emax))
if(v(dst0)>zero)lk=max(lk,(v(dst0)-rad)*w(phipin))
goto 240

! Handle (nearly) singular h + alpha*d**2

! Negate alphak to indicate special case
330 alphak=-alphak
! Allocate storage for scratch vector x
x0=q0+p
x=x0+1

! Use inverse power method with start from lsvmin to obtain
! approximate eigenvector corresponding to smallest eigenvalue
! of (d**-1)*h*(d**-1).

delta=kappa*delta
call lsvmin(p,l,w(x:x+p-1),w(1:p),t)

k=0
! Normalize w
340 w(1:p)=t*w(1:p)
! Complete current inv. power iter. -- replace w by (l**-t)*w.
call litvmu(p,w(1:p),l,w)
t1=one/v2norm(p,w)
t=t1*t
if(t<=delta)goto 370
if(k>30)goto 270
k=k+1

! Start next inv. power iter. by storing normalized w in x.
w(x0+1:x0+p)=t1*w(1:p)
! Compute w = (l**-1)*x.
call livmul(p,w(1:p),l,w(x:))
t=one/v2norm(p,w)
goto 340

370 w(1:p)=t1*w(1:p)

! Now w is the desired approximate (unit) eigenvector and
! t*x = ((d**-1)*h*(d**-1) + alphak*i)*w.

sw=dotprd(p,w(q:),w)
t1=(rad+dst)*(rad-dst)
root=sqrt(sw*sw+t1)
if(sw<zero)root=-root
si=t1/(sw+root)
! Accept current step if adding si*w would lead to a
! further relative reduction in psi of less than v(epslon)/3.
v(preduc)=half*twopsi
t1=zero
t=si*(alphak*sw-half*si*(alphak+t*dotprd(p,w(x:),w)))
if(t<epso6*twopsi)goto 390
v(preduc)=v(preduc)+t
dst=rad
t1=-si
390 do i=1,p
 j=q0+i
 w(j)=t1*w(i)-w(j)
 step(i)=w(j)/d(i)
enddo
v(gtstep)=dotprd(p,dig,w(q:))

! Save values for use in a possible restart

410 v(dstnrm)=dst
v(stppar)=alphak
w(lk0)=lk
w(uk0)=uk
v(rad0)=rad
w(dstsav)=dst

! Restore diagonal of dihdi

j=0
do i=1,p
 j=j+i
 k=diag0+i
 dihdi(j)=w(k)
enddo
goto 450

! Special case -- g = 0

430 v(stppar)=zero
v(preduc)=zero
v(dstnrm)=zero
v(gtstep)=zero
step(1:p)=zero

450 return

END SUBROUTINE gqtstp


SUBROUTINE itsmry (d,iv,p,v,x)
!-----------------------------------------------!
! Print NL2SOL (version 2.2) iteration summary  !
!-----------------------------------------------!
INTEGER,INTENT(in) :: p
INTEGER,INTENT(inout) :: iv(:)
REAL(dp),INTENT(in) :: d(:),v(:),x(:)
INTEGER cov1,g1,i,ii,iv1,i1,m,nf,ng,ol,pu
!CHARACTER(4) :: model1(6) = (/ '    ', '    ', '    ', '    ', '  g ',  &
!                  &'  s ' /), &
!                &model2(6) = (/ ' g  ', ' s  ', 'g-s ', 's-g ', '-s-g',  &
!                  &'-g-s' /)
REAL(dp) nreldf,oldf,preldf,reldf

!  ***  no external functions or subroutines  ***

!  ***  iv subscript values  ***

INTEGER,PARAMETER :: covmat=26,covprt=14,g=28,covreq=15,needhd=39,nfcall=6,&
 &nfcov=40,ngcov=41,ngcall=30,niter=31,outlev=19,prntit=48,prunit=21,&
 &solprt=22,statpr=23,sused=57,x0prt=24

!  ***  v subscript values  ***

INTEGER,PARAMETER :: dstnrm=2,f=10,f0=13,fdif=11,nreduc=6,preduc=7,reldx=17,&
 &size=47,stppar=5

REAL(dp),PARAMETER :: zero=0.d0

pu=iv(prunit)
if(pu<=0)goto 610
iv1=iv(1)
ol=iv(outlev)
if(iv1<2.or.iv1>15)goto 320
if(ol==0)goto 70
if(iv1>=12)goto 70
if(iv1>=10.and.iv(prntit)==0)goto 70
if(iv1>2)goto 10
iv(prntit)=iv(prntit)+1
if(iv(prntit)<abs(ol))goto 610
10 nf=iv(nfcall)-abs(iv(nfcov))
iv(prntit)=0
reldf=zero
preldf=zero
oldf=v(f0)
if(oldf<=zero)goto 20
reldf=v(fdif)/oldf
preldf=v(preduc)/oldf
20 if(ol>0)goto 40

! Print short summary line
iv(needhd)=0
write(pu,60)iv(niter),nf,v(f),reldf,preldf,v(reldx)
goto 70

! Print long summary line
40 continue
iv(needhd)=0
m=iv(sused)
nreldf=zero
if(oldf>zero)nreldf=v(nreduc)/oldf
write(pu,*)'NL2SOL iteration summary :'
write(pu,*)'--------------------------'
write(pu,'(a,1x,i16)')' Iteration                  :',iv(niter)
write(pu,'(a,1x,i16)')' No. function evals.        :',nf
write(pu,'(a,1x,e16.8)')' Function value             :',v(f)
write(pu,'(a,1x,e16.8)')' Relative function diff.    :',reldf
write(pu,'(a,1x,e16.8)')' Predicted rel. func. diff. :',preldf
write(pu,'(a,1x,e16.8)')' Rel. parameter change      :',v(reldx)
select case(m)
case(1)
 write(pu,'(a)')' Quadratic model            :   Gauss-Newton'
case(2)
 write(pu,'(a)')' Quadratic model            :   Augmented'
case(3)
 write(pu,'(a)')' Quadratic models           :   Gauss-Newton + Augmented'
case(4)
 write(pu,'(a)')' Quadratic models           :   Augmented + Gauss-Newton'
case(5)
 write(pu,'(a)')' Quadratic models           :   Gauss-Newton + Augmented +&
  & Gauss-Newton'
case(6)
 write(pu,'(a)')' Quadratic models           :   Augmented + Gauss-Newton +&
  & Augmented'
endselect
write(pu,'(a,1x,e16.8)')' Marquardt parameter        :',v(stppar)
write(pu,'(a,1x,e16.8)')' Sizing factor              :',v(size)
write(pu,'(a,1x,e16.8)')' Norm of d times step       :',v(dstnrm)
write(pu,'(a,1x,e16.8)')' Convergence parameter      :',nreldf
write(pu,*)

60 format(' ', i5, i6, 4F11.3, a3, a4, 4F11.3)

70 select case (iv1)
case (1:2) ; goto 610
case (3)
 write(pu,90) ; goto 370
 90 format(' NL2SOL: parameter convergence.'/)
case (4)
 write(pu,110) ; goto 370
 110 format(' NL2SOL: relative function convergence.'/)
case (5)
 write(pu,130) ; goto 370
 130 format(' NL2SOL: parameter and relative function convergence.'/)
case (6)
 write(pu,150) ; goto 370
 150 format(' NL2SOL: absolute function convergence.'/)
case (7)
 write(pu,170) ; goto 370
 170 format(' NL2SOL: singular convergence.'/)
case (8)
 write(pu,190) ; goto 370
 190 format(' NL2SOL: false convergence.'/)
case (9)
 write(pu,210) ; goto 370
 210 format(' NL2SOL: function evaluation limit.'/)
case (10)
 write(pu,230) ; goto 370
 230 format(' NL2SOL: iteration limit.'/)
case (11)
 write(pu,250) ; goto 370
 250 format(' NL2SOL: halting due to stopx.'/)
case (12) ; goto 340
case (13)
 write(pu,270) ; goto 340
 270 format(' NL2SOL: initial sum of squares overflows.'/)
case (14)
 write(pu,290) ; goto 610
 290 format(' NL2SOL: bad parameters to assess.'/)
case (15)
 write(pu,310)
 310 format(' NL2SOL: Jacobian could not be computed.'/)
 if(iv(niter)>0)goto 420
 goto 340
endselect

320 write(pu,330)iv1
330 format(' NL2SOL: iv(1) is ',i5/)
goto 610

! Initial call on itsmry
340 if(iv(x0prt)/=0)then
 write(pu,*)'NL2SOL initial parameters :'
 write(pu,*)'---------------------------'
 write(pu,*)'    i              x(i)           d(i)'
 write(pu,*)'----- ----------------- --------------'
 write(pu,350)(i,x(i),d(i),i=1,p)
 write(pu,*)'--------------------------------------'
endif
350 format(1x,i5,1x,f17.6,1x,f14.3)
if(iv1>=13)goto 610
iv(needhd)=0
iv(prntit)=0
if(ol==0)goto 610
write(pu,360)v(f)
360 format(' Initial function value: ',e16.8/)
goto 610

! Print various information requested on solution.
370 iv(needhd)=1
if(iv(statpr)==0)goto 420
oldf=v(f0)
preldf=zero
nreldf=zero
if(oldf<=zero)goto 380
preldf=v(preduc)/oldf
nreldf=v(nreduc)/oldf
380 nf=iv(nfcall)-iv(nfcov)
ng=iv(ngcall)-iv(ngcov)
write(pu,*)'NL2SOL final summary :'
write(pu,*)'----------------------'
write(pu,'(a,1x,e16.8)')' Final function value       :',v(f)
write(pu,'(a,1x,e16.8)')' Rel. parameter change      :',v(reldx)
write(pu,'(a,1x,i16)')' Number of function evals.  :',nf
write(pu,'(a,1x,i16)')' Number of gradient evals.  :',ng
write(pu,'(a,1x,e16.8)')' Predicted rel. func. diff. :',preldf
write(pu,'(a,1x,e16.8)')' Convergence parameter      :',nreldf

if(iv(nfcov)>0)write(pu,400)iv(nfcov)
400 format (' Func. evals for covariance : ',i16)
if(iv(ngcov)>0)write(pu,410)iv(ngcov)
410 format (' Grad. evals for covariance : ',i16)
write(pu,*)

420 if(iv(solprt)==0)goto 460
iv(needhd)=1
g1=iv(g)
write(pu,*)'Final parameters :'
write(pu,*)'------------------'
write(pu,*)'    i              x(i)           d(i)           g(i)'
write(pu,*)'----- ----------------- -------------- --------------'
do i=1,p
 write(pu,450)i,x(i),d(i),v(g1)
 g1=g1+1
enddo
write(pu,*)'-----------------------------------------------------'
write(pu,*)
450 format(1x,i5,1x,f17.6,2(1x,f14.3))

460 if(iv(covprt)==0)goto 610
cov1=iv(covmat)
iv(needhd)=1
if(cov1<0.0)then
 goto 470
elseif(cov1==0.0)then
 goto 500
else
 goto 520
endif
470 if(-1==cov1)write(pu,480)
480 format(' NL2SOL: indefinite covariance matrix.'/)
if(-2==cov1)write(pu,490)
490 format(' NL2SOL: oversize steps in computing covariance.'/)
goto 610

500 write(pu,510)
510 format(' NL2SOL: covariance matrix not computed'/)
goto 610

520 i=abs(iv(covreq))
write(pu,*)'NL2SOL print-out of covariance matrix :'
write(pu,*)'---------------------------------------'
if(i<=1)write(pu,530)
530 format(' Cov. matrix = scale * H**-1 * (J**t * J) * H**-1'/)
if(i==2)write(pu,540)
540 format (' Cov. matrix = scale * H**-1'/)
if(i>=3)write(pu,550)
550 format(' Cov. matrix = scale * (J**t * J)**-1'/)
ii=cov1-1
if(ol<=0)goto 580
do i=1,p
 i1=ii+1
 ii=ii+i
 write(pu,'(a,i3,a)')'Row ',i,' is:'
 write(pu,570)v(i1:ii)
enddo
write(pu,*)
570 format(4(1x,e16.8))
goto 610

580 do i=1,p
 i1=ii+1
 ii=ii+i
 write(pu,'(a,i3,a)')'Row ',i,' is:'
 write(pu,600)v(i1:ii)
enddo
write(pu,*)
600 format(4(1x,e16.8))

610 return

END SUBROUTINE itsmry


SUBROUTINE linvrt(n,lin,l)
!------------------------------------------------------------!
! Compute  lin = l**-1,  both  n x n  lower triang. stored   !
! compactly by rows.  lin and l may share the same storage.  !
!------------------------------------------------------------!
INTEGER,INTENT(in) :: n
REAL(dp),INTENT(in) :: l(:)
REAL(dp),INTENT(inout) :: lin(:)
INTEGER i,ii,im1,jj,j0,j1,k,k0,np1
REAL(dp) t
REAL(dp),PARAMETER :: one=1.d0,zero=0.d0

np1=n+1
j0=n*(np1)/2
do ii=1,n
 i=np1-ii
 lin(j0)=one/l(j0)
 if(i <= 1) goto 40
 j1=j0
 im1=i-1
 do jj=1,im1
  t=zero
  j0=j1
  k0=j1-jj
  do k=1,jj
   t=t-l(k0)*lin(j0)
   j0=j0-1
   k0=k0+k-i
  enddo
  lin(j0)=t/l(k0)
 enddo
 j0=j0-1
enddo
40 return

END SUBROUTINE linvrt


SUBROUTINE litvmu (n, x, l, y)
!---------------------------------------------------------------------------!
! Solve  (l**t)*x = y,  where  l  is an  n x n  lower triangular            !
! matrix stored compactly by rows.  x and y may occupy the same storage.    !
!---------------------------------------------------------------------------!
INTEGER,INTENT(in) :: n
REAL(dp),INTENT(in) :: l(:),y(:)
REAL(dp),INTENT(inout) :: x(:)
INTEGER i,ii,ij,im1,i0,j,np1
REAL(dp) xi
REAL(dp),PARAMETER :: zero=0.d0

x(1:n)=y(1:n)
np1=n+1
i0=n*(n+1)/2
do ii=1,n
 i=np1-ii
 xi=x(i)/l(i0)
 x(i)=xi
 if(i<=1)goto 40
 i0=i0-i
 if(xi==zero)cycle
 im1=i-1
 do j=1,im1
  ij=i0+j
  x(j)=x(j)-xi*l(ij)
 enddo
enddo
40 return

END SUBROUTINE litvmu


SUBROUTINE livmul(n,x,l,y)
!-------------------------------------------------------------------------!
! Solve  l*x = y, where  l  is an  n x n  lower triangular                !
! matrix stored compactly by rows.  x and y may occupy the same storage.  !
!-------------------------------------------------------------------------!
INTEGER,INTENT(in) :: n
REAL(dp),INTENT(in) :: l(:),y(:)
REAL(dp),INTENT(inout) :: x(:)
INTEGER i,j,k
REAL(dp) t
REAL(dp),PARAMETER :: zero=0.d0

do k=1,n
 if(y(k)/=zero)goto 20
 x(k)=zero
enddo
goto 40
20 j=k*(k+1)/2
x(k)=y(k)/l(j)
if(k>=n)goto 40
k=k+1
do i=k,n
 t=dotprd(i-1,l(j+1:),x)
 j=j+i
 x(i)=(y(i)-t)/l(j)
enddo
40 return
END SUBROUTINE livmul


SUBROUTINE lmstep(d,g,ierr,ipivot,ka,p,qtr,r,step,v,w)
!-----------------------------------------------------------------------------!
! Compute Levenberg-Marquardt step using More-Hebden technique.
! NL2SOL version 2.2.
!
!     Given the r matrix from the QR decomposition of a Jacobian
!     matrix, j, as well as q-transpose times the corresponding
!     residual vector, resid, this subroutine computes a Levenberg-
!     Marquardt step of approximate length v(radius) by the more-
!     technique.
!
!      d (in)  = the scale vector.
!      g (in)  = the gradient vector (j**t)*r.
!   ierr (i/o) = return code from qrfact or qrfgs -- 0 means r has
!             full rank.
! ipivot (i/o) = permutation array from qrfact or qrfgs, which compute
!             qr decompositions with column pivoting.
!     ka (i/o).  ka<0 on input means this is the first call on
!             lmstep for the current r and qtr.  on output ka con-
!             tains the number of hebden iterations needed to determine
!             step.  ka = 0 means a gauss-newton step.
!      p (in)  = number of parameters.
!    qtr (in)  = (q**t)*resid = q-transpose times the residual vector.
!      r (in)  = the r matrix, stored compactly by columns.
!   step (out) = the levenberg-marquardt step computed.
!      v (i/o) contains various constants and variables described below.
!      w (i/o) = workspace of length p*(p+5)/2 + 4.
!
!  ***  entries in v  ***
!
! v(dgnorm) (i/o) = 2-norm of (d**-1)*g.
! v(dstnrm) (i/o) = 2-norm of d*step.
! v(dst0)   (i/o) = 2-norm of gauss-newton step (for nonsing. j).
! v(epslon) (in) = max. rel. error allowed in twonorm(r)**2 minus
!             twonorm(r - j*step)**2.  (see algorithm notes below.)
! v(gtstep) (out) = inner product between g and step.
! v(nreduc) (out) = half the reduction in the sum of squares predicted
!             for a gauss-newton step.
! v(phmnfc) (in)  = tol. (together with v(phmxfc)) for accepting step
!             (more*s sigma).  the error v(dstnrm) - v(radius) must lie
!             between v(phmnfc)*v(radius) and v(phmxfc)*v(radius).
! v(phmxfc) (in)  (see v(phmnfc).)
! v(preduc) (out) = half the reduction in the sum of squares predicted
!             by the step returned.
! v(radius) (in)  = radius of current (scaled) trust region.
! v(rad0)   (i/o) = value of v(radius) from previous call.
! v(stppar) (i/o) = marquardt parameter (or its negative if the special
!             case mentioned below in the algorithm notes occurs).
!
! note -- see data statement below for values of above subscripts.
!
!  ***  usage notes  ***
!
!     if it is desired to recompute step using a different value of
!     v(radius), then this routine may be restarted by calling it
!     with all parameters unchanged except v(radius).  (this explains
!     why many parameters are listed as i/o).  on an intiial call (one
!     with ka = -1), the caller need only have initialized d, g, ka, p,
!     qtr, r, v(epslon), v(phmnfc), v(phmxfc), v(radius), and v(rad0).
!
!  ***  application and usage restrictions  ***
!
!     this routine is called as part of the nl2sol (nonlinear least-
!     squares) package (ref. 1).
!
!  ***  algorithm notes  ***
!
!     this code implements the step computation scheme described in
!     refs. 2 and 4.  fast givens transformations (see ref. 3, pp. 60-
!     62) are used to compute step with a nonzero marquardt parameter.
!        a special case occurs if j is (nearly) singular and v(radius)
!     is sufficiently large.  in this case the step returned is such
!     that  twonorm(r)**2 - twonorm(r - j*step)**2  differs from its
!     optimal value by less than v(epslon) times this optimal value,
!     where j and r denote the original jacobian and residual.  (see
!     ref. 2 for more details.)
!
!  ***  functions and subroutines called  ***
!
! dotprd - returns inner product of two vectors.
! litvmu - apply inverse-transpose of compact lower triang. matrix.
! livmul - apply inverse of compact lower triang. matrix.
! vcopy  - copies one vector to another.
! v2norm - returns 2-norm of a vector.
!
!  ***  references  ***
!
! 1.  dennis, j.e., gay, d.m., and welsch, r.e. (1981), an adaptive
!             nonlinear least-squares algorithm, acm trans. math.
!             software, vol. 7, no. 3.
! 2.  gay, d.m. (1981), computing optimal locally constrained steps,
!             siam j. sci. statist. computing, vol. 2, no. 2, pp.
!             186-197.
! 3.  lawson, c.l., and hanson, r.j. (1974), solving least squares
!             problems, prentice-hall, englewood cliffs, n.j.
! 4.  more, j.j. (1978), the levenberg-marquardt algorithm, implemen-
!             tation and theory, pp.105-116 of springer lecture notes
!             in mathematics no. 630, edited by g.a. watson, springer-
!             verlag, berlin and new york.
!
!  ***  general  ***
!
!     coded by david m. gay.
!     this subroutine was written in connection with research
!     supported by the national science foundation under grants
!     mcs-7600324, dcr75-10143, 76-14311dss, mcs76-11989, and
!     mcs-7906671.
!
!------------------------------------------------------------------------------
INTEGER,INTENT(in) :: p
INTEGER,INTENT(inout) :: ierr,ipivot(:),ka
REAL(dp),INTENT(in) :: d(:),g(:),qtr(:),r(:)
REAL(dp),INTENT(inout) :: v(:),w(:)
REAL(dp),INTENT(inout) :: step(:)
INTEGER dstsav,i,ip1,i1,j1,k,kalim,l,lk0,phipin,pp1o2,res,res0,rmat,rmat0,uk0
INTEGER,PARAMETER :: dgnorm=1,dstnrm=2,dst0=3,epslon=19,gtstep=4,&
 &nreduc=6,phmnfc=20,phmxfc=21,preduc=7,radius=8,rad0=9,stppar=5
REAL(dp) a,adi,alphak,b,dfacsq,dst,dtol,d1,d2,lk,oldphi,phi,phimax,phimin, &
 &psifac,rad,si,sj,sqrtak,t,twopsi,uk,wl
REAL(dp),PARAMETER :: dfac=256.d0,eight=8.d0,half=0.5d0, &
 &negone=-1.d0,one=1.d0,p001=1.d-3,three=3.d0,ttol=2.5d0,zero=0.d0

psifac=0.d0 ; alphak=0.d0

! For use in recomputing step, the final values of lk and uk,
! the inverse derivative of more*s phi at 0 (for nonsing. j)
! and the value returned as v(dstnrm) are stored at w(lk0),
! w(uk0), w(phipin), and w(dstsav) respectively.
lk0=p+1
phipin=lk0+1
uk0=phipin+1
dstsav=uk0+1
rmat0=dstsav

! A copy of the r-matrix from the qr decomposition of j is
! stored in w starting at w(rmat), and a copy of the residual
! vector is stored in w starting at w(res).  The loops below
! that update the qr decomp. For a nonzero marquardt parameter
! work on these copies.
rmat=rmat0+1
pp1o2=p*(p+1)/2
res0=pp1o2+rmat0
res=res0+1
rad=v(radius)
if(rad>zero)psifac=v(epslon)/((eight*(v(phmnfc)+one)+three)* rad**2)
phimax=v(phmxfc)*rad
phimin=v(phmnfc)*rad

! dtol, dfac, and dfacsq are used in rescaling the fast Givens representation
! of the updated QR decomposition.
dtol=one/dfac
dfacsq=dfac*dfac

! oldphi is used to detect limits of numerical accuracy.  If
! we recompute step and it does not change, then we accept it.
oldphi=zero
lk=zero
uk=zero
kalim=ka+12

! Start or restart, depending on ka

if(ka<0)then
 goto 10
elseif(ka==0)then
 goto 20
else
 goto 310
endif

! Fresh start -- compute v(nreduc)

10 ka=0
kalim=12
k=p
if(ierr/=0)k=abs(ierr)-1
v(nreduc)=half*dotprd(k,qtr,qtr)

! Set up to try initial Gauss-Newton step

20 v(dst0)=negone
if(ierr/=0)goto 50

! Compute Gauss-Newton step

! Note -- the r-matrix is stored compactly by columns in
! r(1), r(2), r(3), ...  it is the transpose of a
! lower triangular matrix stored compactly by rows, and we
! treat it as such when using litvmu and livmul.

call litvmu(p,w(1:p),r,qtr)

! Temporarily store permuted -d*step in step.

do i=1,p
 j1=ipivot(i)
 step(i)=d(j1)*w(i)
enddo
dst=v2norm(p,step)
v(dst0)=dst
phi=dst-rad
if(phi<=phimax)goto 350

! If this is a restart, go to 110
if(ka>0)goto 70

! Gauss-Newton step was unacceptable.  compute l0.

do i=1,p
 j1=ipivot(i)
 step(i)=d(j1)*(step(i)/dst)
enddo
call livmul(p,step(1:p),r,step)
t=one/v2norm(p,step)
w(phipin)=(t/dst)*t
lk=phi*w(phipin)

! Compute u0.

50 w(1:p)=g(1:p)/d(1:p)
v(dgnorm)=v2norm(p,w)
uk=v(dgnorm)/rad
if(uk<=zero)goto 330

! alphak will be used as the current Marquardt parameter. We
! use more*s scheme for initializing it.

alphak=abs(v(stppar))*v(rad0)/rad

! Top of loop -- increment ka, copy r to rmat, qtr to res

70 ka=ka+1
call vcopy(pp1o2,w(rmat:rmat+pp1o2-1),r)
call vcopy(p,w(res:res+p-1),qtr)

! Safeguard alphak and initialize fast givens scale vector.

if(alphak<=zero.or.alphak<lk.or.alphak>=uk)alphak=uk*max(p001,sqrt(lk/uk))
sqrtak=sqrt(alphak)
w(1:p)=one

! Add alphak*d and update qr decomp. using fast givens trans.

do i=1,p
! Generate, apply 1st givens trans. for row i of alphak*d.
! (use step to store temporary row)  ***
 l=i*(i+1)/2+rmat0
 wl=w(l)
 d2=one
 d1=w(i)
 j1=ipivot(i)
 adi=sqrtak*d(j1)
 if(adi >= ABS(wl)) goto 110
 90 a=adi/wl
 b=d2*a/d1
 t=a*b+one
 if(t > ttol) goto 110
 w(i)=d1/t
 d2=d2/t
 w(l)=t*wl
 a=-a
 do j1=i,p
  l=l+j1
  step(j1)=a*w(l)
 enddo
 goto 130

 110 b=wl/adi
 a=d1*b/d2
 t=a*b+one
 if(t>ttol)goto 90
 w(i)=d2/t
 d2=d1/t
 w(l)=t*adi
 do j1=i,p
  l=l+j1
  wl=w(l)
  step(j1)=-wl
  w(l)=a*wl
 enddo

 130 if(i==p)goto 240

! Now use Givens trans. to zero elements of temp. row

 ip1=i+1
 do i1=ip1,p
  l=i1*(i1+1)/2+rmat0
  wl=w(l)
  si=step(i1-1)
  d1=w(i1)

! Rescale row i1 if necessary.

  if(d1>=dtol)goto 150
  d1=d1*dfacsq
  wl=wl/dfac
  k=l
  do j1=i1,p
   k=k+j1
   w(k)=w(k)/dfac
  enddo

! Use Givens trans. to zero next element of temp. row

  150 if(abs(si)>abs(wl))goto 180
  if(si==zero)cycle
  160 a=si/wl
  b=d2*a/d1
  t=a*b+one
  if(t>ttol)goto 180
  w(l)=t*wl
  w(i1)=d1/t
  d2=d2/t
  do j1=i1,p
   l=l+j1
   wl=w(l)
   sj=step(j1)
   w(l)=wl+b*sj
   step(j1)=sj-a*wl
  enddo
  goto 200

  180 b=wl/si
  a=d1*b/d2
  t=a*b+one
  if(t>ttol)goto 160
  w(i1)=d2/t
  d2=d1/t
  w(l)=t*si
  do j1=i1,p
   l=l+j1
   wl=w(l)
   sj=step(j1)
   w(l)=a*wl+sj
   step(j1)=b*sj-wl
  enddo

! Rescale temp. row if necessary.

  200 if(d2>=dtol)cycle
  d2=d2*dfacsq
  step(i1:p)=step(i1:p)/dfac
 enddo
enddo

! Compute step.

240 call litvmu(p,w(res:res+p-1),w(rmat:),w(res:))
! Recover step and store permuted -d*step at w(res).
do i=1,p
 j1=ipivot(i)
 k=res0+i
 t=w(k)
 step(j1)=-t
 w(k)=t*d(j1)
enddo
dst=v2norm(p,w(res:))
phi=dst-rad
if(phi<=phimax.and.phi>=phimin)goto 370
if(oldphi==phi)goto 370
oldphi=phi

! Check for (and handle) special case.

if(phi>zero)goto 270
if(ka>=kalim)goto 370
twopsi=alphak*dst*dst-dotprd(p,step,g)
if(alphak>=twopsi*psifac)goto 270
v(stppar)=-alphak
goto 380

! Unacceptable step -- update lk, uk, alphak, and try again.

260 if(phi<zero)uk=min(uk,alphak)
goto 280
270 if(phi<zero)uk=alphak
280 do i=1,p
 j1=ipivot(i)
 k=res0+i
 step(i)=d(j1)*(w(k)/dst)
enddo
call livmul(p,step(1:p),w(rmat:),step)
step(1:p)=step(1:p)/sqrt(w(1:p))
t=one/v2norm(p,step)
alphak=alphak+t*phi*t/rad
lk=max(lk,alphak)
goto 70

! Restart.

310 lk=w(lk0)
uk=w(uk0)
if(v(dst0)>zero.and.v(dst0)-rad<=phimax)goto 20
alphak=abs(v(stppar))
dst=w(dstsav)
phi=dst-rad
t=v(dgnorm)/rad
if(rad>v(rad0))goto 320

! Smaller radius
uk=t
if(alphak<=zero)lk=zero
if(v(dst0)>zero)lk=max(lk,(v(dst0)-rad)*w(phipin))
goto 260

! Bigger radius
320 if(alphak<=zero.or.uk>t)uk=t
lk=zero
if(v(dst0)>zero)lk=max(lk,(v(dst0)-rad)*w(phipin))
goto 260

! Special case -- rad <= 0 or (g = 0 and j is singular).

330 v(stppar)=zero
dst=zero
lk=zero
uk=zero
v(gtstep)=zero
v(preduc)=zero
step(1:p)=zero
goto 390

! Acceptable Gauss-Newton step -- recover step from w.

350 alphak=zero
do i=1,p
 j1=ipivot(i)
 step(j1)=-w(i)
enddo

! Save values for use in a possible restart.

370 v(stppar)=alphak
380 v(gtstep)=dotprd(p,step,g)
v(preduc)=half*(alphak*dst*dst-v(gtstep))
390 v(dstnrm)=dst
w(dstsav)=dst
w(lk0)=lk
w(uk0)=uk
v(rad0)=rad

END SUBROUTINE lmstep


SUBROUTINE lsqrt (n1, n, l, a, irc)
!---------------------------------------------------------------------!
! Compute rows n1 through n of the Cholesky factor  l  of             !
! a = l*(l**t),  where  l  and the lower triangle of  a  are both     !
! stored compactly by rows (and may occupy the same storage).         !
! irc = 0 means all went well.  irc = j means the leading             !
! principal  j x j  submatrix of  a  is not positive definite --      !
! and  l(j*(j+1)/2)  contains the (nonpos.) reduced j-th diagonal.    !
!---------------------------------------------------------------------!
INTEGER,INTENT(in) :: n1, n
INTEGER,INTENT(out) :: irc
REAL(dp),INTENT(in) :: a(:)
REAL(dp),INTENT(inout) :: l(:)
INTEGER i,ij,ik,im1,i0,j,jk,jm1,j0,k
REAL(dp) :: t,td,zero=0.d0

i0=n1*(n1-1)/2
do i=n1,n
 td=zero
 if(i==1)goto 40
 j0=0
 im1=i-1
 do j=1,im1
  t=zero
  if(j==1)goto 20
  jm1=j-1
  do k=1,jm1
   ik=i0+k
   jk=j0+k
   t=t+l(ik)*l(jk)
  enddo
  20 ij=i0+j
  j0=j0+j
  t=(a(ij)-t)/l(j0)
  l(ij)=t
  td=td+t*t
 enddo
 40 i0=i0+i
 t=a(i0)-td
 if(t<=zero)goto 60
 l(i0)=sqrt(t)
enddo

irc=0
goto 70

60 l(i0)=t
irc=i

70 return

END SUBROUTINE lsqrt


SUBROUTINE lsvmin (p,l,x,y,fn_val)
!------------------------------------------------------------------!
! Estimate smallest sing. value of packed lower triang. matrix l   !
!------------------------------------------------------------------!
INTEGER,INTENT(in) :: p
REAL(dp),INTENT(in) :: l(:)
REAL(dp),INTENT(inout) :: x(:), y(:)
REAL(dp),INTENT(out) :: fn_val
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     This function returns a good over-estimate of the smallest
!     singular value of the packed lower triangular matrix l.
!
!  p (in)  = the order of l.  l is a  p x p  lower triangular matrix.
!  l (in)  = array holding the elements of  l  in row order, i.e.
!             l(1,1), l(2,1), l(2,2), l(3,1), l(3,2), l(3,3), etc.
!  x (out) if lsvmin returns a positive value, then x is a normalized
!             approximate left singular vector corresponding to the
!             smallest singular value.  this approximation may be very
!             crude.  if lsvmin returns zero, then some components of x
!             are zero and the rest retain their input values.
!  y (out) if lsvmin returns a positive value, then y = (l**-1)*x is an
!             unnormalized approximate right singular vector correspond-
!             ing to the smallest singular value.  this approximation
!             may be crude.  if lsvmin returns zero, then y retains its
!             input value.  the caller may pass the same vector for x
!             and y (nonstandard fortran usage), in which case y over-
!             writes x (for nonzero lsvmin returns).

!  ***  application and usage restrictions  ***

!     there are no usage restrictions.

!  ***  algorithm notes  ***

!     the algorithm is based on (1), with the additional provision that
!     lsvmin = 0 is returned if the smallest diagonal element of l
!     (in magnitude) is not more than the unit roundoff times the
!     largest.  the algorithm uses a random number generator proposed
!     in (4), which passes the spectral test with flying colours -- see
!     (2) and (3).

!  ***  subroutines and functions called  ***

!        v2norm - function, returns the 2-norm of a vector.

!  ***  references  ***

!     (1) cline, a., moler, c., stewart, g., and wilkinson, j.h.(1977),
!         an estimate for the condition number of a matrix, report
!         tm-310, applied math. div., argonne national laboratory.

!     (2) hoaglin, d.c. (1976), theoretical properties of congruential
!         random-number generators --  an empirical view,
!         memorandum ns-340, dept. of statistics, harvard univ.

!     (3) knuth, d.e. (1969), the art of computer programming, vol. 2
!         (seminumerical algorithms), addison-wesley, reading, mass.

!     (4) smith, c.s. (1971), multiplicative pseudo-random number
!         generators with prime modulus, j. assoc. comput. mach. 18,
!         pp. 586-593.

!  ***  history  ***

!     designed and coded by david m. gay (winter 1977/summer 1978).

!  ***  general  ***

!     this subroutine was written in connection with research
!     supported by the national science foundation under grants
!     mcs-7600324, dcr75-10143, 76-14311dss, and mcs76-11989.

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

INTEGER i,ii,j,ji,jj,jjj,jm1,j0,pplus1
REAL(dp) b,psj,sminus,splus,t,xminus,xplus

REAL(dp),PARAMETER :: half=0.5d0,one=1.d0,r9973=9973.d0,zero=0.d0
INTEGER,SAVE :: ix=2

! First check whether to return lsvmin = 0 and initialize x.

ii=0
do i=1,p
 x(i)=zero
 ii=ii+i
 if(l(ii)==zero)goto 100
enddo
if(mod(ix,9973)==0)ix=2
pplus1=p+1

! Solve (l**t)*x = b, where the components of b have randomly
! chosen magnitudes in (.5,1) with signs chosen to make x large.

! do j = p to 1 by -1...
do jjj=1,p
 j=pplus1-jjj
! Determine x(j) in this iteration. note for i = 1,2,...,j
! that x(i) holds the current partial sum for row i.
 ix=mod(3432*ix,9973)
 b=half*(one+real(ix,dp)/r9973)
 xplus=(b-x(j))
 xminus=(-b-x(j))
 splus=abs(xplus)
 sminus=abs(xminus)
 jm1=j-1
 j0=j*jm1/2
 jj=j0+j
 xplus=xplus/l(jj)
 xminus=xminus/l(jj)
 if(jm1 == 0) goto 30
 do i=1,jm1
  ji=j0+i
  splus=splus+abs(x(i)+l(ji)*xplus)
  sminus=sminus+abs(x(i)+l(ji)*xminus)
 enddo
 30 if(sminus>splus)xplus=xminus
 x(j)=xplus
! Update partial sums.
 if(jm1==0)cycle
 do i=1,jm1
  ji=j0+i
  x(i)=x(i)+l(ji)*xplus
 enddo
enddo

! Normalize x.

t=one/v2norm(p,x)
x(1:p)=t*x(1:p)

! Solve l*y = x and return svmin = 1/twonorm(y)

do j=1,p
 psj=zero
 jm1=j-1
 j0=j*jm1/2
 if(jm1==0)goto 80
 do i=1,jm1
  ji=j0+i
  psj=psj+l(ji)*y(i)
 enddo
 80 jj=j0+j
 y(j)=(x(j)-psj)/l(jj)
enddo

fn_val=one/v2norm(p,y)
goto 110

100 fn_val=zero
110 return
END SUBROUTINE lsvmin


SUBROUTINE ltsqar(n,a,l)
!----------------------------------------------------------------!
! Set a to lower triangle of (l**t) * l                          !
!                                                                !
! l = n x n lower triang. matrix stored rowwise.                 !
! a is also stored rowwise and may share storage with l.         !
!----------------------------------------------------------------!
INTEGER,INTENT(in) :: n
REAL(dp),INTENT(in) :: l(:)
REAL(dp),INTENT(inout) :: a(:)
INTEGER i,ii,iim1,i1,j,k,m
REAL(dp) lii,lj

ii=0
do i=1,n
 i1=ii+1
 ii=ii+i
 m=1
 if(i==1)goto 30
 iim1=ii-1
 do j=i1,iim1
  lj=l(j)
  do k=i1,j
   a(m)=a(m)+lj*l(k)
   m=m+1
  enddo
 enddo
 30 lii=l(ii)
 do j=i1,ii
  a(j)=lii*l(j)
 enddo
enddo

END SUBROUTINE ltsqar


SUBROUTINE parchk(iv,n,nn,p,v)
!-----------------------------------------------------------------!
! Check NL2SOL (version 2.2) parameters, print changed values.    !
!-----------------------------------------------------------------!
INTEGER,INTENT(in) :: n,nn,p
INTEGER,INTENT(inout) :: iv(:)
REAL(dp),INTENT(inout) :: v(:)

! dfault -- supplies dfault parameter values.
! rmdcon -- returns machine-dependent constants.
! vcopy  -- copies one vector to another.

INTEGER i,iv1,jtolp,k,l,m,pu
! CHARACTER(4) :: which(3)
REAL(dp) machep,vk

! iv and v subscripts

INTEGER,PARAMETER :: nvdflt=27
REAL(dp),PARAMETER :: zero=0.d0

INTEGER,PARAMETER :: dtype=16,dtype0=29,d0init=37,epslon=19,inits=25, &
 &jtinit=39,jtol0=86,jtol1=87,oldn=45,oldnn=46,oldp=47,parprt=20,parsv1=51,&
 &prunit=21

REAL(dp),SAVE :: big=0.d0,tiny=1.d0
!!$CHARACTER(4),PARAMETER :: vn(2,27) = RESHAPE( (/  &
!!$                                          &'epsl','on..', 'phmn','fc..',  &
!!$                                          &'phmx','fc..', 'decf','ac..',  &
!!$                                          &'incf','ac..', 'rdfc','mn..',  &
!!$                                          &'rdfc','mx..', 'tune','r1..',  &
!!$                                          &'tune','r2..', 'tune','r3..',  &
!!$                                          &'tune','r4..', 'tune','r5..',  &
!!$                                          &'afct','ol..', 'rfct','ol..',  &
!!$                                          &'xcto','l...', 'xfto','l...',  &
!!$                                          &'lmax','0...', 'dltf','dj..',  &
!!$                                          &'d0in','it..', 'dini','t...',  &
!!$                                          &'jtin','it..', 'dltf','dc..',  &
!!$                                          &'dfac','....', 'rlim','it..',  &
!!$                                          &'cosm','in..', 'delt','a0..',  &
!!$                                          &'fuzz','....' /), (/ 2, 27 /) )

REAL(dp) :: vm(27) = (/ 1.0D-3, -0.99d0, 1.0D-3, 1.0D-2, &
              &1.2d0, 1.d-2, 1.2d0, 0.d0, 0.d0, 1.d-3, -1.d0, &
              &0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,    &
              &0.d0, -10.d0, 0.d0, 0.d0, 0.d0, 1.d10, 0.d0, &
              &0.d0, 1.01d0 /)

REAL(dp) :: vx(27) = (/ 0.9d0, -1.d-3, 1.d1, 0.8d0,    &
              &1.d2, 0.8d0, 1.d2, 0.5d0, 0.5d0, 1.d0, 1.d0, &
              &1.d0, 1.d0, 0.1d0, 1.d0, 1.d0, 1.d0, 1.d0,   &
              &1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0,    &
              &1.d0, 1.d2 /)

!!$CHARACTER(4),PARAMETER :: cngd(3) = (/ '---c', 'hang', 'ed v' /), &
!!$                                &dflt(3) = (/ 'nond', 'efau', 'lt v' /)

!.......................................................................

if(iv(1)==0)call dfault(iv,v)
pu=iv(prunit)
iv1=iv(1)
if(iv1/=12)goto 30
if(nn>=n.and.n>=p.and.p>=1)goto 20
iv(1)=16
if(pu>=0)write(pu,10)nn,n,p
10 format(' NL2SOL: Bad nn, n, or p. nn =', i5, ', n =', i5, ', p =', i5/)

goto 300
20 k=iv(21)
call dfault(iv(21:),v(33:))
iv(21)=k
iv(dtype0)=iv(dtype+20)
iv(oldn)=n
iv(oldnn)=nn
iv(oldp)=p
! which(1)=dflt(1)
! which(2)=dflt(2)
! which(3)=dflt(3)
goto 80
30 if(n==iv(oldn).and.nn==iv(oldnn).and.p==iv(oldp))goto 50
iv(1)=17
if(pu>=0)write(pu,40)iv(oldnn),iv(oldn),iv(oldp),nn,n,p
!40 format('0///// (nn,n,p) changed from (', i5, ',', i5, ',', i3,  &
!    &') to (' , i5, ',', i5, ',', i3, ').')
40 format(' NL2SOL: (nn,n,p) changed from (', i5, ',', i5, ',', i3,')'/&
         &' to (' , i5, ',', i5, ',', i3, ').'/)
goto 300

50 if(iv1<=11.and.iv1>=1)goto 70
iv(1)=50
if(pu>=0)write(pu,60)iv1
60 format(' NL2SOL: iv(1) =', i5, ' should be between 0 and 12.')
goto 300

70 continue ! which(1)=cngd(1)
! which(2)=cngd(2)
! which(3)=cngd(3)

80 if(big>tiny)goto 90
tiny=rmdcon(1)
machep=rmdcon(3)
big=rmdcon(6)
vm(12)=machep
vx(12)=big
vm(13)=tiny
vx(13)=big
vm(14)=machep
vm(17)=tiny
vx(17)=big
vm(18)=machep
vx(19)=big
vx(20)=big
vx(21)=big
vm(22)=machep
vx(24)=rmdcon(5)
vm(25)=machep
vm(26)=machep
90 m=0
if(iv(inits)>=0.and.iv(inits)<=2)goto 110
m=18
if(pu>=0)write(pu,100)iv(inits)
100 format (' NL2SOL: iv(25) =', i4, ' should be between 0 and 2.')
110 k=epslon

do i=1,nvdflt
 vk=v(k)
 if(vk>=vm(i).and.vk<=vx(i))goto 130
 m=k
 if(pu>=0)write(pu,120)k,vk,vm(i),vx(i)
 120 format (' NL2SOL: v(',i2,') =',e12.4, ' should be between',e12.4,' and',&
  &e12.4)
 130 k=k+1
enddo

if(iv1==12.and.v(jtinit)>zero)goto 170

! Check jtol values.

jtolp=jtol0+p
do i=jtol1,jtolp
 if(v(i)>zero)cycle
 k=i-jtol0
 if(pu>=0)write(pu,150)k,i,v(i)
 150 format(' NL2SOL: jtol(',i3, ') = v(',i3, ') =',f11.3,&
  &' should be positive.')
 m=i
enddo

170 if(m==0)goto 180
iv(1)=m
goto 300

180 if(pu<0.or.iv(parprt)==0)goto 300
if(iv1/=12.or.iv(inits)==0)goto 200
m=1
write(pu,190)iv(inits)
190 format(' iv(25) =',i3)
200 if(iv(dtype)==iv(dtype0))goto 220
!if(m==0)write(pu,230)which
m=1
write(pu,210)iv(dtype)
!210 format(' dtype..... iv(16) =', i3)
210 format(' iv(16) =', i3)
220 k=epslon
l=parsv1
do i=1,nvdflt
 if(v(k)==v(l))goto 250
 m=1
 write(pu,240)k,v(k)
 240 format (' v(',i2,') =',e12.4)
 250 k=k+1
 l=l+1
enddo
iv(dtype0)=iv(dtype)
call vcopy (nvdflt,v(parsv1:parsv1+nvdflt-1),v(epslon:))
if(iv1/=12)goto 300
if(v(jtinit)>zero)goto 280
jtolp=jtol0+p
write(pu,*)'Initial jtol array:'
write(pu,270)v(jtol1:jtolp)
!270 format('0(initial) jtol array...'/ (' ', 6F12.3))
270 format(6(1x,e12.4))
280 if(v(d0init)>zero)goto 300
k=jtol1+p
l=k+p-1
write(pu,*)'Initial d0 array:'
write(pu,290)v(k:l)
!290 format ('0(initial) d0 array...'/' ', 6F12.3)
290 format (6(1x,e12.4))

300 return
END SUBROUTINE parchk


SUBROUTINE qapply(n,p,j,r,ierr)
INTEGER,INTENT(in) :: n, p, ierr
REAL(dp),INTENT(in) :: j(:,:)
REAL(dp),INTENT(inout) :: r(:)

!     ..................................................................

!     *****purpose.
!     this subroutine applies to r the orthogonal transformations
!     stored in j by qrfact

!     *****parameter description.
!     on input.

!        nn is the row dimension of the matrix j as declared in
!             the calling program dimension statement

!        n is the number of rows of j and the size of the vector r

!        p is the number of columns of j and the size of sigma

!        j contains on and below its diagonal the column vectors
!             u which determine the householder transformations
!             ident - u*u.transpose

!        r is the right hand side vector to which the orthogonal
!             transformations will be applied

!        ierr if non-zero indicates that not all the transformations
!             were successfully determined and only the first
!             abs(ierr) - 1 transformations will be used

!     on output.

!        r has been overwritten by its transformed image

!     *****application and usage restrictions.
!     none

!     *****algorithm notes.
!     the vectors u which determine the householder transformations
!     are normalized so that their 2-norm squared is 2.  the use of
!     these transformations here is in the spirit of (1).

!     *****subroutines and functions called.

!     dotprd - function, returns the inner product of vectors

!     *****references.
!     (1) businger, p. a., and golub, g. h. (1965), linear least squares
!        solutions by householder transformations, numer. math. 7,
!        pp. 269-276.

!     *****history.
!     designed by david m. gay, coded by stephen c. peters (winter 1977)

!     *****general.

!     this subroutine was written in connection with research
!     supported by the national science foundation under grants
!     mcs-7600324, dcr75-10143, 76-14311dss, and mcs76-11989.

!     ..................................................................

INTEGER k,l,nl1
REAL(dp) t

k=p
if(ierr/=0)k=abs(ierr)-1

do l=1,k
 nl1=n-l+1
 t=-dotprd(nl1,j(l:,l),r(l:))
 r(l:n)=r(l:n)+t*j(l:n,l)
enddo

END SUBROUTINE qapply


SUBROUTINE qrfact(m,n,qr,alpha,ipivot,ierr,nopivk,sum)
!------------------------------------------------------------!
! Compute the qr decomposition of the matrix stored in qr.   !
!------------------------------------------------------------!

INTEGER,INTENT(in) :: m, n, nopivk
INTEGER,INTENT(inout) :: ierr, ipivot(:)
REAL(dp),INTENT(inout) :: qr(:,:)
REAL(dp),INTENT(inout) :: alpha(:),sum(:)
INTEGER i,j,jbar,k,k1,minum,mk1
REAL(dp) alphak,beta,qrkk,qrkmax,sigma,temp,rktol1,sumj
REAL(dp),SAVE :: rktol=0._dp,ufeta=0._dp

! dotprd... returns inner product of two vectors.
! rmdcon... returns machine-dependent constants.
! vaxpy... computes scalar times one vector plus another.
! vscopy... sets all elements of a vector to a scalar.
! v2norm... returns the 2-norm of a vector.

REAL(dp),PARAMETER :: one=1.d0,p01=0.01d0,p99=0.99d0,zero=0.d0

!-----------------------------------------------------------------------!

!     *****purpose.

!     this subroutine does a qr-decomposition on the m x n matrix qr,
!        with an optionally modified column pivoting, and returns the
!        upper triangular r-matrix, as well as the orthogonal vectors
!        used in the transformations.

!     *****parameter description.
!     on input.

!        nm must be set to the row dimension of the two dimensional
!             array parameters as declared in the calling program
!             dimension statement.

!        m must be set to the number of rows in the matrix.

!        n must be set to the number of columns in the matrix.

!        qr contains the real rectangular matrix to be decomposed.

!     nopivk is used to control pivotting.  columns 1 through
!        nopivk will remain fixed in position.

!        sum is used for temporary storage for the subroutine.

!     on output.

!        qr contains the non-diagonal elements of the r-matrix
!             in the strict upper triangle. the vectors u, which
!             define the householder transformations   i - u*u-transp,
!             are in the columns of the lower triangle. these vectors u
!             are scaled so that the square of their 2-norm is 2.0.

!        alpha contains the diagonal elements of the r-matrix.

!        ipivot reflects the column pivoting performed on the input
!             matrix to accomplish the decomposition. the j-th
!             element of ipivot gives the column of the original
!             matrix which was pivoted into column j during the
!             decomposition.

!        ierr is set to.
!             0 for normal return,
!             k if no non-zero pivot could be found for the k-th
!                  transformation, or
!             -k for an error exit on the k-th transformation.
!             if an error exit was taken, the first (k - 1)
!             transformations are correct.


!     *****applications and usage restrictions.
!     this may be used when solving linear least-squares problems --
!     see subroutine qr1 of rosepack.  it is called for this purpose
!     by llsqst in the nl2sol (nonlinear least-squares) package.

!     *****algorithm notes.
!     this version of qrfact tries to eliminate the occurrence of
!     underflows during the accumulation of inner products.  rktol1
!     is chosen below so as to insure that discarded terms have no
!     effect on the computed two-norms.

!     adapted from the algol routine solve (1).

!     *****references.
!     (1)     businger,p. and golub,g.h., linear least squares
!     solutions by householder transformations, in wilkinson,j.h.
!     and reinsch,c.(eds.), handbook for automatic computation,
!     volume ii. linear algebra, springer-verlag, 111-118 (1971).
!     prepublished in numer.math. 7, 269-276 (1965).

!     *****history.
!     this amounts to the subroutine qr1 of rosepack with rktol1 used
!     in place of rktol below, with v2norm used to initialize (and
!     sometimes update) the sum array, and with calls on dotprd and
!     vaxpy in place of some loops.

!     *****general.

!     development of this program supported in part by
!     national science foundation grant gj-1154x3 and
!     national science foundation grant dcr75-08802
!     to national bureau of economic research, inc.

!     ..........  ufeta is the smallest positive floating point number
!        s.t. ufeta and -ufeta can both be represented.

!     ..........  rktol is the square root of the relative precision
!        of floating point arithmetic (machep).
!------------------------------------------------------------------------!
if(ufeta>zero)goto 10
ufeta=rmdcon(1)
rktol=rmdcon(4)
10 ierr=0
rktol1=p01*rktol

do j=1,n
 sum(j)=v2norm(m,qr(:,j))
 ipivot(j)=j
enddo

minum=min(m,n)

do k=1,minum
 mk1=m-k+1
! k-th Householder transformation.
 sigma=zero
 jbar=0
! Find largest column sum.
 if(k<=nopivk)goto 50
 DO j=k,n
  if(sigma>=sum(j))cycle
  sigma=sum(j)
  jbar=j
 enddo

 if(jbar == 0) goto 120
 if(jbar == k) goto 50
! Column interchange.
 i=ipivot(k)
 ipivot(k)=ipivot(jbar)
 ipivot(jbar)=i
 sum(jbar)=sum(k)
 sum(k)=sigma

 do i=1,m
  sigma=qr(i,k)
  qr(i,k)=qr(i,jbar)
  qr(i,jbar)=sigma
 enddo
! End of column interchange.
! Second inner product.
 50 qrkmax=zero

 do i=k,m
  if(abs(qr(i,k))>qrkmax)qrkmax=abs(qr(i,k))
 enddo

 if(qrkmax<ufeta)goto 110
 alphak=v2norm(mk1,qr(k:,k))/qrkmax
 sigma=alphak**2

! End second inner product.
 qrkk=qr(k,k)
 if(qrkk>=zero)alphak=-alphak
 alpha(k)=alphak*qrkmax
 beta=qrkmax*sqrt(sigma-(qrkk*alphak/qrkmax))
 qr(k,k)=qrkk-alpha(k)
 qr(k:m,k)=qr(k:m,k)/beta
 k1=k+1
 if(k1>n)cycle

 do j=k1,n
  temp=-dotprd(mk1,qr(k:,k),qr(k:,j))

! Set qr(i,j) = qr(i,j) + temp*qr(i,k), i = k,...,m.

  call vaxpy (mk1, qr(k:k+mk1-1,j), temp, qr(k:,k), qr(k:,j))

  if(k1>m)cycle
  sumj=sum(j)
  if(sumj<ufeta)cycle
  temp=abs(qr(k,j)/sumj)
  if(temp<rktol1)cycle
  if(temp>=p99)goto 80
  sum(j)=sumj*sqrt(one-temp**2)
  cycle
  80 sum(j)=v2norm(m-k,qr(k1:,j))
 enddo
! Eend of k-th Householder transformation.
enddo

goto 150
! Error exit on k-th transformation.
110 ierr=-k
goto 130
! No non-zero acceptable pivot found.
120 ierr=k
130 do i=k,n
 alpha(i)=zero
 if(i>k)call vscopy(i-k,qr(k:k+i-k-1,i),zero)
enddo

150 return

END SUBROUTINE qrfact


FUNCTION reldst(p,d,x,x0) RESULT(fn_val)
!---------------------------------------------------------!
! Compute and return relative difference between x and x0 !
! NL2SOL version 2.2                                      !
!---------------------------------------------------------!
INTEGER,INTENT(in) :: p
REAL(dp),INTENT(in) :: d(:),x(:),x0(:)
REAL(dp) fn_val
INTEGER i
REAL(dp) :: emax,t,xmax,zero=0.d0

emax=zero
xmax=zero
do i=1,p
 t=abs(d(i)*(x(i)-x0(i)))
 if(emax<t)emax=t
 t=d(i)*(abs(x(i))+abs(x0(i)))
 if(xmax<t)xmax=t
enddo
fn_val=zero
if(xmax>zero)fn_val=emax/xmax

END FUNCTION reldst


SUBROUTINE rptmul(func,ipivot,j,p,rd,x,y,z)
!---------------------------------------------------------------------!
!  func = 1... set  y = rmat * (perm**t) * x.                         !
!  func = 2... set  y = perm * (rmat**t) * rmat * (perm**t) * x.      !
!  func = 3... set  y = perm * (rmat**t) x.                           !
!                                                                     !
!  perm = matrix whose i-th col. is the ipivot(i)-th unit vector.     !
!  rmat is the upper triangular matrix whose strict upper triangle    !
!       is stored in  j  and whose diagonal is stored in rd.          !
!  z is a scratch vector.                                             !
!  x and y may share storage.                                         !
!---------------------------------------------------------------------!
INTEGER,INTENT(in) :: func,p,ipivot(:)
REAL(dp),INTENT(in) :: j(:,:),rd(:)
REAL(dp),INTENT(inout) :: x(:),y(:),z(:)
INTEGER i,im1,k,km1
REAL(dp) zk

if(func>2)goto 50

! First set  z = (perm**t) * x

do i=1,p
 k=ipivot(i)
 z(i)=x(k)
enddo

! Now set  y = rmat * z

y(1)=z(1)*rd(1)
if(p<=1)goto 40
do k=2,p
 km1=k-1
 zk=z(k)
 y(1:km1)=y(1:km1)+j(1:km1,k)*zk
 y(k)=zk*rd(k)
enddo

40 if(func<=1)goto 110
goto 70

50 y(1:p)=x(1:p)

! Set z = (rmat**t) * y

70 z(1)=y(1)*rd(1)
if(p==1)goto 90
do i=2,p
 im1=i-1
 z(i)=y(i)*rd(i)+dotprd(im1,j(1:,i),y)
enddo

! Now set  y = perm * z

90 do i=1,p
 k=ipivot(i)
 y(k)=z(i)
enddo

110 return
END SUBROUTINE rptmul


SUBROUTINE slupdt(a,cosmin,p,size,step,u,w,wchmtd,wscale,y)
!--------------------------------------------!
! Update symmetric  a  so that  a * step = y !
! (lower triangle of  a  stored rowwise)     !
!--------------------------------------------!
INTEGER,INTENT(in) :: p
REAL(dp),INTENT(in) :: cosmin,size,step(:),wchmtd(:),y(:)
REAL(dp),INTENT(inout) :: a(:),u(:),w(:)
REAL(dp),INTENT(out) :: wscale
INTEGER i,j,k
REAL(dp) denmin,sdotwm,t,ui,wi
REAL(dp),PARAMETER :: half=0.5d0,one=1.d0,zero=0.d0

sdotwm=dotprd(p,step,wchmtd)
denmin=cosmin*v2norm(p,step)*v2norm(p,wchmtd)
wscale=one
if(denmin/=zero)wscale=min(one,abs(sdotwm/denmin))
t=zero
if(sdotwm/=zero)t=wscale/sdotwm
w(1:p)=t*wchmtd(1:p)
call slvmul(p,u(1:p),a,step)
t=half*(size*dotprd(p,step,u)-dotprd(p,step,y))
u(1:p)=t*w(1:p)+y(1:p)-size*u(1:p)

! Set  a = a + u*(w**t) + w*(u**t)

k=1
do i=1,p
 ui=u(i)
 wi=w(i)
 do j=1,i
  a(k)=size*a(k)+ui*w(j)+wi*u(j)
  k=k+1
 enddo
enddo

END SUBROUTINE slupdt


SUBROUTINE slvmul(p,y,s,x)
!-----------------------------------------------!
! Set  y = s * x,  s = p x p symmetric matrix.  !
! lower triangle of  s  stored rowwise.         !
!-----------------------------------------------!
INTEGER,INTENT(in) :: p
REAL(dp),INTENT(in) :: s(:), x(:)
REAL(dp),INTENT(inout) :: y(:)
INTEGER i,im1,j,k
REAL(dp) xi

j=1
do i=1,p
 y(i)=dotprd(i,s(j:),x)
 j=j+i
enddo

if(p<=1)goto 40
j=1
do i=2,p
 xi=x(i)
 im1=i-1
 j=j+1
 do k=1,im1
  y(k)=y(k)+s(j)*xi
  j=j+1
 enddo
enddo

40 return
END SUBROUTINE slvmul


FUNCTION stopx() RESULT(fn_val)
!-------------------------------------------------------------------!
! This function may serve as the stopx (asynchronous interruption)  !
! function for the NL2SOL (nonlinear least-squares) package at      !
! those installations which do not wish to implement a              !
! dynamic stopx.                                                    !
!                                                                   !
! At installations where the NL2SOL system is used                  !
! interactively, this dummy stopx should be replaced by a           !
! function that returns .true. if and only if the interrupt         !
! (break) key has been pressed since the last call on stopx.        !
!-------------------------------------------------------------------!
LOGICAL fn_val

fn_val=stop_nl2sol
stop_nl2sol=.false.

END FUNCTION stopx


SUBROUTINE vaxpy(p,w,a,x,y)
!-------------------------------------------------------!
! Set w = a*x + y  --  w, x, y = p-vectors, a = scalar  !
!-------------------------------------------------------!
INTEGER,INTENT(in) :: p
REAL(dp),INTENT(in) :: a,x(:),y(:)
REAL(dp),INTENT(inout) :: w(:)

w(1:p)=a*x(1:p)+y(1:p)

END SUBROUTINE vaxpy


SUBROUTINE vcopy (p, y, x)
!------------------------------------------!
! Set y = x, where x and y are p-vectors.  !
!------------------------------------------!
INTEGER,INTENT(in) :: p
REAL(dp),INTENT(in) :: x(:)
REAL(dp),INTENT(inout) :: y(:)

y(1:p)=x(1:p)

END SUBROUTINE vcopy


SUBROUTINE vscopy (p, y, s)
!------------------------------!
! Set p-vector y to scalar s.  !
!------------------------------!
INTEGER,INTENT(in) :: p
REAL(dp),INTENT(in) :: s
REAL(dp),INTENT(inout) :: y(:)

y(1:p)=s

END SUBROUTINE vscopy


FUNCTION v2norm (p, x) RESULT(fn_val)
!----------------------------------------------!
! Return the 2-norm of the p-vector x, taking  !
! care to avoid the most likely underflows.    !
!----------------------------------------------!
INTEGER,INTENT(in) :: p
REAL(dp),INTENT(in) :: x(:)
REAL(dp) fn_val
INTEGER i,j
REAL(dp) r,scale,t,xi

REAL(dp),PARAMETER :: one=1.d0,zero=0.d0
REAL(dp),SAVE :: sqteta=0.d0

if(p>0)goto 10
fn_val=zero
goto 70
10 do i=1,p
 if(x(i)/=zero)goto 30
enddo
fn_val=zero
goto 70

30 scale=abs(x(i))
if(i<p)goto 40
fn_val=scale
goto 70
40 t=one
if(sqteta==zero)sqteta=rmdcon(2)

! sqteta is (slightly larger than) the square root of the
! smallest positive floating point number on the machine.
! the tests involving sqteta are done to prevent underflows.

j=i+1
do i=j,p
 xi=abs(x(i))
 if(xi>scale)goto 50
 r=xi/scale
 if(r>sqteta)t=t+r*r
 cycle
 50 r=scale/xi
 if(r<=sqteta)r=zero
 t=one+t*r*r
 scale=xi
enddo

fn_val=scale*sqrt(t)
70 return
END FUNCTION v2norm


END MODULE toms573


! The following taken from Alan Miller's documentation regarding NL2SOL.
! It seems pgf90 might be the code he refers to.

!** WARNING **
!-------------
!Use this code with care.
!
!Some compilers use `copy in, copy out' where the INTENT of an array
!is IN OUT.   A call such as:
!
!CALL l7vml(p, v(g1:), v(rmat1:), v(qtr1:))
!
!which occurs in routine RN2G, can result in routine l7vml copying the
!whole of the array v from location g1 to one place, the whole of the
!same array from location rmat1 to another place, and the whole of the
!array from location qtr1 to another place.   On exit from l7vml, the
!contents of these temporary workplaces are then written back to array
!v with partial overwriting.
!
!This can be overcome either by replacing the call with something like:
!
!CALL l7vml(p, v(g1:g2), v(rmat1:rmat2), v(qtr1:qtr2))
!
!where g2, rmat2 and qtr2 are the appropriate upper limits to the parts
!of array v which are to be used by routine l7vml, or the large array v
!can be replaced by a number of smaller arrays so that the call becomes:
!
!CALL l7vml(p, g1, rmat1, qtr1)
!
!Either of these requires a very substantial amount of work, as there are
!many instances where different parts of array v are passed to the same
!subroutine.
!
!Try the test programs T_SHORT and T_LONG with your compiler.
!If it crashes, then you cannot use this code.   Most compilers do not
!use `copy in, copy out', so this code works.   At least one major compiler
!does give trouble.   The problem is not with the compiler but with this
!code.
!
! *** 07.2004 *** nl2sol.f90 modified (based on these guidelines) so that
! pgf90 compiler worked.
! Problem seemed occur for arrays of intent OUT mainly, and subroutine dfault
! had to be altered so that arrays iv() and v() are treated as IN OUT, which
! solves a problem with this subroutine when called using
! `call dfault(iv(21:),v(33:))' [from parchk; note there is no possible
! overlapping here].  It looks like the `copy in, copy out' issue is
! somewhat more delicate; probably what pgf does is to store the output array
! back on the original one by overwriting even non-modified elements - which
! is not a bad thing to do, actually. Changes performed according to this idea
! seem to fix the pgf issue: now, all calls to subroutines have explicit upper
! bounds for all arrays whose intent is OUT or INOUT.
