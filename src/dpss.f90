!     spsmooth
!     An extension for mgcv and gam.
!     Copyright (C) 2012 Wesley Burr
!
!     This file is part of the spsmooth package for R.
!
!     The spsmooth package is free software: you can redistribute it 
!     and/or modify it under the terms of the GNU General Public License as 
!     published by the Free Software Foundation, either version 2 of the 
!     License, or any later version.
!
!     The spsmooth package is distributed in the hope that it will be 
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with spsmooth.  If not, see <http://www.gnu.org/licenses/>.
!
!     If you wish to report bugs please contact the author. 
!     Wesley Burr <wburr@mast.queensu.ca>
!     239 Jeffery Hall, Queen's University, Kingston Ontario
!     Canada, K7L 3N6

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     dpss.f90 calculate dpss's using lapack dstebz and dstein
!     using the tridiagonal method.
!      
!     Note the subroutine expects memory for the matrix v and the vector 
!     ev to be allocated by the calling program. In the spsmooth R package, 
!     the allocations occur in dpss.R.
!    
!     Code trivially adapted from dpss.R and dpss.f90 from the 'multitaper'
!     package, also available on CRAN. Code is included in this package
!     to prevent users from having to install multitaper just to get
!     access to these sequences.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Subroutine dpss
!
! Calculates discrete prolate spheroidal sequences using the 
! tridiagonal formulation given in 
!    Spectral Analysis for Physical Applications
!    D.B. Percival and A.T. Walden
!    Cambridge University Press, 1993
!      Chapter 9
! using LAPACK in place of EISPACK. Code also makes use of the 
! David Slepian trick of reducing the problem by splitting the 
! symmetric tridiagonal matrix into two half-sized matrices based
! on parity. This trick was first mentioned in 1977 in a Bell Labs
! Technical Memo.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dpss (n, k, nw, v, ev)

  implicit none
  integer :: n, k,  oddK, evenK, nOdd, nEven, i, j, &
       oddK_2, evenK_2M1, iTest
  double precision :: nw, w, ctpw, pi, twopi, sr2, &
       dlamch, abstol, sqrtsumsq
  parameter(pi=3.141592653589793d0,twopi=2.0d0*pi)
  double precision, pointer :: d(:), e(:), work(:), &
       blockDbleMem(:), evLocal(:)!, vlocal(:,:)
  double precision :: ev(k), v(n,k)

  integer, pointer :: blockIntMem(:)
  logical :: is_evenN
  character :: cmach

  interface
     subroutine tridiagMatrixEigen(n, k, d, e, v, ldv, ev, &
          abstol, blockIntMem, work)
       implicit none
       character :: range, order, cmach
       integer :: nsplit, m, info
       integer, target :: blockIntMem(5*n+k)
       integer,  pointer :: iblock(:), isplit(:), iwork(:), ifail(:)
       integer :: k, ldv, n, il
       double precision :: vl, vu, abstol
       double precision :: d(n), e(n-1), v(ldv,k),  ev(n), work(5*n)
     end subroutine tridiagMatrixEigen
  end interface
  
  w = nw/dble(n)
  ctpw =  dcos(twopi*w)
  oddK = k/2
  evenK = k - oddK
  sr2 = dsqrt(2.0d0)

  nOdd = n/2
  nEven = n - nOdd
  
  ! Allocate memory and use pointers to reduce malloc calls.
  ! The values 5 and 8 provide the required memory space for the two 
  ! LAPACK calls, see the documentation for dstebz and dstein.
  ! Memory used in this procedure is allocated in the following two calls 
  ! and then appropriate blocks of memory are accessed using pointers. 
  !
  ! Note: memory space for the matrix used by the calling program must 
  ! be allocated in the calling program!

  allocate(blockIntMem(5*nEven+k))
  allocate(blockDbleMem(8*nEven-1))

  ! Allocate space for even and odd parity blocks
  d => blockDbleMem(1:nEven)
  e => blockDbleMem((nEven+1):(2*nEven-1))
  work => blockDbleMem((2*nEven):(7*nEven-1))
  evLocal => blockDbleMem((7*nEven):(8*nEven-1))

  i = 0
  d = (/ (((n-1-2*i) / 2.0d0)**2 * ctpw, i=0, nEven-1) /)

  ! Convert n and i to double before the multiplication; if this
  ! is not done, very large choices of N cause overflow.
  e = (/ ((dble(i) * dble(n - i)) / 2.0d0, i = 1, nEven-1) /)

  is_evenN  = (modulo(n, 2) .eq. 0)
  
  ! Ensure integer values are converted to double before multiplication
  ! to avoid overflow.
  if(is_evenN)  then 
     d(nEven) = ((n+1-2*nEven)/2.0d0)**2 * ctpw + dble(nEven) * dble(nOdd)/2.0d0 
  else 
     e(nEven -1) =  e(nEven -1) * sr2
  end if
  
  cmach = 'S'
  abstol = 2.0d0 * dlamch(cmach)
  
  ! Set ldv (ldz) to 2*n to force dstein to skip odd columns which
  ! will be used for odd eigenvectors. The matrix v is considered to 
  ! be of a different shape in the call to tridiagMatrixEigen2
  call tridiagMatrixEigen(nEven, evenK, d, e, v, 2*n, evlocal, &
       abstol, blockIntMem, work)
  
  if(.not. is_evenN) then
    do i = 1, k, 2
       v(nEven, i) =  v(nEven, i) * sr2
    end do
  end if
  
  do i = 1, k, 2
    v((nEven+1):n ,i) = v(nOdd:1:-1, i)
  end do
  
  ! Reorder eigenvalues and eigenfunction columns for even
  evenK_2M1 = evenK * 2 - 1
  j = 1
  ev((/ (i, i=evenK_2M1, 1, -2) /)) = evLocal((/ (j, j=1, evenK, 1) /)) 
  v(:, (/ (i, i=1, evenK_2M1, 2) /)) =  & 
       v(:,(/ (j, j=evenK_2M1, 1, -2) /))
  
  if(k > 1)  then
     
     ! Odd eigenfunctions (if any) follow a similar procedure to the 
     ! even functions. If n is odd, nOdd < nEven. Again, ensure that
     ! integer values are converted to double before multiplication
     ! to avoid overflow.

     d(1:nOdd) = (/ (((n-1-2*i) / 2.0d0)**2 * ctpw, i=0, nOdd-1) /)
     
     ! Convert n and i to double before the multiplication. 
     e(1:(nOdd-1)) = (/ ((dble(i) * dble(n - i)) / 2.0d0, i = 1, nOdd-1) /)
     
     if(is_evenN) then 
        ! convert n and i to double before the multiplication
        d(nOdd) = ((n+1-2*nEven)/2.0d0)**2 * ctpw - dble(nEven) * dble(nOdd)/2.0d0  
     end if
     
     call tridiagMatrixEigen(nOdd, oddK, d, e, v(1,2), 2*n, evlocal, &
          abstol, blockIntMem, work)
     
     if(.not. is_evenN) then
        do i= 2, k, 2
           v(nEven,i) = 0.0d0
        end do
     end if
     
     do i = 2, k, 2
        v((nEven+1):n ,i) = - v(nOdd:1:-1, i)
     end do
     
     ! Reorder eigenvalues and eigenfunction columns for odd parity.
     oddK_2 = oddK * 2
     ev((/ (i, i=oddK_2, 2, -2) /)) = evLocal((/ (j, j=1, oddK, 1) /))
     v(:, (/ (i, i=2, oddK_2, 2) /)) =  & 
          v(:,(/ (j, j=oddK_2, 2, -2) /))
     
  end if
  
  ! Normalize the tapers.
  iTest = 0
  if (mod(n,2) .eq. 0) then
     iTest = n/2 +1
  else 
     iTest = n/2 +2
  end if

  ! This sets the polarity to that of Slepians 1978 paper
  !    Prolate Spheroidal Wave Functions, Fourier Analysis, and
  !    Uncertainty -- V: The Discrete Case
  !    David Slepian
  !    Bell Systems Technical Journal, 1978
  !    http://www.alcatel-lucent.com/bstj/vol57-1978/articles/bstj57-5-1371.pdf
  ! which sets the sequences to have positive slope at t=0, where
  ! t = -(N-1)/2, ..., -1, 0, 1, ..., (N-1)/2.

  do j = 1 , k
     sqrtsumsq = dsqrt(sum( v(:,j)**2 ))
     
     ! set polarity to Slepian 78 
     ! differs from Percival and Walden
     ! dpss slope up at centre, this  agrees with 
     ! Thomson 82
     if(v(iTest,j) < 0.0d0) then
        sqrtsumsq = -1.0d0 * sqrtsumsq
     end if
     v(:,j) = v(:,j) / sqrtsumsq
  end do

  nullify(d,e,work,evLocal)
  deallocate(blockDbleMem)
  deallocate(blockIntMem)
  
end subroutine dpss


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Subroutine tridiagMatrixEigen
!
! Wrapper for LAPACK routines dstein and dstebz.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine tridiagMatrixEigen(n, k, d, e, v, ldv, ev, &
     abstol, blockIntMem, work)

  ! Assumes memory is allocated by the calling subroutine.
  implicit none
  
  character :: range, order, cmach
  integer :: nsplit, m, info
  integer, target :: blockIntMem(5*n+k)
  integer,  pointer :: iblock(:), isplit(:), iwork(:), ifail(:)
  integer :: k, ldv, n, il
  double precision :: vl, vu, abstol
  double precision :: d(n), e(n-1), v(ldv,k),  ev(n), work(5*n)

  range = 'I'
  order = 'E'
  cmach = 'S'

  il = n - k + 1
  m = k
     
  iblock => blockIntMem(1:n)
  isplit => blockIntMem((n+1):(2*n))
  iwork => blockIntMem((2*n+1):(5*n))
  ifail => blockIntMem((5*n+1):(5*n+k))

  call dstebz(range, order, n, vl, vu, il, n, &   
       abstol, d, e, m, nsplit, ev, iblock, isplit, &
       work, iwork, info)

  call dstein(n, d, e, m, ev, iblock, isplit, v, ldv, &
           work, iwork, ifail, info)

  nullify(iblock, isplit, iwork, ifail)

end subroutine tridiagMatrixEigen
