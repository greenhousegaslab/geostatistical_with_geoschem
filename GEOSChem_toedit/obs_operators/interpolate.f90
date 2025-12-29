    module interpolate
!
!  SYNOPSIS
!     use interpolate
!
!  DESCRIPTION
!    INTERPOLATE contains routines for binary lookup of
!    tables of unique values in ascending and descending order and
!    for index interpolation.
!
!  ROUTINES IMPLEMENTED
!    pure function index_interp(val,table,n)
!    pure subroutine lookup(val,table,n,ibracket)
!    pure subroutine lookdown(val,table,n,ibracket)
!  
!  REQUIRED PACKAGES
!    None.
!
!  MODIFICATION HISTORY
!    01/12/05 - first implementation.
!
!  PROGRAMMER
!    Eric Kostelich, Dept. of Mathematics and Statistics,
!    Arizona State University, Tempe, AZ 85287-1804, kostelich@asu.edu
!
    contains
! ------------------------------------------------------------------- 
    pure real function index_interp(val,table,n)
!  INDEX_INTERP - Given a table of length N in sorted order and a value VAL,
!  apply linear interpolation to find the "subscript" R such TABLE(R)=VAL.
!  For example, if VAL is halfway between TABLE(J) and TABLE(J+1), then
!  INDEX_INTERP returns the floating-point number J+0.5.
!
    implicit none
    integer,intent(in)::n
    real,intent(in)::val,table(n)
!
!  Local variables
!
    integer::b(2)
!
    if(n.le.1) then
       index_interp=real(n)
       return
    endif
!
    if(table(1).lt.table(2)) then  ! assume ascending order
       call lookup(val,table(1:n),n,b)
    else
       call lookdown(val,table(1:n),n,b)
    endif
    if(b(1).eq.0) then
       index_interp=merge(1.0,0.0,val.eq.table(b(2)))
    else if(b(2).gt.n) then
       index_interp=real(n+1)
    else
       index_interp=b(1)+(val-table(b(1)))/(table(b(2))-table(b(1)))
    endif
    return
    end function index_interp
! ------------------------------------------------------------------- 
    pure subroutine lookup(val,table,n,ibracket)
!  LOOKUP - does binary search for the real value VAL in the real
!  table TABLE of length N, which is assumed sorted in ascending order
!  (this is not checked).  LOOKUP returns the 2-vector IBRACKET such that
!  TABLE(IBRACKET(1)) < VAL <= TABLE(IBRACKET(2)).
!  If VAL < TABLE(1) then IBRACKET = [0,1];
!  if VAL > TABLE(N) then IBRACKET = [N,N+1].
!
    implicit none
    integer,intent(in)::n
    real,intent(in)::val,table(n)
    integer,intent(out)::ibracket(2)
!
!  Local variables
!
    integer::lo,hi,mid
!
    if(val.le.table(1)) then
       ibracket=(/0,1/)
    else if(val.gt.table(n)) then
       ibracket=(/n,n+1/)
    else
       lo=1
       hi=n
       do while(hi-lo.gt.1)
          mid=(hi+lo)/2
          if(val.le.table(mid)) then
             hi=mid
          else
             lo=mid
          endif
       enddo
       ibracket=(/lo,hi/)
    endif
    return
    end subroutine lookup
! ------------------------------------------------------------------- 
    pure subroutine lookdown(val,table,n,ibracket)
!  LOOKDOWN - does binary search for the real value VAL in the real
!  table TABLE of length N, which is assumed sorted in descending order
!  (this is not checked).  LOOKDOWN returns the 2-vector IBRACKET such that
!  TABLE(IBRACKET(1)) > VAL >= TABLE(IBRACKET(2)).
!  If VAL > TABLE(1) then IBRACKET = [0,1];
!  if VAL < TABLE(N) then IBRACKET = [N,N+1].
!
    implicit none
    integer,intent(in)::n
    real,intent(in)::val,table(n)
    integer,intent(out)::ibracket(2)
!
!  Local variables
!
    integer::lo,hi,mid
!
    if(val.ge.table(1)) then
       ibracket=(/0,1/)
    else if(val.lt.table(n)) then
       ibracket=(/n,n+1/)
    else
       lo=1
       hi=n
       do while(hi-lo.gt.1)
          mid=(hi+lo)/2
          if(val.ge.table(mid)) then
             hi=mid
          else
             lo=mid
          endif
       enddo
       ibracket=(/lo,hi/)
    endif
    return
    end subroutine lookdown
! ------------------------------------------------------------------- 
    end module interpolate
