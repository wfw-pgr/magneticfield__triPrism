module bTriPrismMod
contains
  
  
  ! ====================================================== !
  ! === calculate magnetic field from prism            === !
  ! ====================================================== !
  subroutine magneticField__triPrism( vertex, bfdpos, mvec, zp, bfield, depth, nDiv_z )
    implicit none
    integer         , parameter   :: nVert=3, dim=3
    integer         , parameter   :: xc_=1, yc_=2, zc_=3, S_=4
    integer         , parameter   :: v1_=1, v2_=2, v3_=3, lo_=1, hi_=2
    integer         , intent(in)  :: depth, nDiv_z
    double precision, intent(in)  :: zp(2)
    double precision, intent(in)  :: mvec(dim), bfdpos(dim), vertex(dim,nVert)
    double precision, intent(out) :: bfield
    double precision, allocatable :: subdivE(:,:,:), triangles(:,:)
    integer                       :: is, iz, nSubE
    double precision              :: xm, ym, zm, delz, zBase, volu

    nSubE    =     4**( depth )
    allocate( subdivE(dim,nVert,nSubE), triangles(4,nSubE) )
    
    ! ------------------------------------------------------ !
    ! --- [1] subdivision of a trianle ( recursive  )    --- !
    ! ------------------------------------------------------ !
    call single__subdivide( vertex(:,v1_), vertex(:,v2_), vertex(:,v3_), subdivE, depth, nSubE )
    
    ! ------------------------------------------------------ !
    ! --- [2] inspect triangle characteristic            --- !
    ! ------------------------------------------------------ !
    call inspect__triangles_Area_CoGs( subdivE, triangles, nSubE )
    
    ! ------------------------------------------------------ !
    ! --- [3] calculate magnetic moment bfield           --- !
    ! ------------------------------------------------------ !
    bfield  =  0.d0
    delz    = ( zp(hi_) - zp(lo_) ) / dble( nDiv_z )
    zBase   =   zp(lo_) + 0.5d0*delz
    do is=1, nSubE
       xm   = triangles(xc_,is)
       ym   = triangles(yc_,is)
       volu = delz * triangles(S_,is)
       do iz=1, nDiv_z
          zm      =  zBase + delz*dble( iz-1 )
          bfield  = bfield + volu*( + MagneticMoment( xm, ym,  zm, bfdpos, mvec ) &
               &                    + MagneticMoment( xm, ym, -zm, bfdpos, mvec ) )
       enddo
    enddo

    ! ------------------------------------------------------ !
    ! --- [4] deallocation                               --- !
    ! ------------------------------------------------------ !
    deallocate( subdivE, triangles )
    
    return
  end subroutine magneticField__triPrism


  ! ====================================================== !
  ! === calculate center of gravity & area of triangle === !
  ! ====================================================== !
  subroutine inspect__triangles_Area_CoGs( subdivE, triangles, nSubE )
    implicit none
    integer         , parameter   :: dim=3, nVert=3
    integer         , intent(in)  :: nSubE
    double precision, intent(in)  :: subdivE(dim,nVert,nSubE)
    double precision, intent(out) :: triangles(4,nSubE)
    integer                       :: is, iv, ic
    double precision              :: v21(2), v31(2)
    double precision, parameter   :: onethird = 1.d0/3.d0
    integer         , parameter   :: xc_=1, yc_=2, zc_=3, S_=4
    integer         , parameter   :: v1_=1, v2_=2, v3_=3

    ! ------------------------------------------------------ !
    ! --- [1] calculate Center of Gravity for triangle   --- !
    ! ------------------------------------------------------ !
    triangles(:,:) = 0.d0
    do is=1, nSubE
       do iv=1, nVert
          do ic=xc_, zc_
             triangles(ic,is) = triangles(ic,is) + subdivE(ic,iv,is)
          enddo
       enddo
    enddo
    triangles(xc_:zc_,:) = onethird * triangles(xc_:zc_,:)

    ! ------------------------------------------------------ !
    ! --- [2] calculate Area for triangle                --- !
    ! ------------------------------------------------------ !
    do is=1, nSubE
       ! -- assume that   z=const.   for all 3 nodes :: 0th-element -- !
       v21(:)           = subdivE(xc_:yc_,v2_,is) - subdivE(xc_:yc_,v1_,is)
       v31(:)           = subdivE(xc_:yc_,v3_,is) - subdivE(xc_:yc_,v1_,is)
       triangles(S_,is) = 0.5d0 * abs( v21(xc_)*v31(yc_) - v21(yc_)*v31(xc_) )
    enddo

    return
  end subroutine inspect__triangles_Area_CoGs

  
  ! ====================================================== !
  ! === recursive subroutine to refine triangles       === !
  ! ====================================================== !
  recursive subroutine single__subdivide( node1, node2, node3, node_ret, depth, node_len )
    implicit none
    integer         , parameter     :: dim=3, nVert=3
    integer         , intent(in)    :: depth, node_len
    double precision, intent(inout) :: node_ret(dim,nVert,node_len)
    double precision, intent(in)    :: node1 (3), node2 (3), node3 (3)
    double precision                :: node12(3), node23(3), node31(3)
    integer                         :: sect_len
    integer         , parameter     :: v1_=1, v2_=2, v3_=3
    
    if      ( depth.eq.0 ) then
       node_ret(:,v1_,1) = node1(:)
       node_ret(:,v2_,1) = node2(:)
       node_ret(:,v3_,1) = node3(:)
       
    else if ( depth.ge.1 ) then
       node12(:) = 0.5d0 * ( node1(:) + node2(:) )
       node23(:) = 0.5d0 * ( node2(:) + node3(:) )
       node31(:) = 0.5d0 * ( node3(:) + node1(:) )
       sect_len  = node_len / 4
       call single__subdivide( node1 , node12, node31, node_ret(:,:,(           1):(1*sect_len) ), depth-1, sect_len )
       call single__subdivide( node2 , node23, node12, node_ret(:,:,(  sect_len+1):(2*sect_len) ), depth-1, sect_len )
       call single__subdivide( node3 , node31, node23, node_ret(:,:,(2*sect_len+1):(3*sect_len) ), depth-1, sect_len )
       call single__subdivide( node12, node23, node31, node_ret(:,:,(3*sect_len+1):(4*sect_len) ), depth-1, sect_len )
    endif
    
    return
  end subroutine single__subdivide

  
  ! ====================================================== !
  ! === calculate MagnetiMoment from position          === !
  ! ====================================================== !
  function MagneticMoment( xm, ym, zm, bfdpos, mvec )
    implicit none
    double precision, intent(in) :: xm, ym, zm, bfdpos(3), mvec(3)
    double precision             :: rabs, mdotr, rvec3Inv, rvec(3), rhat(3)
    double precision             :: MagneticMoment
    double precision, parameter  :: fourpi = 16.d0*atan(1.d0)
    integer         , parameter  :: x_=1, y_=2, z_=3

    rvec(x_)          = bfdpos(x_) - xm
    rvec(y_)          = bfdpos(y_) - ym
    rvec(z_)          = bfdpos(z_) - zm
    rabs              = sqrt( rvec(1)**2 + rvec(2)**2 + rvec(3)**2 )
    rhat(:)           = rvec(:) / rabs
    rvec3Inv          = 1.d0 / ( fourpi * rabs**3 )
    mdotr             = mvec(1)*rhat(1) + mvec(2)*rhat(2) + mvec(3)*rhat(3)
    MagneticMoment    = rvec3Inv * ( 3.d0 * mdotr * rhat(z_) - mvec(z_) )
    return
  end function MagneticMoment
  
  
end module bTriPrismMod


