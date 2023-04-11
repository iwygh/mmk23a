program midplane
implicit none
character*20,allocatable :: name_list(:)
real*8,allocatable :: vtx_list(:),vty_list(:),vtz_list(:)
real*8 vtx,vty,vtz,nx,ny,nz,theta_mid_deg,phi_mid_deg
integer iobj,nobj,line_count
character*50 infile,outfile,dummytext
infile = 'vtvec_clones.csv'
outfile = 'vm17_mean_plane_clones_noclones.csv'
infile = trim(adjustl(infile))
outfile = trim(adjustl(outfile))
call get_line_count(infile,line_count)
nobj = line_count-1 ! one-line header
if (allocated(vtx_list)) then
    deallocate(vtx_list)
end if
if (allocated(vty_list)) then
    deallocate(vty_list)
end if
if (allocated(vtz_list)) then
    deallocate(vtz_list)
end if
allocate( vtx_list(1:nobj),vty_list(1:nobj),vtz_list(1:nobj) )
open(99,file=infile)
read(99,*) ! skip first line (header)
do iobj = 1,nobj
    read(99,*) dummytext,vtx,vty,vtz
    vtx_list(iobj) = vtx
    vty_list(iobj) = vty
    vtz_list(iobj) = vtz
end do
close(99)
open(90,file=outfile)
write(90,*) 'theta_mid_deg ','phi_mid_deg '
call fit_plane_vxyz(nobj,vtx_list,vty_list,vtz_list, & ! inputs
                & theta_mid_deg,phi_mid_deg) ! outputs
write(90,*) theta_mid_deg,phi_mid_deg
write(*,*) theta_mid_deg,phi_mid_deg
infile = 'vtvec_noclones.csv'
infile = trim(adjustl(infile))
call get_line_count(infile,line_count)
nobj = line_count-1 ! one-line header
if (allocated(vtx_list)) then
    deallocate(vtx_list)
end if
if (allocated(vty_list)) then
    deallocate(vty_list)
end if
if (allocated(vtz_list)) then
    deallocate(vtz_list)
end if
allocate( vtx_list(1:nobj),vty_list(1:nobj),vtz_list(1:nobj) )
open(99,file=infile)
read(99,*) ! skip first line (header)
do iobj = 1,nobj
    read(99,*) dummytext,vtx,vty,vtz
    vtx_list(iobj) = vtx
    vty_list(iobj) = vty
    vtz_list(iobj) = vtz
end do
close(99)
call fit_plane_vxyz(nobj,vtx_list,vty_list,vtz_list, & ! inputs
                & theta_mid_deg,phi_mid_deg) ! outputs
write(90,*) theta_mid_deg,phi_mid_deg
close(90)
write(*,*) theta_mid_deg,phi_mid_deg
end program




subroutine fit_plane_vxyz(nobj,vx_list,vy_list,vz_list, & ! inputs
                    & theta_mid_deg,phi_mid_deg) ! outputs
implicit none
integer iobj,nobj
real*8 thetamin_deg,thetamax_deg,dtheta_deg,phimin_deg,phimax_deg,dphi_deg
real*8 thetamin,thetamax,dtheta,phimin,phimax,dphi,theta_mid_deg,phi_mid_deg
real*8 sm,smmin,theta_mid,phi_mid,delta,nx,ny,nz,theta,phi
real*8 vx_list(1:nobj),vy_list(1:nobj),vz_list(1:nobj)
real*8, parameter::pi = 2.0d0*dasin(1.0d0)
thetamin_deg = 0.02d0
thetamax_deg = 40.0d0
dtheta_deg = 0.01d0
phimin_deg = 0.0d0
phimax_deg = 359.98d0
dphi_deg = 0.01d0
thetamin = thetamin_deg * pi / 180.0d0
thetamax = thetamax_deg * pi / 180.0d0
dtheta = dtheta_deg * pi / 180.0d0
phimin = phimin_deg * pi / 180.0d0
phimax = phimax_deg * pi / 180.0d0
dphi = dphi_deg * pi / 180.0d0
smmin = 1.0d9
theta_mid = 1000.0d0
phi_mid = 1000.0d0
theta = thetamin
phi = phimin
do while (theta .lt. thetamax)
    do while (phi .lt. phimax)
        nx = dsin(theta) * dcos(phi)
        ny = dsin(theta) * dsin(phi)
        nz = dcos(theta)
        sm = 0.0d0
        do iobj = 1,nobj
            delta = vx_list(iobj)*nx + vy_list(iobj)*ny + vz_list(iobj)
            sm = sm + dabs(delta)
        end do
        if (sm .lt. smmin) then
            smmin = sm
            theta_mid = theta
            phi_mid = phi
        end if
        phi = phi + dphi
    end do
    theta = theta + dtheta
    phi = phimin
end do
theta_mid_deg = theta_mid * 180.0d0 / pi
phi_mid_deg = phi_mid * 180.0d0 / pi
return
end subroutine




subroutine get_line_count(infile,& ! inputs
                  & line_count ) ! outputs
implicit none
integer line_count,io
character(len=*) infile
open(99,file=infile,iostat=io)
if (io/=0) stop 'get_line_count: Cannot open file!'
line_count = 0
do
    read(99,*,iostat=io)
    if (io/=0) exit
    line_count = line_count + 1
end do
close(99)
return
end subroutine
