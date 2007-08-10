program pargbcyl
use mpi
use mpicelllist
use io
use nrtype
use particle
use particlearray
implicit none
character(len=30) :: statefile
integer :: N, Nrelax, Nprod, Nratio, anchor, voltyp, seed, allign, debug, &
          &domainw 
real(dp),dimension(:), pointer :: xs,ys,zs,uxs,uys,uzs
real(dp) :: rcyl,hcyl, T, pres, Kw, epses, eps0, rsphere, spmyy, epssphere, &
           &sigma0, siges, cutoff, maxdr 
logical,dimension(:),pointer :: rods
type(particledat), dimension(:), pointer :: parray,myparray
integer :: rc, mynp

call MPI_INIT(rc)
call ReadParams(statefile,Nrelax,Nprod,Nratio,T,pres,anchor,voltyp,Kw,&
               &seed,epses,eps0,rsphere,spmyy,epssphere,sigma0,siges,&
               &allign,debug,cutoff,maxdr,domainw) 
call readState(statefile,N,xs,ys,zs,uxs,uys,uzs,rods,rcyl,hcyl)
call xyzToParticle(N,xs,ys,zs,uxs,uys,uzs,rods,parray)
call initCellList(parray,N,domainw,cutoff,maxdr,-rcyl,rcyl,-rcyl,rcyl,rcyl )
call myParticles(myparray,mynp)
!write(*,*) mynp
call writestateMPI(rcyl,hcyl,N,N,0,myparray,mynp)
call MPI_BARRIER(MPI_COMM_WORLD,rc)
call MPI_FINALIZE(rc)
end program pargbcyl
