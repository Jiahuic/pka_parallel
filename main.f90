! This is main drive for compute pka by parallelization
program pka_tabi_parallel
use molecule
use comdata
use bicg
use treecode
use treecode3d_procedures
use MPI_var

include 'mpif.h'

integer i,j,isrc,ierr,is,ie
integer, dimension(:), allocatable :: rcounts,displs
real*8 :: soleng
real*8, dimension(:), allocatable :: sitepotential, sitepartpot
character(100) fhead
integer status(MPI_STATUS_SIZE)

call MPI_INIT(ierr)

call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, num_process, ierr)
ALLOCATE(displs(num_process),STAT=ierr)
ALLOCATE(rcounts(num_process),STAT=ierr)
!??????????????????????????????????????????????????????????????
!PARAMETERS: now can be read from usrdate.in file

!PB equation
!eps0=1.d0;            !the dielectric constant in molecule
!eps1=80.d0;           !the dielectric constant in solvent
!bulk_strength=0.15d0  !ion_strength with units (M)$I=\sum\limits_{i=1}^nc_iz_i^2$

!Treecode
!order=3               !The order of taylor expansion
!maxparnode=500        !maximum particles per leaf
!theta=0.8d0           !MAC, rc/R<MAC, the bigger MAC, the more treecode
open(101,file="usrdata.in")
READ(101,*,IOSTAT = MEOF) fhead, protein
READ(101,*,IOSTAT = MEOF) fhead, den
READ(101,*,IOSTAT = MEOF) fhead, eps0
READ(101,*,IOSTAT = MEOF) fhead, eps1
READ(101,*,IOSTAT = MEOF) fhead, bulk_strength
READ(101,*,IOSTAT = MEOF) fhead, order
READ(101,*,IOSTAT = MEOF) fhead, maxparnode
READ(101,*,IOSTAT = MEOF) fhead, theta
close(101)
!print *,den,eps0,eps1,bulk_strength,order,maxparnode,theta
!???????????????????????????????????????????????????????????????

if (my_rank == 0) then  
  lenfname = len(protein)
  do while (protein(lenfname:lenfname) .eq. ' ')
    lenfname = lenfname - 1
  enddo
  rslt=system('python filenames.py '//protein(1:lenfname))
  nsite = 0
  OPEN(102,FILE=protein(1:lenfname)//".names")
  do
    read(102,*,IOSTAT = MEOF) sitename
  
      if ((MEOF .eq. 0)) then
        nsite=nsite+1
      endif
      IF(MEOF .LT. 0) EXIT
  enddo
  CLOSE(102)
  !nsite = 34
endif
call MPI_Bcast(nsite, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

ALLOCATE(sitelength(nsite),STAT=ierr)
ALLOCATE(sitelist(nsite),STAT=ierr)
if (my_rank == 0) then
  ALLOCATE(sitepotential(nsite),STAT=ierr)
  sitepotential = 0.D0
endif

! my_rank = 0 read the sitation names and boardcast to all
if (my_rank == 0) then
  OPEN(102,FILE=protein(1:lenfname)//".names")
  do i=1,nsite
    read(102,*,IOSTAT = MEOF) sitename
      sitelength(i) = len(sitename)
      do while (sitename(sitelength(i):sitelength(i)) .eq. ' ')
        sitelength(i) = sitelength(i) - 1
      enddo
  
      sitelist(i) = sitename
  enddo
  CLOSE(102)
  rslt=system('rm '//protein(1:lenfname)//".names")
endif

call MPI_Bcast(sitelength, nsite, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
do i=1,nsite
  call MPI_Bcast(sitelist(i), 100, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
enddo

! job assignment
is = my_rank*nsite/num_process+1
ie = (my_rank+1)*nsite/num_process
ALLOCATE(sitepartpot(ie-is+1),STAT=ierr)
sitepartpot = 0.D0

do i = 1,num_process
  displs(i) = (i-1)*nsite/num_process
  rcounts(i) = i*nsite/num_process-displs(i)
!  if (my_rank == 0) then
!    print *,displs(i),rcounts(i)
!  endif
enddo

!print *,'my_rank',is,my_rank,ie-is+1

isrc = 1
do i = is,ie
  fname = sitelist(i)
  soleng = 0.D0
  call tabipb(soleng)
  sitepartpot(isrc) = soleng
  !sitepartpot(isrc) = DBLE(i)
  isrc = isrc+1
enddo
call MPI_GATHERV(sitepartpot,ie-is+1,MPI_REAL8,sitepotential,rcounts,displs,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

if (my_rank == 0) then
  OPEN(102,FILE='free_energy.txt')
  do i=1,nsite
    fname = sitelist(i)
    write(102,*) fname(1:sitelength(i)),sitepotential(i)
  enddo
  CLOSE(102)
  DEALLOCATE(sitepotential)
endif
DEALLOCATE(sitelength)
DEALLOCATE(sitelist)
DEALLOCATE(sitepartpot)
DEALLOCATE(displs)
DEALLOCATE(rcounts)

call MPI_FINALIZE
end program pka_tabi_parallel
