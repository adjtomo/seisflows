program kernels_sum

implicit none

integer, parameter :: NMAX = 9999999

integer :: isrc, NSRC, NLOCAL, j
double precision, allocatable, dimension(:) :: kernel1, kernel2, kernel3
double precision, allocatable, dimension(:) :: kernel_sum1, kernel_sum2, kernel_sum3
double precision, allocatable, dimension(:,:) :: coord
character(len=150) :: dir, filename, arg
integer :: ios
double precision :: dummy1, dummy2, dummy3, dummy4, dummy5
character(len=150) :: dummy_character

! read input arguments
if( iargc() .ne. 2 ) stop 'USAGE: kernels_sum NSRC DIR'
j=1; call getarg(j, arg); read(arg,*,iostat=ios) NSRC;
if (ios /= 0) stop 'Error reading NSRC'
j=2; call getarg(j, dir);
if (ios /= 0) stop 'Error reading DIR'

! dynamically allocate arrays
filename = trim(dir) // '/' // '000000'
open(unit=3,file=filename,status='old',action='read')
NLOCAL = 0
do j=1,NMAX
   read(3,*,iostat=ios) dummy1, dummy2, dummy3, dummy4, dummy5
   if (ios /= 0) exit
   NLOCAL=NLOCAL+1
enddo
close(3)

allocate( coord(2,NLOCAL) )
allocate( kernel1(NLOCAL) )
allocate( kernel2(NLOCAL) )
allocate( kernel3(NLOCAL) )
allocate( kernel_sum1(NLOCAL) )
allocate( kernel_sum2(NLOCAL) )
allocate( kernel_sum3(NLOCAL) )

kernel1(:) = 0.0d0
kernel2(:) = 0.0d0
kernel3(:) = 0.0d0
kernel_sum1(:) = 0.0d0
kernel_sum2(:) = 0.0d0
kernel_sum3(:) = 0.0d0

! read first kernel
write(dummy_character,'(i6.6)') 0
filename = trim(dir) // '/' // trim(ADJUSTL(dummy_character))
open(unit=3,file=filename,status='old',action='read')
do j = 1,nlocal
   read(3,*) coord(1,j), coord(2,j), kernel1(j), kernel2(j), kernel3(j)
enddo
close(3)
kernel_sum1 = kernel1
kernel_sum2 = kernel2
kernel_sum3 = kernel3

! read subsequent kernels
do isrc = 2, NSRC
   write(dummy_character,'(i6.6)') isrc-1
   filename = trim(dir) // '/' // trim(ADJUSTL(dummy_character))
   open(unit=3,file=filename,status='old',action='read')
   do j = 1,nlocal
      read(3,*) dummy1, dummy2, kernel1(j), kernel2(j), kernel3(j)
   enddo
   close(3)
   kernel_sum1(:)=kernel_sum1(:)+kernel1(:)
   kernel_sum2(:)=kernel_sum2(:)+kernel2(:)
   kernel_sum3(:)=kernel_sum3(:)+kernel3(:)
enddo

! save result
filename = trim(dir) // '/sum'
open(unit=4,file=filename,status='unknown',action='write')
do j = 1, NLOCAL
      write(4,'(5e11.3)') coord(1,j),coord(2,j),kernel_sum1(j),kernel_sum2(j),kernel_sum3(j)
enddo
close(4)

end program kernels_sum

