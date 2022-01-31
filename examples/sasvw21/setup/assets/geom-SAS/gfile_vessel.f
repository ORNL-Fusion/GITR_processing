      program gfile_vessel
      implicit none

      integer i,j,nw,nh,nlim,nbbs
      real gr1
      real, dimension(1:1000):: xlim,ylim
      character(len=10) gc1,gc2,gc3,gc4,gc5,gc6
     
      open(unit=7,file='gfile',status='old')
      open(unit=8,file='vvfile',status='unknown')

      read(7,2000) gc1,gc2,gc3,gc4,gc5,gc6,i,nw,nh
      read(7,2020) gr1
      read(7,2020) gr1
      read(7,2020) gr1
      read(7,2020) gr1
      read(7,2020) (gr1,i=1,nw)
      read(7,2020) (gr1,i=1,nw)
      read(7,2020) (gr1,i=1,nw)
      read(7,2020) (gr1,i=1,nw)
      read(7,2020) ((gr1,i=1,nw),j=1,nh)
      read(7,2020) (gr1,i=1,nw)
      read(7,2022) nbbs,nlim
      read(7,2020) (gr1,gr1,i=1,nbbs)
      read(7,2020) (xlim(i),ylim(i),i=1,nlim)
      do i=1,nlim
         write(8,*) xlim(i)*1000.,ylim(i)*1000.
      end do

 2000 format (6a8,3i4)
 2020 format (5e16.9)
 2022 format (2i5)
      end
