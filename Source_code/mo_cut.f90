          program read_data
          implicit none
          real,allocatable                  :: amatin(:,:),amat(:,:)
          integer                           :: istat
          integer                           :: ii,jj,kk,ll,mm,nn,oo
          integer                           :: mlp,ndum,ndim
          character(len=25)                 :: file_in,file_out


          write(*,'(1x,a)', advance="no") "enter name of input file:"
          read *, file_in
          
          write(*,'(1x,a)', advance="no") "enter name of output file:"
          read *, file_out

          open(10,file=file_in,status="old",action="read",iostat=istat)
          open(20,file=file_out,status="unknown",iostat=istat)

          write(*,'(1x,a)', advance="no") "enter size of the matrix:"
          read *, ndim



          allocate (amatin(ndim,ndim), amat(ndim,ndim), stat=istat)
          if (istat /= 0) stop "*****allocation failed*****"


               do ii=1,ndim 
               read(10,*)  amatin(ii,:)
               end do
                  
 
                               
           do ii=1,ndim 
              amat(:,ii)= amatin(ii,:)
           end do

           do jj=1,ndim
             write(20,400) amat(jj,:)
           end do


  400 format(5e15.7)
      close(10)
      close(20)

                end program  read_data

