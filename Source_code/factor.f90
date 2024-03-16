      module factors

      contains
      subroutine combination (n_n,r_r,com)      
      implicit none
      integer :: n_n,&
                 r_r,&
                 com,&
                 factor,&
                 fact1,&
                 fact2,&
                 fact3
  
      if (r_r .gt. n_n) then
      write(*,*) "error in combination function"
      end if
        factor = (n_n-r_r)

        call factorial(n_n,fact1) 
        call factorial(r_r,fact2)
        call factorial(factor,fact3)
       com=fact1/(fact2*fact3)        
       end subroutine combination


      subroutine dfactorial (a,dfact)   
      implicit none

      integer :: a,&
                 i,&
                 dfact
      dfact=1
      do i=a,1,-2
        dfact=dfact*i
      end do
      end subroutine dfactorial
   

      subroutine factorial(a,fact)
      implicit none
      integer :: a,&
                 i,&
                 fact
      fact=1
      do i=1,a
          fact =fact*i
      end do
      end subroutine factorial

      end module factors
     
      

           
        









