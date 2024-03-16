       module basis_fun
       contains
       subroutine position_wav(l,m,n,alpha,x,y,z,r,theta,phi,xa,ya,za,res,grad_x,grad_y,grad_z)
       use factors
       implicit none
       real*8 , parameter  ::    zero=0.d0,&
                                 pi= 3.14159265358979323846264338d0
 

       integer  ::    l,&
                      m,&
                      n,&
                      dfact1,&
                      dfact2,&
                      dfact3

       real*8   ::    x,&
                      y,&                                            !!-------------POSITION-WEB
                      z,&
                      r,&
                      theta,&
                      phi,&
                      xa,&
                      ya,&
                      za,&
                      alpha
       

        real*8  ::   norm,&
                     x_a,&
                     y_a,&
                     z_a,&
                     p,&
                     np,&
                     dmp,&
                     res,&
                     diff_x,&
                     diff_y,&
                     diff_z,&
                     grad_x,&
                     grad_y,&
                     grad_z,&
                     r_r
          
           p=(l+m+n)
          ! write(*,*) p 
           np=(4*alpha)**(p)
            !write(*,*) np
           call dfactorial(((2*l)-1),dfact1)
           call dfactorial(((2*m)-1),dfact2)
           call dfactorial(((2*n)-1),dfact3)
           dmp=dfact1*dfact2*dfact3
           !write(*,*) dmp 
           norm =((2.d0*alpha/pi)**0.75d0)*((np/dmp)**0.5d0)
           !!write(*,*) norm,alpha
           x_a=x-xa 
           y_a=y-ya
           z_a=z-za
           !write(*,*) x,y,z,r

          res=norm*(x_a**l)*(y_a**m)*(z_a**n)*&
          & exp(-alpha*x_a*x_a)*exp(-alpha*y_a*y_a)*exp(-alpha*z_a*z_a)

      !!diff_x=l+m+n/r-(2.d0*alpha*r)
      !!diff_y=((((l+m)*cos(theta))/sin(theta))-((n*sin(theta))/cos(theta))) !! FOR SP.GRID 
      !!diff_z=(((m*cos(phi))/sin(phi))-((l*sin(phi))/cos(phi)))

        diff_x=(l/x_a)-(2.0d0*alpha*x_a)
        diff_y=(m/y_a)-(2.0d0*alpha*y_a) !! TEST FOR CCG GRID
        diff_z=(n/z_a)-(2.0d0*alpha*z_a)




        
          grad_x=res*diff_x
          grad_y=res*diff_y
          grad_z=res*diff_z
         
         !! write(*,*)  grad_x,grad_y,grad_z
 




        !! write(*,*) res 
         end subroutine  position_wav



!!------------------------------------------MOMENTUM-WEB-------------------------------------------------------------------------------------


          
        subroutine momentum_wav (l,m,n,alpha,px,py,pz,xa,ya,za,res,res_diff_x,&
         &res_diff_y,res_diff_z)
        use factors
        implicit none
        real*8 , parameter :: zero=0.d0,&
                              pi= 3.14159265358979323846264338d0

 
        integer ::  l,&
                    m,&
                    n,&
                    i,&
                    j,&
                    k,&
                    ll,&
                    fact1,&
                    fact2,&
                    fact3,&
                    dfact1,&
                    dfact2,&
                    dfact3,&
                    fact_i,&
                    fact_l_i,&
                    diff_fact_l_i,&
                    fact_j,&
                    fact_m_j,&
                    diff_fact_m_j,&
                    fact_k,&
                    fact_n_k,&
                    diff_fact_n_k
        
        real*8  ::  alpha,&
                    px,&
                    py,&
                    pz,&
                    xa,&
                    ya,&
                    za,&
                    temp1,&
                    temp2,&
                    temp3,&
                    diff_temp1,&
                    diff_temp2,&
                    diff_temp3,&
                    r1,&
                    r2,&
                    r3,&
                    r4,&
                    psq,&
                    pda,&
                    np

        
      complex*16 :: t=(0.d0,1),&
                    temp,&
                    res,&
                    res_diff_x,&
                    res_diff_y,&
                    res_diff_z,&
                    c1,&
                    z,&
                    c2,&
                    tempwn 
  
     
 
          temp1=zero
          temp2=zero
          temp3=zero
          diff_temp1=zero
          diff_temp2=zero
          diff_temp3=zero
          res=zero
          res_diff_x=zero
          res_diff_y=zero
          res_diff_z=zero



          ll=l+m+n
          c1=(t**ll)  
          pda=(px*xa)+(py*ya)+(pz*za)
          psq=(px*px)+(py*py)+(pz*pz)
          c2=exp((-1)*(psq/(4*alpha))+(t*pda))
         
          call factorial(l,fact1)
          call factorial(m,fact2)
          call factorial(n,fact3)

          r1=fact1*fact2*fact3
          r2=(2*pi*alpha)**(0.75d0)

          call dfactorial(((2*l)-1),dfact1)
          call dfactorial(((2*m)-1),dfact2)
          call dfactorial(((2*n)-1),dfact3)
          
          r3=((dfact1*dfact2*dfact3)**0.5d0)
          temp=((c1*r1*c2))/(r2*r3)



          do i=0,int(l/2)
            call factorial(i,fact_i)
            call factorial((l-(2*i)),fact_l_i)
            temp1=temp1+(((-1)**i)*((px/sqrt(alpha))**(l-(2*i)))/((fact_i)*(fact_l_i)))
             !!write(*,*)(px/sqrt(alpha))**(l-(2*i)) 
          end do


            
          do i=0,int(l/2)
            call factorial(i,fact_i)
            call factorial((l-(2*i)),diff_fact_l_i)
            diff_temp1=diff_temp1+(((-1)**i)*&
            &(((1.d0/sqrt(alpha))**(l-(2*i)))*(l-(2*i))*(px**(l-(2*i)-1)))/((fact_i)*(diff_fact_l_i)))



            !! write(*,*) (px/sqrt(alpha))**(l-(2*i)-1)
          end do
             !! write(*,*) diff_temp1,temp1
              


         
          do j=0,int(m/2)
             call factorial(j,fact_j)
             call factorial((m-(2*j)),fact_m_j)
             temp2=temp2+(((-1)**j)*((py/sqrt(alpha))**(m-(2*j)))/((fact_j)*(fact_m_j)))      
          end do



         do j=0,int(m/2)
             call factorial(j,fact_j)
             call factorial((m-(2*j)),diff_fact_m_j)
             !diff_temp2=diff_temp2+(((-1)**j)*((py/sqrt(alpha))**(m-(2*j)-1))/((fact_j)*(diff_fact_m_j)))
             diff_temp2=diff_temp2+(((-1)**j)*&
             &(((1.d0/sqrt(alpha))**(m-(2*j)))*(m-(2*j))*(py**(m-(2*j)-1)))/((fact_j)*(diff_fact_m_j)))


              !! write(*,*) diff_fact_m_j       
          end do

               



          do k=0,int(n/2)
           call factorial(k,fact_k)
           call factorial((n-(2*k)),fact_n_k)
            temp3=temp3+(((-1)**k)*((pz/sqrt(alpha))**(n-(2*k)))/((fact_k)*(fact_n_k)))      
          end do
         ! write(*,*) px,py,pz


          do k=0,int(n/2)
           call factorial(k,fact_k)
           call factorial((n-(2*k)),diff_fact_n_k)
           !diff_temp3=diff_temp3+(((-1)**k)*((pz/sqrt(alpha))**(n-(2*k)-1))/((fact_k)*(diff_fact_n_k)))
            
           diff_temp3=diff_temp3+(((-1)**k)*&
           &(((1.d0/sqrt(alpha))**(n-(2*k)))*(n-(2*k))*(pz**(n-(2*k)-1)))/((fact_k)*(diff_fact_n_k)))
      
          end do
           !! write(*,*) diff_temp3

          
          res=temp*temp1*temp2*temp3


         
           res_diff_x= ((-px/(2*alpha))*res)+(t*xa*res)+(temp*diff_temp1*&
           &temp2*temp3)


           res_diff_y=((-py/(2*alpha))*res)+(t*ya*res)+(temp*temp1*&
           &diff_temp2*temp3)


          
           res_diff_z=((-pz/(2*alpha))*res)+(t*za*res)+(temp*temp1*&
           &temp2*diff_temp3)
 


  

           !! write(*,*) real(res_diff_x),real(res_diff_y),real(res_diff_z),real(res)














  400     format(8f20.16,2x,i10)   
          end subroutine momentum_wav

          end module basis_fun



          
























