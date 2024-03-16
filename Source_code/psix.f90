        module psix_fun
        contains
         subroutine psix(method,rpflag,noc,noa,npg,n_orb,npg1,ci,kj1,z1,ll,mm,nn,alphaj1,xxa,yya,zza,coef_a,coef_b,&
               & xx,yy,zz,r_r,theta,phi,n_a,n_b,r1,r1_diff_x,r1_diff_y,r1_diff_z,r2,r2_diff_x,r2_diff_y,r2_diff_z)     
        use basis_fun
       
       real*8 , parameter :: zero=0.0000d0
       integer,intent(in) :: noa,&
                             noc,&
                             npg,&
                             n_orb,&
                             n_b,&
                             rpflag,&
                             method

       integer :: ki,i,j,temp,now,big
         
       integer , intent(inout) :: n_a


      integer, allocatable,intent(inout) :: npg1(:),&
                                         ll(:),&
                                         mm(:),&
                                         nn(:)


 


      real*8 ,intent(inout) ::    xx,&
                                  yy,&
                                  zz,&
                                  r_r,&
                                  theta,&
                                  phi


      real*8, allocatable,intent (inout) :: ci(:),&
                                         z1(:),&
                                         xxa(:),&
                                         yya(:),&
                                         zza(:),&
                                         coef_a(:,:),&
                                         kj1(:),&
                                         alphaj1(:),&
                                         coef_b(:,:)


      complex*16   :: func,&
                      res1,&
                      ires_a,&
                      ires_a_diff_x,&
                      ires_a_diff_y,&
                      ires_a_diff_z,&
                      mo_ires_a_diff_x,&   
                      mo_ires_a_diff_y,& 
                      mo_ires_a_diff_z,&
                      iires_a,&
                      iires_a_diff_x,&
                      iires_a_diff_y,&
                      iires_a_diff_z,&
                      mo_iires_a_diff_x,&  
                      mo_iires_a_diff_y,&
                      mo_iires_a_diff_z,&
                     
                      mo_iires_a_diff_x2,& 
                      mo_iires_a_diff_y2,&
                      mo_iires_a_diff_z2,&


                      ires_b,&
                      ires_b_diff_x,&
                      ires_b_diff_y,&
                      ires_b_diff_z,&
              
                      mo_ires_b_diff_x,&
                      mo_ires_b_diff_y,&
                      mo_ires_b_diff_z,&
 
                      


                      iires_b,&
                      iires_b_diff_x,&
                      iires_b_diff_y,&
                      iires_b_diff_z,&
                        
                      mo_iires_b_diff_x,&
                      mo_iires_b_diff_y,&
                      mo_iires_b_diff_z,&

       
                      mo_iires_b_diff_x2,&  
                      mo_iires_b_diff_y2,&
                      mo_iires_b_diff_z2,&


                      func_diff_x,&
                      func_diff_y,&
                      func_diff_z,&
                      mo_func_diff_x,&
                      mo_func_diff_y,&
                      mo_func_diff_z,&
                      res1_diff_x,&
                      res1_diff_y,&
                      res1_diff_z,&
                      mo_res1_diff_x,&
                      mo_res1_diff_y,&
                      mo_res1_diff_z

      complex*16,intent(inout) :: r1,&
                                  r2,&
                                  r1_diff_x,&
                                  r1_diff_y,& 
                                  r1_diff_z,&
                                  r2_diff_x,&
                                  r2_diff_y,&
                                  r2_diff_z
 


      complex*16, allocatable  :: psi(:),&
                                  res(:),&
                                  res_diff_x(:),&
                                  res_diff_y(:),&
                                  res_diff_z(:),&
                                  mo_res_diff_x(:),&
                                  mo_res_diff_y(:),&
                                  mo_res_diff_z(:),&
                                  psi_diff_x(:),&
                                  psi_diff_y(:),&
                                  psi_diff_z(:),&
                                  mo_psi_diff_x(:),&
                                  mo_psi_diff_y(:),&
                                  mo_psi_diff_z(:),&
                                  mo_grad_x(:),&
                                  mo_grad_y(:),&
                                  mo_grad_z(:)



       
       real*8 ,allocatable :: res_po (:),&
                               grad_x(:),&
                               grad_y(:),&
                               grad_z(:)                           

                       allocate (res(npg),&
                                 res_po(npg),&
                                 grad_x(npg),&
                                 grad_y(npg),&
                                 res_diff_x(npg),&
                                 res_diff_y(npg),&
                                 res_diff_z(npg),&
                                 mo_res_diff_x(npg),&
                                 mo_res_diff_y(npg),&
                                 mo_res_diff_z(npg),&
                                 grad_z(npg),&
                                 mo_grad_x(npg),&
                                 mo_grad_y(npg),&
                                 mo_grad_z(npg),&
                                 psi_diff_x(noc),&
                                 psi_diff_y(noc),&
                                 mo_psi_diff_x(noc),&
                                 mo_psi_diff_y(noc),&
                                 mo_psi_diff_z(noc),&
                                 psi_diff_z(noc),&
                                 psi(noc) , stat=istat)
                      if(istat /=0) stop "*******allocation fails******"
            psi(:)=zero
            temp=0
            now=0
            func=zero
            func_diff_x=zero  
            func_diff_y=zero
            func_diff_z=zero
            mo_func_diff_x=zero
            mo_func_diff_y=zero
            mo_func_diff_z=zero
            

            ires_a=zero
            ires_a_diff_x=zero
            ires_a_diff_y=zero
            ires_a_diff_z=zero
            
            mo_ires_a_diff_x=zero
            mo_ires_a_diff_y=zero
            mo_ires_a_diff_z=zero


 

            iires_a_diff_x=zero
            iires_a_diff_y=zero
            iires_a_diff_z=zero
            mo_iires_a_diff_x=zero 
            mo_iires_a_diff_y=zero
            mo_iires_a_diff_z=zero
      
            mo_iires_a_diff_x2=zero      
            mo_iires_a_diff_y2=zero
            mo_iires_a_diff_z2=zero

            ires_b_diff_x=zero
            ires_b_diff_y=zero
            ires_b_diff_z=zero
 
            mo_ires_b_diff_x=zero 
            mo_ires_b_diff_y=zero
            mo_ires_b_diff_z=zero

            iires_b_diff_x=zero
            iires_b_diff_y=zero
            iires_b_diff_z=zero
            mo_iires_b_diff_x=zero
            mo_iires_b_diff_y=zero
            mo_iires_b_diff_z=zero

            mo_iires_b_diff_x2=zero
            mo_iires_b_diff_y2=zero
            mo_iires_b_diff_z2=zero

            ires_b=zero
            iires_a=zero
            iires_b=zero
            res(:)=zero
            res_po(:)=zero
            grad_x(:)=zero      
            grad_y(:)=zero   
            grad_z(:)=zero
            res1_diff_x=zero
            res1_diff_y=zero
            res1_diff_z=zero
            mo_res1_diff_x=zero
            mo_res1_diff_y=zero
            mo_res1_diff_z=zero
            psi_diff_x(:)=zero
            psi_diff_y(:)=zero
            psi_diff_z(:)=zero
            mo_psi_diff_x(:)=zero
            mo_psi_diff_y(:)=zero
            mo_psi_diff_z(:)=zero
            res_diff_x(:)=zero
            res_diff_y(:)=zero
            res_diff_z(:)=zero
            mo_res_diff_x(:)=zero
            mo_res_diff_y(:)=zero
            mo_res_diff_z(:)=zero
            mo_grad_x(:)=zero
            mo_grad_y(:)=zero
            mo_grad_z(:)=zero

           r1_diff_x=zero
           r1_diff_y=zero
           r1_diff_z=zero

           r2_diff_x=zero
           r2_diff_y=zero
           r2_diff_z=zero



          do i=1,noa
             now=now+z1(i)
          end do
      
          if (method ==1) then
          n_a=now/2
          end if

          do i=1,noc
            res1=zero
            res(:)=zero

            res1_diff_x=zero
            res1_diff_y=zero
            res1_diff_z=zero

            res_diff_x(:)=zero
            res_diff_y(:)=zero
            res_diff_z(:)=zero




            mo_res1_diff_x=zero
            mo_res1_diff_y=zero
            mo_res1_diff_z=zero
          
           mo_res_diff_x(:)=zero
           mo_res_diff_y(:)=zero
           mo_res_diff_z(:)=zero


          do j=1+temp,npg1(i)+temp
          if ( rpflag==1) then
           call position_wav(ll(j),mm(j),nn(j),alphaj1(j),xx,yy,zz,r_r,theta,phi,xxa(j),yya(j),&
          &zza(j),res_po(j),grad_x(j),grad_y(j),grad_z(j))
           res(j)=res_po(j)
           res_diff_x(j)=grad_x(j)
           res_diff_y(j)=grad_y(j)
           res_diff_z(j)=grad_z(j)
           !! write(*,*) real(grad_x(j)),real(grad_y(j)),real(grad_z(j)),res_po(j) 

           end if

           
           if (rpflag==2) then
          call momentum_wav (ll(j),mm(j),nn(j),alphaj1(j),xx,yy,zz,xxa(j),yya(j),&
          &zza(j),res(j),mo_grad_x(j),mo_grad_y(j),mo_grad_z(j))

           mo_res_diff_x(j)=mo_grad_x(j)
           mo_res_diff_y(j)=mo_grad_y(j)
           mo_res_diff_z(j)=mo_grad_z(j)

           ! write(*,*) mo_res_diff_x(j),mo_res_diff_y(j),mo_res_diff_z(j)

           end if


          func=res(j)
          res1 = res1+kj1(j)*func

          if(rpflag==1) then 
          func_diff_x=res_diff_x(j)
          func_diff_y=res_diff_y(j)
          func_diff_z=res_diff_z(j)
           
         !! write(*,*) real(func_diff_x),real(func)

          res1_diff_x= res1_diff_x + kj1(j)*func_diff_x
          res1_diff_y= res1_diff_y + kj1(j)*func_diff_y
          res1_diff_z= res1_diff_z + kj1(j)*func_diff_z

          !! write(*,*) real(res1_diff_x),real(res1),kj1(j)
          end if

          if (rpflag==2) then
            mo_func_diff_x=mo_res_diff_x(j)
            mo_func_diff_y=mo_res_diff_y(j)
            mo_func_diff_z=mo_res_diff_z(j)

            mo_res1_diff_x=mo_res1_diff_x + kj1(j)*mo_func_diff_x
            mo_res1_diff_y=mo_res1_diff_y + kj1(j)*mo_func_diff_y
            mo_res1_diff_z=mo_res1_diff_z + kj1(j)*mo_func_diff_z

           !! write(*,*) real( mo_res1_diff_x),real(mo_res1_diff_y),real(mo_res1_diff_z)
          end if           
 




          end do
          temp=temp+npg1(i)
          psi(i)=res1

          if(rpflag==1) then 
          psi_diff_x(i)=res1_diff_x
          psi_diff_y(i)=res1_diff_y
          psi_diff_z(i)=res1_diff_z
          
          end if

          if(rpflag==2) then
             mo_psi_diff_x(i)=mo_res1_diff_x        
             mo_psi_diff_y(i)=mo_res1_diff_y
             mo_psi_diff_z(i)=mo_res1_diff_z
            !! write(*,*) real(mo_psi_diff_x(i)),real(mo_psi_diff_y(i)),real(mo_psi_diff_z(i))
             
          end if








          end do 

         if (n_orb > n_a) then
             big=n_orb
          else
            big=n_a
        end if 


        do ki=1,big
            ires_a = zero
            if(rpflag==1) then
              ires_a_diff_x=zero
              ires_a_diff_y=zero
              ires_a_diff_z=zero
             end if                                 
 
            if (rpflag==2) then
               mo_ires_a_diff_x=zero
               mo_ires_a_diff_y=zero
               mo_ires_a_diff_z=zero
             end if

 
 
           do i=1,noc
              ires_a=ires_a+(coef_a(ki,i)*ci(i)*psi(i))
         
             if (rpflag==1) then
             ires_a_diff_x=ires_a_diff_x+(coef_a(ki,i)*ci(i)*psi_diff_x(i))
             ires_a_diff_y=ires_a_diff_y+(coef_a(ki,i)*ci(i)*psi_diff_y(i))
             ires_a_diff_z=ires_a_diff_z+(coef_a(ki,i)*ci(i)*psi_diff_z(i))
             !write(*,*) ires_a_diff_x
             end if

             if (rpflag==2) then
             mo_ires_a_diff_x=mo_ires_a_diff_x+(coef_a(ki,i)*ci(i)*mo_psi_diff_x(i))  
             mo_ires_a_diff_y=mo_ires_a_diff_y+(coef_a(ki,i)*ci(i)*mo_psi_diff_y(i))  
             mo_ires_a_diff_z=mo_ires_a_diff_z+(coef_a(ki,i)*ci(i)*mo_psi_diff_z(i))
           
            !! write(*,*) real(mo_ires_a_diff_x),real(mo_ires_a_diff_y),real(mo_ires_a_diff_z)
             end if           
 



           end do
           if( ki==n_orb) then
            r1=ires_a 
           if(rpflag==1) then
            r1_diff_x=ires_a_diff_x
            r1_diff_y=ires_a_diff_y 
            r1_diff_z=ires_a_diff_z
            end if
                                                        !! NOT APPLICABLE
           if(rpflag==2) then
            r1_diff_x=mo_ires_a_diff_x 
            r1_diff_y=mo_ires_a_diff_y
            r1_diff_z=mo_ires_a_diff_z
           end if

           end if
           
              !! write(*,*) real(ires_a_diff_x),real(ires_a)
           

           if ( ki <= n_a) then
              if (method ==1) then 
             iires_a = iires_a + 2.d0*(conjg(ires_a)*ires_a)
                 
             if (rpflag==1) then
             iires_a_diff_x=iires_a_diff_x +2.d0*(conjg(ires_a)*ires_a_diff_x)
             iires_a_diff_y=iires_a_diff_y +2.d0*(conjg(ires_a)*ires_a_diff_y)
             iires_a_diff_z=iires_a_diff_z +2.d0*(conjg(ires_a)*ires_a_diff_z)
             end if
               !! write(*,*) real(iires_a_diff_x),real( iires_a )

             if (rpflag==2) then
              mo_iires_a_diff_x=mo_iires_a_diff_x + 2.d0*(conjg(ires_a)*mo_ires_a_diff_x)
              mo_iires_a_diff_y=mo_iires_a_diff_y + 2.d0*(conjg(ires_a)*mo_ires_a_diff_y)
              mo_iires_a_diff_z=mo_iires_a_diff_z + 2.d0*(conjg(ires_a)*mo_ires_a_diff_z)

              mo_iires_a_diff_x2=mo_iires_a_diff_x2 + 2.d0*(ires_a*conjg(mo_ires_a_diff_x))
              mo_iires_a_diff_y2=mo_iires_a_diff_y2 + 2.d0*(ires_a*conjg(mo_ires_a_diff_y))
              mo_iires_a_diff_z2=mo_iires_a_diff_z2 + 2.d0*(ires_a*conjg(mo_ires_a_diff_z))



             end if
           


             
      !!!----------------------------------------------------------------------------------------------------FOR UNRESTRICTED CALCULATION------------------------------------------- 




              else if ( method == 2) then
             iires_a = iires_a + (conjg(ires_a)*ires_a)

              if(rpflag==1) then
              iires_a_diff_x=iires_a_diff_x +(conjg(ires_a)*ires_a_diff_x)                  
              iires_a_diff_y=iires_a_diff_y +(conjg(ires_a)*ires_a_diff_y)
              iires_a_diff_z=iires_a_diff_z +(conjg(ires_a)*ires_a_diff_z)
              end if

              if(rpflag==2) then
             
              mo_iires_a_diff_x=mo_iires_a_diff_x + (conjg(ires_a)*mo_ires_a_diff_x)
              mo_iires_a_diff_y=mo_iires_a_diff_y + (conjg(ires_a)*mo_ires_a_diff_y)
              mo_iires_a_diff_z=mo_iires_a_diff_z + (conjg(ires_a)*mo_ires_a_diff_z)
                                                                                          
              mo_iires_a_diff_x2=mo_iires_a_diff_x2 + (ires_a*conjg(mo_ires_a_diff_x))
              mo_iires_a_diff_y2=mo_iires_a_diff_y2 + (ires_a*conjg(mo_ires_a_diff_y))
              mo_iires_a_diff_z2=mo_iires_a_diff_z2 + (ires_a*conjg(mo_ires_a_diff_z))
           
              end if

                   


   
              end if
           end if 







        end do
               if(rpflag==1) then
              ! write(*,*) real(iires_a_diff_x),real(iires_a_diff_y),real(iires_a_diff_z),real(iires_a)
              end if

 


       !!!----------------------------------------------------------------------------------------------------FOR UNRESTRICTED CALCULATION------------------------------------------- 


           if ( method == 2) then
         
           do ki=1,n_b
            ires_b = zero

            if(rpflag==1) then
             ires_b_diff_x=zero
             ires_b_diff_y=zero
             ires_b_diff_z=zero
             end if

           if (rpflag==2) then 
             mo_ires_b_diff_x=zero 
             mo_ires_b_diff_y=zero
             mo_ires_b_diff_z=zero

           end if



 
           do i=1,noc
              ires_b=ires_b+(coef_b(ki,i)*ci(i)*psi(i))
             
            if (rpflag==1) then
             ires_b_diff_x=ires_b_diff_x+(coef_b(ki,i)*ci(i)*psi_diff_x(i))
             ires_b_diff_y=ires_b_diff_y+(coef_b(ki,i)*ci(i)*psi_diff_y(i))
             ires_b_diff_z=ires_b_diff_z+(coef_b(ki,i)*ci(i)*psi_diff_z(i))
            !!  write(*,*) ires_b_diff_x,ires_b
             end if


            if (rpflag==2) then
             mo_ires_b_diff_x=mo_ires_b_diff_x+(coef_b(ki,i)*ci(i)*mo_psi_diff_x(i))
             mo_ires_b_diff_y=mo_ires_b_diff_y+(coef_b(ki,i)*ci(i)*mo_psi_diff_y(i))
             mo_ires_b_diff_z=mo_ires_b_diff_z+(coef_b(ki,i)*ci(i)*mo_psi_diff_z(i))
            end if

           end do



           if( ki==n_orb) then
            r1=ires_b

            if(rpflag==1) then
            r1_diff_x=ires_b_diff_x
            r1_diff_y=ires_b_diff_y 
            r1_diff_z=ires_b_diff_z
            end if                                       !! NOT APPLICABLE
            
            if(rpflag==2) then
            r1_diff_x=mo_ires_b_diff_x 
            r1_diff_y=mo_ires_b_diff_y 
            r1_diff_z=mo_ires_b_diff_z 
            end if 
 
           end if
 
              !! write(*,*) ires_b,ires_b_diff_z

           if ( ki <= n_b) then
               iires_b = iires_b + (conjg(ires_b)*ires_b)
           
              if(rpflag==1) then
              iires_b_diff_x=iires_b_diff_x +(conjg(ires_b)*ires_b_diff_x)                  
              iires_b_diff_y=iires_b_diff_y +(conjg(ires_b)*ires_b_diff_y)               
              iires_b_diff_z=iires_b_diff_z +(conjg(ires_b)*ires_b_diff_z)
              end if

              if (rpflag==2) then
              mo_iires_b_diff_x=mo_iires_b_diff_x +(conjg(ires_b)*mo_ires_b_diff_x) 
              mo_iires_b_diff_y=mo_iires_b_diff_y +(conjg(ires_b)*mo_ires_b_diff_y)
              mo_iires_b_diff_z=mo_iires_b_diff_z +(conjg(ires_b)*mo_ires_b_diff_z)

              mo_iires_b_diff_x2=mo_iires_b_diff_x2 +(ires_b*conjg(mo_ires_b_diff_x))
              mo_iires_b_diff_y2=mo_iires_b_diff_y2 +(ires_b*conjg(mo_ires_b_diff_y))
              mo_iires_b_diff_z2=mo_iires_b_diff_z2 +(ires_b*conjg(mo_ires_b_diff_z))

              end if



           end if

           end do
           end if 
          


 



          
           deallocate ( res,&
                        psi)      

        r2 = iires_a+iires_b                      !! FINAL DENSITY 

        if(rpflag==1) then
           r2_diff_x=2.d0*(iires_a_diff_x+iires_b_diff_x)
           r2_diff_y=2.d0*(iires_a_diff_y+iires_b_diff_y)   !! FINAL GRADIENT OF DENSITY IN POSITION SPACE
           r2_diff_z=2.d0*(iires_a_diff_z+iires_b_diff_z)

         !! write(*,*) real(r2_diff_x)


         end if
 
        if(rpflag==2)then
          r2_diff_x=mo_iires_a_diff_x+ mo_iires_a_diff_x2+ mo_iires_b_diff_x+ mo_iires_b_diff_x2 
          r2_diff_y=mo_iires_a_diff_y+ mo_iires_a_diff_y2+ mo_iires_b_diff_y+ mo_iires_b_diff_y2  
          r2_diff_z=mo_iires_a_diff_z+ mo_iires_a_diff_z2+ mo_iires_b_diff_z+ mo_iires_b_diff_z2

         !! write(*,*) real(r2_diff_x),real(r2_diff_y),real(r2_diff_z)   !! FINAL GRADIENT OF DENSITY IN MOMENTUM SPACE
        end if 







       
       end subroutine psix
       end module  psix_fun
                





                 

  
