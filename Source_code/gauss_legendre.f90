         module Property_cal 
         contains
         subroutine property(rpflag,method,n_a,n_b,max_num_thread,num_threads,r_max,r_min,npvec,ntheta,nphi,noa,npg,noc,n_orb,&
          &ci,z1,xxa,yya,zza,coef_a,coef_b,kj1,alphaj1,npg1,ll,mm,nn,final_res_1,final_res_2,final_res_3,&
         &final_res_4,final_res_5,final_res_6,final_res_7,final_res_8,final_res_0,final_res_10,final_res_11,&
          &final_res_15,final_res_13,final_res_16,final_res_14,final_res_17,final_res_18,final_res_19,final_res_20,&
          & c_i_s_position_n,c_i_s_position_1,c_e_s_position_n,c_e_s_position_1, c_i_s_momentum_n,c_i_s_momentum_1,&
          & c_e_s_momentum_n, c_e_s_momentum_1)
 
         use omp_lib
         use psix_fun
         real*8, parameter   ::  pi=4.d0*atan(1.d0),&
                                 twopi=2.d0*pi,&
                                 zero=0.d0,&
                                 one=1.d0


  
         integer,intent(in) :: noa,&
                               npg,&
                               noc,&
                               n_orb,&
                               max_num_thread,&
                               num_threads,&
                               npvec,&
                               ntheta,&
                               nphi,&
                               rpflag                 
                              


         real*8,intent (in) :: r_max,&
                               r_min    

      real*8, allocatable,intent(inout) :: ci(:),&
                                           z1(:),&
                                           xxa(:),&
                                           yya(:),&
                                           zza(:),&
                                           coef_a(:,:),&
                                           coef_b(:,:),&
                                           kj1(:),&
                                           alphaj1(:)


      integer, allocatable,intent(inout) :: npg1(:),&
                                            ll(:),&
                                            mm(:),&
                                            nn(:)

       integer, intent(inout) :: n_a,&
                                 n_b,&
                                 method

      real*8, allocatable :: rgd(:),&
                             tgd(:),&
                             pgd(:),&
                             x_grd(:),&
                             y_grd(:),&
                             z_grd(:),&
                             mwr(:),&
                             mwt(:),&
                             mwp(:),&
                             func(:),&
                             func2(:),&
                             m1(:),&
                             grad_rho_x(:),&
                             grad_rho_y(:),&
                             grad_rho_z(:),&
                             mod_grad_sq_x(:),&
                             mod_grad_sq_y(:),&
                             mod_grad_sq_z(:),&
                             m2(:),&
                             m3(:),&
                             diff_fact_x(:),&
                             diff_fact_y(:),&
                             diff_fact_z(:),&
                             px_diff(:),&
                             py_diff(:),&
                             pz_diff(:),&
                             p_sqr(:)
                              


       complex*16 ,allocatable :: density(:,:),&
                                  grad_density_x(:,:),&
                                  grad_density_y(:,:),&
                                  grad_density_z(:,:),&
                                  modsq_grad_density(:,:),&
                                  I_integrand(:,:)


       complex*16 :: res1,&
                     res1_diff_x,&
                     res1_diff_y,&
                     res1_diff_z,&
                     res2,&
                     res2_diff_x,&
                     res2_diff_y,&
                     res2_diff_z
                       
                             
       integer :: total_count,&
                  ip,&
                  ith,&
                  iph,&       
                  tnum,&  
                  ntnp,&
                  r_index,&
                  theta_index,&
                  phi_index,&
                  idx,&
                 ! n_a,&
                 ! n_b,&
                 ! method,&
                  ii,&
                  jj,&
                  kk
         
       real*8  :: theta_max,&
                  theta_min,&
                  phi_max,&
                  phi_min,&
                  px,&
                  py,&
                  pz,&
                  r_r,&
                  theta,&
                  phi,&
                  resth,&
                  resphi,&
                  sresphi,&
                  sresth,&
                  I_phi,&
                  I_theta,&
                  I_r,&
                  I_p,&
                  I_p_f,&
                  I_p_theta,&
                  I_p_phi,&
                  I_phi_2,&
                  I_theta_2,&
                  I_r_2,&
                  I_phi_3,&
                  I_theta_3,&
                  I_r_3,&
                  I_r_position,&
                  I_r_momentum,&
                  oni_r,&
                  oni_theta,&
                  oni_phi,&
                  oni_p,&
                  oni_theta_p,&
                  oni_phi_p,&
                  rm1,&
                  rm2,&
                  rm3,&
                  rm4,&
                  rm5,&
                  rm6,&
                  rm7,&
                  rm8,&
                  j_q,&
                  t1,&
                  t2,&
                  time_diff,&
                  time_diff_1,&
                  consnt_fact_1st_orde,&
                  h_theta,&
                  sum_diff_x,&
                  sum_diff_y,&
                  sum_diff_z,&
                  hx,&
                  h_phi,&
                  consnt_fact_1st_orde_theta,&
                  consnt_fact_1st_orde_phi,&
                  dx,&
                  dy,&
                  dz,&
                  xmn,&
                  ymn,&
                  zmn
                  



        real*8,intent(inout) :: final_res_1,&
                                final_res_2,&
                                final_res_3,&
                                final_res_4,&
                                final_res_5,&
                                final_res_6,&
                                final_res_7,&
                                final_res_8,&
                                final_res_10,& 
                                final_res_11,&
                                final_res_0,&
                                final_res_15,&
                                final_res_13,&
                                final_res_16,&
                                final_res_14,&
                                final_res_17,&
                                final_res_18,&
                                final_res_19,&
                                final_res_20,&
                                c_i_s_position_n,&
                                c_e_s_position_n,&
                                c_i_s_position_1,&
                                c_e_s_position_1,&
                                c_i_s_momentum_n,&
                                c_e_s_momentum_n,& 
                                c_i_s_momentum_1,&
                                c_e_s_momentum_1






                     

           integer*4 :: time(3),&    
                        date(3)


           real*8 :: final_res_12 


             total_count=npvec*ntheta*nphi
            
                 allocate ( rgd(npvec),&
                  x_grd(npvec),&
                  y_grd(ntheta),&
                  z_grd(nphi),&
                  mwr(npvec),&
                  m1(total_count),&
                  m2(total_count),&
                  m3(total_count),&
                  grad_rho_x(total_count),&
                  grad_rho_y(total_count),&
                  grad_rho_z(total_count),&
                  mod_grad_sq_x(total_count),&
                  mod_grad_sq_y(total_count),&
                  mod_grad_sq_z(total_count),&
                  tgd(ntheta),&                   
                  mwt(ntheta),&
                  pgd(nphi),&
                  density(0:1,(total_count)),&
                  grad_density_x(0:1,(total_count)),&
                  grad_density_y(0:1,(total_count)),&
                  grad_density_z(0:1,(total_count)),&
                  modsq_grad_density(0:1,(total_count)),&
                  func(total_count),&
                  px_diff(total_count),&
                  py_diff(total_count),&
                  pz_diff(total_count),&
                  diff_fact_x(total_count),&
                  diff_fact_y(total_count),&
                  diff_fact_z(total_count),&
                  p_sqr(total_count),&
                  I_integrand(0:1,(total_count)),&
                  func2(total_count),&
                  mwp(nphi),  stat=istat)
                  if(istat /=0) stop "******allocation error *****"

        
       
       theta_min=zero
       theta_max=pi

       phi_min=zero
        if( rpflag==1)then
       phi_max= 2.d0*pi
        end if

        if (rpflag==2) then 
        phi_max=pi 
        end if 


         

       rgd(:)=zero
       x_grd(:)=zero
       y_grd(:)=zero
       z_grd(:)=zero
       mwr(:)=zero
       tgd(:)=zero
       mwt(:)=zero
       pgd(:)=zero
       mwp(:)=zero
       func(:)=zero
       func2(:)=zero
       sum_diff_x=zero
       sum_diff_y=zero
       sum_diff_z=zero
       diff_fact_x(:)=zero
       diff_fact_y(:)=zero
       diff_fact_z(:)=zero
       grad_rho_x(:)=zero
       grad_rho_y(:)=zero
       grad_rho_z(:)=zero
       mod_grad_sq_x(:)=zero
       mod_grad_sq_y(:)=zero
       mod_grad_sq_z(:)=zero
       I_integrand(:,:)=zero
       p_sqr(:)=zero


       px=zero
       py=zero
       pz=zero
       r_r=zero 
       theta=zero
       phi=zero 
       ntnp=ntheta*nphi
       density(:,:)=zero
       grad_density_x(:,:)=zero
       grad_density_y(:,:)=zero
       grad_density_z(:,:)=zero
       modsq_grad_density(:,:)=zero

       I_integrand=zero
       I_r=zero
       I_p=zero
       I_p_f=zero
       I_p_theta=zero
       I_p_phi=zero
       rm1=zero
       rm2=zero
       rm3=zero
       rm4=zero
       rm5=zero
       rm6=zero
       rm7=zero
       rm8=zero
       j_q=zero
       m1=zero
 
       final_res_0=zero
       final_res_1=zero
       final_res_2=zero
       final_res_3=zero
       final_res_4=zero
       final_res_5=zero
       final_res_6=zero
       final_res_7=zero
       final_res_8=zero
       final_res_10=zero
       final_res_11=zero       
       final_res_12=zero 
       final_res_13=zero
       final_res_14=zero
       final_res_15=zero
       final_res_16=zero
       final_res_17=zero
       final_res_18=zero
       final_res_19=zero
       final_res_20=zero
       c_i_s_position_n=zero 
       c_i_s_position_1=zero 
       c_e_s_position_n=zero 
       c_e_s_position_1=zero 
       c_i_s_momentum_n=zero
       c_e_s_momentum_n=zero
       c_i_s_momentum_1=zero
       c_e_s_momentum_1=zero




         call gauleg(r_min,r_max,npvec,rgd,mwr)
         call gauleg(theta_min,theta_max,ntheta,tgd,mwt)      
         call gauleg(phi_min,phi_max,nphi,pgd,mwp)


      !$ call omp_set_num_threads (num_threads)      
         
       call cpu_time (t1)   


     !$omp parallel do private(i,ip,tnum,ith,iph,res1,res1_diff_x,res1_diff_y,res1_diff_z,res2, res2_diff_x,&                                                        
     !$omp res2_diff_y,res2_diff_z,px,py,pz,r_r,theta,phi) shared (noc,total_count,method,noa,npg,n_orb,npg1,ci,kj1,z1,ll,mm,nn,&
     !$omp alphaj1,xxa,yya,zza,coef_a,coef_b,n_a,n_b,rgd,tgd,pgd,density,grad_density_x,grad_density_y,grad_density_z,rpflag,&
     !$omp modsq_grad_density, I_integrand)  
          do i=1,total_count
          ip=(i/ntnp) +1
          if ( mod(i,ntnp)==0) then
          ip=ip-1
          end if
          tnum=(i/nphi)
          if ( mod(i,nphi)==0) then
          tnum=tnum-1
          end if
          ith=mod(tnum,ntheta)+1
          iph=mod(i,nphi)
          if ( mod(i,nphi).eq.0) then              
          iph=nphi
          end if
         !! write(*,*) ip,ith,iph,i
         !! write(*,*) ll(i)
 
        px = rgd(ip)*sin(tgd(ith))*cos(pgd(iph))
        py = rgd(ip)*sin(tgd(ith))*sin(pgd(iph))
        pz = rgd(ip)*cos(tgd(ith))








      call psix(method,rpflag,noc,noa,npg,n_orb,npg1,ci,kj1,z1,ll,mm,nn,alphaj1,xxa,yya,zza,coef_a,coef_b,&
           & px,py,pz,r_r,theta,phi,n_a,n_b,res1,res1_diff_x,res1_diff_y,res1_diff_y,res2,&
           &res2_diff_x,res2_diff_y,res2_diff_z)

           density(0,i)=res2
           density(1,i)=res1

           !if(rpflag==1) then
           grad_density_x(0,i)=res2_diff_x 
           grad_density_y(0,i)=res2_diff_y 
           grad_density_z(0,i)=res2_diff_z
 
           grad_density_x(1,i)=res1_diff_x
           grad_density_y(1,i)=res1_diff_y
           grad_density_z(1,i)=res1_diff_z




       modsq_grad_density(0,i)=abs(grad_density_x(0,i)*grad_density_x(0,i))+&
       &abs(grad_density_y(0,i)*grad_density_y(0,i))+&                             !! TEST FOR CCG GRID
       &abs(grad_density_z(0,i)*grad_density_z(0,i))
           



          
           I_integrand(0,i)=modsq_grad_density(0,i)/density(0,i)

           if (rpflag==2) then
            
            !!write(*,*) real(I_integrand(0,i)) 
           end if
           


      end do

       !$omp end parallel do
       call cpu_time (t2)
       time_diff = t2-t1
    
       !$ time_diff=time_diff/num_threads
         
          time_diff_1= (time_diff/60.d0)

      !!________________________________________________________________________
      !! ONE CAN SEE THE TIME EXECUTION FOR BOTH SPACE.
      !! write(*,*) 
      !! write(*,*)
      !! if ( num_threads > 1) then
      !! write(*,*) "------------------------------------"
      !! write(*,30) int(time_diff_1),mod(time_diff,60.d0)
      !! else 
       !!write(*,*) "---------------------------------------"
       !!write(*,35) time_diff
       !!end if 
      !!_________________________________________________________________________

              if(rpflag==1) then
              I_r=zero 
              oni_r=zero
              idx=1 
              do r_index=1,(npvec)
               resth=zero

              I_theta=zero
              sresth=zero
              oni_theta=zero

              
  
              do theta_index=1 ,ntheta   
              I_phi=zero
              sresphi=zero
              resphi=zero
              oni_phi=zero


              do phi_index=1,nphi
              I_phi=I_phi+(mwp(phi_index)*real(I_integrand(0,idx)))
              sresphi=sresphi-1.d0*(mwp(phi_index)*real(density(0,idx))*log(real(density(0,idx))))
              resphi=resphi+(mwp(phi_index)*real(density(0,idx)))
              oni_phi=oni_phi+(mwp(phi_index)*real(density(0,idx)**2))
              idx=idx+1
              end do

              I_theta=I_theta + (mwt(theta_index)*sin(tgd(theta_index))*I_phi)
             
              sresth=sresth + (mwt(theta_index)*sin(tgd(theta_index))*sresphi)

              resth=resth + (mwt(theta_index)*sin(tgd(theta_index))*resphi)
               
              oni_theta=oni_theta+(mwt(theta_index)*sin(tgd(theta_index))*oni_phi)

              end do
            
              I_r= I_r+(mwr(r_index)*I_theta*(rgd(r_index)**2.d0))
              rm8=rm8+(mwr(r_index)*(rgd(r_index)**2.d0)*sresth)
              rm3=rm3+(mwr(r_index)*(rgd(r_index)**2.d0)*resth)
              oni_r=oni_r+(mwr(r_index)*oni_theta*(rgd(r_index)**2.d0))

             
              end do
                final_res_15=I_r
                final_res_0 = rm8
                final_res_12 =rm3
                final_res_10 = (final_res_0/final_res_12) + log(final_res_12)
                final_res_13= final_res_15/ final_res_12
                final_res_17= oni_r
                final_res_18= oni_r/(final_res_12**2)

                c_i_s_position_n= final_res_15*exp((2.0d0* final_res_0)/3.0d0 )
                c_i_s_position_1= final_res_13*exp((2.0d0* final_res_10)/3.0d0 )
                c_e_s_position_n= final_res_17*exp(final_res_0)
                c_e_s_position_1= final_res_18*exp(final_res_10)


                end if



       !! SPHERICAL POLAR GRID 



            if (rpflag==2) then
   
             oni_p=zero
             idx=1 
           do r_index=1,npvec
           resth=zero
           sresth=zero
           
           I_p_theta=zero
           oni_theta_p=zero
                  
 




           do theta_index=1 ,ntheta   
           resphi=zero
           sresphi=zero
           I_p_phi=zero
           oni_phi_p=zero


           do phi_index=1,nphi
           resphi=resphi+(mwp(phi_index)*real(density(0,idx)))
           I_p_phi=I_p_phi+(mwp(phi_index)*real(I_integrand(0,idx)))
           sresphi=sresphi-1.d0*(mwp(phi_index)*real(density(0,idx))*log(real(density(0,idx))))
           oni_phi_p=oni_phi_p+(mwp(phi_index)*real(density(0,idx)**2))
           idx=idx+1
           end do
 
         resth=resth + (mwt(theta_index)*sin(tgd(theta_index))*resphi)
         I_p_theta=I_p_theta+(mwt(theta_index)*sin(tgd(theta_index))*I_p_phi)
         sresth=sresth + (mwt(theta_index)*sin(tgd(theta_index))*sresphi)
         oni_theta_p=oni_theta_p+(mwt(theta_index)*sin(tgd(theta_index))*oni_phi_p)
         end do

            rm1=rm1+(mwr(r_index)*resth)
            rm2=rm2+(mwr(r_index)*(rgd(r_index)**1.d0)*resth)
            rm3=rm3+(mwr(r_index)*(rgd(r_index)**2.d0)*resth)
            rm4=rm4+(mwr(r_index)*(rgd(r_index)**3.d0)*resth)
            rm5=rm5+(mwr(r_index)*(rgd(r_index)**4.d0)*resth)
            rm6=rm6+(mwr(r_index)*(rgd(r_index)**5.d0)*resth)
            rm7=rm7+(mwr(r_index)*(rgd(r_index)**6.d0)*resth)
            I_p= I_p+(mwr(r_index)*I_p_theta*(rgd(r_index)**2.d0)) 
            rm8=rm8+(mwr(r_index)*(rgd(r_index)**2.d0)*sresth) 
            oni_p=oni_p+(mwr(r_index)*oni_theta_p*(rgd(r_index)**2.d0))
          end do
          
          end if       
  

          deallocate ( rgd,&
                       mwr,&
                       tgd,&                   
                       mwt,&
                       pgd,&
                   density,&
                        mwp )

         
        if (rpflag==2) then
        final_res_1 =(2.0d000)*rm1 
        final_res_2 =(2.0d000)*rm2 
        final_res_3 =(2.0d000)*rm3 
        final_res_4 =(2.0d000)*rm4 
        final_res_5 =(2.0d000)*rm5 
        final_res_6 =(2.0d000)*rm6 
        final_res_7 =(2.0d000)*rm7 
        final_res_8 =(2.0d000)*rm8
        final_res_11=(final_res_8/final_res_3) + log(final_res_3)
        final_res_16=(2.0d00*I_p)        
        final_res_14=final_res_16/final_res_3
        final_res_19=(2.0d00*oni_p)
        final_res_20= final_res_19/(final_res_3**2)
        c_i_s_momentum_n= final_res_16*exp((2.0d0* final_res_8)/3.0d0 )
        c_i_s_momentum_1= final_res_14*exp((2.0d0* final_res_11)/3.0d0 )
        c_e_s_momentum_n= final_res_19*exp(final_res_8) 
        c_e_s_momentum_1= final_res_20*exp(final_res_11)

 

        end if 

 
          call itime (time)
          call idate (date)

  30  format(2x,'EXECUTION TIME FOR MUTICORE PROCESSING:'i10.2,1x,'MIN','', f10.2, ' SEC')
  35  format(2x,'EXECUTION TIME FOR SERIAL PROCESSING:',1x,f10.2,' SEC')
  40   format (f15.16)
       end subroutine Property



      subroutine gauleg(x1,x2,n,x,w)
      real*8 ,parameter ::EPS=3.d-14
      integer :: n,&
                 i,&
                 j,& 
                 m

      real*8  :: x1,&
                 x2,&
                 x(n),&
                 w(n),&
                 p1,&
                 p2,&
                 p3,&
                 pp,&
                 xm,&
                 z,&
                 z1

        m=(n+1)/2
        xm=0.5d0*(x2+x1)
        x1=0.5d0*(x2-x1)
        
        do i=1,m
             z=cos(3.141592654d0*(i-0.25d0)/(n+0.5d0))
     1  continue

        p1=1.d0
        p2=0.d0   
        
        do j=1,n
             p3=p2
             p2=p1
             p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
        end do

        pp=n*(z*p1-p2)/(z*z-1.d0)
        z1=z
        z=z1-p1/pp
       
        if(abs(z-z1).gt.EPS) goto 1 
        x(i)=xm-x1*z
        x(n+1-i)=xm+x1*z
        w(i)=2.d0*x1/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
        end do
        return
        end subroutine gauleg



        end module Property_cal
