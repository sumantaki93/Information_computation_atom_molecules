            program main_cal
            use  Property_cal
            use contraction_norm


           real*8, parameter   :: std_com_pa=1.d-4,&
                                  zero=0.d0


           real*8             ::   r_max,&
                                   r_min,&
                                   comp_grid,&
                                   res,&
                                   res1


        real*8, allocatable :: ci(:),&
                             z1(:),&
                             xxa(:),&
                             yya(:),&
                             zza(:),&
                             coef_a(:,:),&
                             coef_b(:,:),&
                             kj1(:),&
                             alphaj1(:)
 



         integer, allocatable :: npg1(:),&
                                 ll(:),&
                                 mm(:),&
                                 nn(:)





        integer :: noa,&
                   npg,&
                   noc,&
                   n_orb,&
                   rpflag,& 
                   j,&
                   k,&
                   temp,&
                   p


           integer :: npvec,&  
                      ntheta,&
                      nphi,&
                      i,&
                      method,&
                      n_a,&
                      n_b,&
                      end_q,&
                      num_threads,&
                      max_num_thread,&
                      r_min_index
               


            integer*4 ::  time(3),&    
                          date(3)

 
            integer :: hostnm, status 

             character*32 name
             !!status = hostnm( name )  


            real*8, allocatable   :: final_res_1(:),&
                                     final_res_2(:),&
                                     final_res_3(:),&
                                     final_res_4(:),&
                                     final_res_5(:),&
                                     final_res_6(:),&
                                     final_res_7(:),&
                                     final_res_8(:),&
                                     final_res_0(:),&
                                     final_res_10(:),&
                                     final_res_11(:),&
                                     final_res_15(:),&
                                     final_res_13(:),&
                                     final_res_16(:),&
                                     final_res_14(:),&
                                     final_res_17(:),&
                                     final_res_18(:),&
                                     final_res_19(:),&
                                     final_res_20(:),&
                                     c_i_s_position_n(:),&
                                     c_e_s_position_n(:),&
                                     c_i_s_position_1(:),&
                                     c_e_s_position_1(:),&
                                     c_i_s_momentum_n(:),&
                                     c_e_s_momentum_n(:),&
                                     c_i_s_momentum_1(:),&
                                     c_e_s_momentum_1(:)

                 open(25,file='noa_npg_noc_norb.dat',status='unknown')
                 read(25,*) noa,npg,noc,n_orb
                 close(25)


                 open(2,file='grid_point.dat',status='unknown')          
                 read(2,*) npvec,ntheta,nphi
                 close(2)

                 open(2,file='CPU.text',status='unknown')          
                 read(2,*)num_threads,max_num_thread 
                 close(2)
                 


                 

 

              end_q=2
   
                            allocate ( ci(noc),&
                                            z1(noa),&
                                           xxa(npg),&
                                           yya(npg),&
                                           zza(npg),&
                                    coef_a(noc,noc),&
                                    coef_b(noc,noc),&    
                                           kj1(npg),&
                                      alphaj1(npg),&
                                         npg1(noc),&
                                           ll(npg),& 
                                           mm(npg),&
                                           nn(npg),&
                               final_res_1(end_q ),&
                               final_res_2(end_q),&
                               final_res_3(end_q),&
                               final_res_4(end_q),&
                                 final_res_5(end_q),&
                               final_res_6(end_q),&
                              final_res_7(end_q),&
                              final_res_8(end_q),&
                              final_res_15(end_q),&
                              final_res_13(end_q),&
                              final_res_16(end_q),&
                              final_res_14(end_q),&
                              final_res_17(end_q),& 
                              final_res_18(end_q),&
                              final_res_19(end_q),&
                              final_res_20(end_q),&
                              final_res_10(end_q),&                
                              final_res_11(end_q),&
                              c_i_s_position_n(end_q),&
                              c_e_s_position_n(end_q),&
                              c_i_s_position_1(end_q),&
                              c_e_s_position_1(end_q),&
                              c_i_s_momentum_n(end_q),&
                              c_e_s_momentum_n(end_q),&
                              c_i_s_momentum_1(end_q),&
                              c_e_s_momentum_1(end_q),&
                              final_res_0(end_q) , stat=istat)

               if(istat /=0) stop "*******allocation failes*****"
  
              final_res_1(:)=zero
              final_res_2(:)=zero
              final_res_3(:)=zero
              final_res_4(:)=zero
              final_res_5(:)=zero
              final_res_6(:)=zero
              final_res_7(:)=zero
              final_res_8(:)=zero
              final_res_0(:)=zero 
              final_res_10(:)=zero
              final_res_11(:)=zero
              final_res_15(:)=zero
              final_res_13(:)=zero
              final_res_16(:)=zero
              final_res_14(:)=zero
              final_res_17(:)=zero 
              final_res_18(:)=zero
              final_res_19(:)=zero
              final_res_20(:)=zero
              c_i_s_position_n(:)=zero
              c_e_s_position_n(:)=zero
              c_i_s_position_1(:)=zero
              c_e_s_position_1(:)=zero
              c_i_s_momentum_n(:)=zero
              c_e_s_momentum_n(:)=zero
              c_i_s_momentum_1(:)=zero
              c_e_s_momentum_1(:)=zero




              ci(:)=zero
              z1(:)=zero
              xxa(:)=zero   
              yya(:)=zero
              zza(:)=zero
              coef_a(:,:)=zero
              coef_b(:,:)=zero
              kj1(:)=zero
              alphaj1(:)=zero
              npg1(:)=zero
              ll(:) =zero
              mm(:) =zero
              nn(:) =zero
              r_max=55.0d0     
              r_min=zero 
              res=zero  
              res1=zero
              temp=0


             status = hostnm( name ) 
  
            ! write(*,*) '___________________________________'
            ! write(*,*) 'JOB IS RUNNING ON HOST', ' ', name
            ! write(*,*) 'UNDER THE OPERATING SYSTEM LINUX'
            ! write(*,*) '____________________________________'


             open(3,file='atomic_number.dat',status='unknown')
             do i=1,noa
             read(3,*) z1(i) 
             end do
             close(3)

             open(10,file = 'npg1.dat',status='unknown')
             do i=1,noc
             read(10,*) npg1(i)
             end do
             close(10)
 
             open(9,file = 'geo.dat',status = 'unknown')
             do i=1,npg
             read(9,*) xxa(i),yya(i),zza(i)
             end do 
             close(9)


             open(11,file = 'lmn.dat',status= 'unknown') 
             do i=1,npg       
             read(11,*) ll(i),mm(i),nn(i)
             end do
             close(11)                     

      
             open(12,file = 'exponent.dat',status= 'unknown')
             do i=1,npg
             read(12,*) alphaj1(i)
             end do  
             close(12)

             open(13,file = 'precoef.dat',status= 'unknown')
             do i=1,npg
             read(13,*) kj1(i)
             end do
             close(13)
          
        
             open(14,file ='alpha_beta_method.dat',status ='unknown')
             read(14,*) n_a,n_b,method
             close(14)

             open(15,file = 'cmatup.dat',status = 'unknown')
             do i=1,noc
             read(15,*) coef_a(i, :)
             end do
             close(15) 
 
             open (16,file ='cmatdw.dat',status='unknown')
             do i=1,noc
             read(16,*) coef_b(i,:)
             end do 
             close(16)

      write(*,*)
      if( num_threads < 8) then
      write(*,*)"_____________________________________________"
      write(*,*)" MAXIMUM THREADS OF MACHINE:"
      !!write(*,*)"_____________________________________________"
      write(*,*) max_num_thread
      write(*,*) "NUMBER OF THREADS LOUNCH :"
      write(*,*) num_threads
      write(*,*)"______________________________________________"

      else if (num_threads == max_num_thread) then
      !!write(*,*)"____________________________________________________" 
     !! write(*,*) "WARNNING :: YOU ARE USING MAXIMUM NUMBER OF THREADS"
      !!write(*,*)"____________________________________________________ "
      end if
      write(*,*) 
      write(*,*)  
      if (num_threads > 1) then 
      !!write(*,*)"________PROGRAM RUNNING IN PARALLEL MODE :________"
      else
      write(*,*)"___________________________________________________"
      write(*,*) "PROGRAM RUNNING IN SERIAL MODE:"
      write(*,*)"___________________________________________________" 
      end if


             
             p=1
          do i=1,noc
             res1=zero
              do j=(1+temp),(npg1(i)+temp) 
                   
                do k=(1+temp),(npg1(i)+temp)
                 call func22(kj1(j),kj1(k),&
                 &xxa(j),yya(j),zza(j),alphaj1(j),ll(j),mm(j),nn(j),&
                 &xxa(k),yya(k),zza(k),alphaj1(k),ll(k),mm(k),nn(k),res)
                 res1=res1+res 
                 !!write(*,*) res1
                end do
              end do
                temp=temp+npg1(i)
                ci(i)=1.0d0/sqrt(res1)
               !! write(*,*) ci(i)
             end do
      


          
              !!SPACIAL CASE FOR H,He & H2(BEFORE THAT RANGE DENSITY WILL BE EXPLICITELY ZERO)
              if (n_a==1 .and. n_b==0) then 
              r_max=55.0d0
              end if

              if (n_a==1 .and. n_b==1) then 
              r_max=35.0d0
              end if
              !!__________________END______________________________________________             
 


            do i=1,2
                rpflag=i 
                if (rpflag==2) then


                 r_max=401.0d0
                
                if (n_a==1 .and. n_b==0) then 
                r_max=201.0d0
                end if

                if (n_a==1 .and. n_b==1) then 
                r_max=351.0d0
                end if
                !!__________________END______________________________________________   

                if (noa > 1) then
                r_max=901
                end if


                if (noa > 1 .and. n_a==1 .and. n_b==1) then
                r_max=201.0d0 
                end if

               

                r_min=zero

              end if
              call property(rpflag,method,n_a,n_b,max_num_thread,num_threads,r_max,r_min,npvec,ntheta,nphi,noa,npg,noc,n_orb,&
              &ci,z1,xxa,yya,zza,coef_a,coef_b,kj1,alphaj1,npg1,ll,mm,nn,final_res_1(i),final_res_2(i),final_res_3(i),&
              &final_res_4(i),final_res_5(i),final_res_6(i),final_res_7(i),final_res_8(i),final_res_0(i),final_res_10(i),&
              &final_res_11(i),final_res_15(i),final_res_13(i),final_res_16(i),final_res_14(i),final_res_17(i),&
              &final_res_18(i),final_res_19(i),final_res_20(i),c_i_s_position_n(i),c_i_s_position_1(i),c_e_s_position_n(i),&
              &c_e_s_position_1(i),c_i_s_momentum_n(i),c_i_s_momentum_1(i),c_e_s_momentum_n(i), c_e_s_momentum_1(i))   
                                                                       
                  
             if( rpflag==1) then

             open(20,file='information.dat',status='unknown')
             open(35,file='complexity.dat',status='unknown')
             write(20,11)
             write(20,*)
             write(20,41)
             write(20,400) final_res_0(i),final_res_10(i),final_res_15(i),final_res_13(i),final_res_17(i),final_res_18(i)
             write(35,13)
             write(35,*)
             write(35,47)
             write(35,450)c_i_s_position_n(i),c_i_s_position_1(i),c_e_s_position_n(i),c_e_s_position_1(i)   
             !close(20)
             !close(35)

             end if

              if ( rpflag==2) then
              
              !open(27,file='momentum_information.dat',status='unknown')
              !open(45,file='momentum_complexity.dat',status='unknown')
              write(20,*)
              write(20,*)
              write(20,12)
              write(20,*)
              write(20,41) 
              write(20,400) final_res_8(i),final_res_11(i),final_res_16(i),final_res_14(i),final_res_19(i),final_res_20(i)
              write(35,*)
              write(35,*)
              write(35,14)           
              write(35,*)
              write(35,47)
              write(35,450) c_i_s_momentum_n(i),c_i_s_momentum_1(i),c_e_s_momentum_n(i), c_e_s_momentum_1(i)
              close(35)
              close(20)

              open(50,file='Moments.dat',status='unknown')
              write(50,15)
              write(50,*)
              write(50,42)
              write(50,*)                 
              write(50,400)final_res_1(i),final_res_2(i),final_res_3(i),final_res_4(i),&
               &final_res_5(i),final_res_6(i),final_res_7(i)
              close(50)   
   


              end if 


 
                            

           end do  
 

             deallocate ( ci,&
                        z1,&
                        xxa,&
                        yya,&
                        zza,&
                        coef_a,&
                        coef_b,&    
                        kj1,&
                        alphaj1,&
                        npg1,&
                        ll,& 
                        mm,&
                        nn)
 
         call itime (time)
         call idate (date)

        !!write(*,*)
        !!write(*,*)
        !!write(*,50) "......................................."
        !!write(*,100) "PROGRAM RUN SUCCESSFULLY:"
        !!write(*,100) "PROGRAM DEVOLOPED BY......."
        !!write(*,100) "SUMAN HAZRA"
        !!write(*,200) date
        !!write(*,100) "PLACE: IISER-KOLKATA"
        !!write(*,300) time 
        !!write(*,50) "----------THANK YOU-----------"
  50    format (40x,25a) 
 200    format (50x,'DATE :',' ',i2.2, '/', i2.2, '/', i4.4)      
 100    format (50x,25a)
 300    format (50x, 'TIME :',' ' ,i2.2, ':', i2.2, ':', i2.2)
  55    format (60x,i50)
  30    format (2x,'TIME TAKEN FOR PARALLEL COMPUTING:'i10.2,1x,'MIN','', f10.2, ' SEC')
  35    format (2x,'TIME TAKEN FOR SERIAL COMPUTING:',1x,f10.2,' SEC')
  41    format (8x,'SHANNON(N)',10x,'SHANNON(N=1)',9x,'FISHER(N)',10x,&
         &'FISHER(N=1)',9x,'ONICESCU(N)',8X,'ONICESCU(N=1)')
  42    format (8x,'<P^-2>',14x,'<P^-1>',14x,'N/<P^0>',14x,&
         &'<P^1>',14x,'<P^2>/2*K.E',9X,'<P^3>',14x,'<P^4>')
 


  47    format (8x,'C(FS)/N',13x,'C(FS)/(N=1)',10x,'C(ES)/N',12x,&
         &'C(ES)/(N=1)')

 
  400    format(8e20.7,i10,i10,i10,i10)
  450    format(4e20.7) 
  106    format(f10.7,1x,f10.7)
   13    format(10x '___________ ATOMIC/MOLECULAR COMPLEXITY IN POSITION SPACE____________')
   14    format(10x '___________ ATOMIC/MOLECULAR COMPLEXITY IN MOMENTUM SPACE____________')
   15    format(42x '___________MOMENTUM MOMENTS FOR ATOMS OR MOLECULES____________')


   11    format (25x '___________INFORMATION THEORETIC MEASURE FOR POSITION SPACE____________')
   12    format (25x '___________INFORMATION THEORETIC MEASURE FOR MOMENTUM SPACE____________')


         end program main_cal
    



                                 
                                                                         
