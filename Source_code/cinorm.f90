             module contraction_norm
             contains

             subroutine func22(cn1,cn2,xa,ya,za,alpha1,l1,m1,n1,xb,yb,zb,&
             &alpha2,l2,m2,n2,res)
             use factors

             implicit none
             
             real*8,parameter :: zero=0.d0,&
                                 pi= 3.14159265358979323846264338d0


             integer :: l1,&
                        m1,&
                        n1,&
                        l2,&
                        m2,&
                        n2,&
                        i,&
                        ltl,&
                        ltm,&
                        ltn,&
                        dm1,&
                        dm2,&
                        dfactl1,&
                        dfactm1,&
                        dfactn1,&
                        dfactl2,&
                        dfactm2,&
                        dfactn2,&
                        dfacti

             real*8  :: cn1,&
                        cn2,&
                        xa,&
                        ya,&
                        za,&
                        alpha1,&
                        xb,&
                        yb,&
                        zb,&
                        alpha2,&
                        norm_1,&
                        norm_2,&
                        p1,&
                        p2,&
                        nm1,&
                        nm2,&           
                        f1i(0:11),&
                        f2i(0:11),&
                        f3i(0:11),&
                        sigma,&
                        ab,&
                        px,&
                        py,&
                        pz,&
                        apx,&
                        apy,&
                        apz,&
                        bpx,&
                        bpy,&
                        bpz,&
                        res,&
                        res1,&
                        res2,&
                        res3

               do i=0,10  
                   f1i(i)=zero
                   f2i(i)=zero
                   f3i(i)=zero
                end do

                res1=zero
                res2=zero
                res3=zero

                 

                p1=l1+m1+n1
                p2=l2+m2+n2
 
                nm1=(4.0d0*alpha1)**p1
                nm2=(4.0d0*alpha2)**p2
            
                call dfactorial(((2*l1)-1), dfactl1)
                call dfactorial(((2*m1)-1), dfactm1)
                call dfactorial(((2*n1)-1), dfactn1)

                call dfactorial(((2*l2)-1), dfactl2)
                call dfactorial(((2*m2)-1), dfactm2)
                call dfactorial(((2*n2)-1), dfactn2)


                dm1=dfactl1*dfactm1*dfactn1
                dm2=dfactl2*dfactm2*dfactn2

                norm_1=(((2.0d0*alpha1)/pi)**0.75d0)*((nm1/dm1)**0.50d0)
                norm_2=(((2.0d0*alpha2)/pi)**0.75d0)*((nm2/dm2)**0.50d0)

                ltl=(l1+l2)/2
                ltm=(m1+m2)/2
                ltn=(n1+n2)/2

                sigma=1.0d0/(alpha1+alpha2)

                ab=((xb-xa)**2.0d0)+((yb-ya)**2.0d0)+((zb-za)**2.0d0)
                px = sigma*((alpha1*xa)+(alpha2*xb)) 
                py = sigma*((alpha1*ya)+(alpha2*yb))
                pz = sigma*((alpha1*za)+(alpha2*zb))

                apx=px-xa
                bpx=px-xb
                apy=py-ya
                bpy=py-yb
                apz=pz-za
                bpz=pz-zb

                call func21(l1,l2,apx,bpx,f1i)

                 
                do i=0,ltl
                  call dfactorial(((2*i)-1), dfacti)
                  res1=res1+(f1i(2*i)*dfacti*((0.5d0*sigma)**i))
                end do



               call func21(m1,m2,apy,bpy,f2i)

                 do i=0,ltm
                  call dfactorial(((2*i)-1), dfacti)
                  res2=res2+(f2i(2*i)*dfacti*((0.5d0*sigma)**i))
                 end do


                call func21(n1,n2,apz,bpz,f3i)

                do i=0,ltn
                  call dfactorial(((2*i)-1), dfacti)
                  res3=res3+(f3i(2*i)*dfacti*((0.5d0*sigma)**i))
                end do


                res=cn1*cn2*norm_1*norm_2*(pi**1.50d0)*(sigma**1.50d0)*exp((-1)*alpha1*alpha2*sigma*ab)*&
                &res1*res2*res3

                end subroutine func22 


             
             subroutine func21(l,m,xa,ya,fp)
             use factors 
             implicit none
             real*8 , parameter :: zero=0.0d00                 
             integer :: l,&
                        m,&
                        i,&
                        r,&
                        k,&
                        l12,&
                        p,&
                        com1,&
                        com2
 
             real*8  :: xa,&
                        ya,&
                        fp(0:11),&
                        a1(0:11),&
                        b1(0:11),&
                        res(0:11,0:11),&
                        fi(0:11)
 
            p=0
            l12=l+m

            do i=0,10
               
             fi(i)=zero
            end do

             fi(:)=zero
             
            do r=0,l
               call combination (l,r,com1) 
               a1(r)=com1*(xa**(l-r))
            end do
               

             do k=0,m
               call combination (m,k,com2) 
               b1(k)=com1*(xa**(m-k))
            end do
 
             do r=0,l
             do k=0,m
               res(r,k)=a1(r)*b1(k)
             end do
             end do


            do r=0,l
            do k=0,m
              p=r+k
              fi(p)=fi(p)+res(r,k)
            end do
            end do


           do i=0,l12
               fp(i)=fi(i)
           end do 
           end subroutine func21 
 
 
           end module contraction_norm          



