#SUBROUTINE compac

  ## It calculates sediment compaction
  #USE variables
  #IMPLICIT NONE

## local variables
  #REAL(KIND=8) a(6),b(6,5),ch(6),ct(6),f(6)
  #REAL(KIND=8) column(jti),oldColumn(jti),sstor(nTotalMat)
  #REAL(KIND=8) mate(nTotalMat),dPor(nTotalMat)
  #REAL(KIND=8) sume,avedummy,consum,z,avepo,dp,avepor
  #REAL(KIND=8) x,h,hold,te,temax,ytemp,xtemp,dpdt
## counter
  #INTEGER j,k,mat,node,k1,l


      #a(2)=0.222222222222222      #2.0/9.0
      #a(3)=0.333333333333333      #1.0/3.0
      #a(4)=0.75                   #3.0/4.0
      #a(5)=1.0
      #a(6)=0.833333333333334      #5.0/6.0
      #b(2,1)=0.222222222222222    #2.0/9.0
      #b(3,1)=0.083333333333333    #1.0/12.0
      #b(3,2)=0.25                 #1.0/4.0
      #b(4,1)=0.5390625            #69.0/128.0
      #b(4,2)=-1.8984375           #-243.0/128.0
      #b(4,3)=2.109375             #135.0/64.0
      #b(5,1)=-1.4166666666666667  #-17.0/12.0
      #b(5,2)=6.75                 #27.0/4.0
      #b(5,3)=-5.4                 #-27.0/5.0
      #b(5,4)=1.0666666666666667   #16.0/15.0
      #b(6,1)=0.15046296296296297  #65.0/432.0
      #b(6,2)=-0.3125              #-5.0/16.0
      #b(6,3)=0.8125               #13.0/16.0
      #b(6,4)=0.14814814814814814  #4.0/27.0
      #b(6,5)=0.034722222222222224 #5.0/144.0
      #ch(1)=0.10444444444444445   #47.0/450.0
      #ch(2)=0.0
      #ch(3)=0.48                  #12.0/25.0
      #ch(4)=0.14222222222222222   #32.0/225.0
      #ch(5)=0.03333333333333333   #1.0/30.0
      #ch(6)=0.24                  #6.0/25.0
      #ct(1)=-0.006666666666666667 #-1.0/150.0
      #ct(2)=0.0
      #ct(3)=0.03                  #3.0/100.0
      #ct(4)=-0.21333333333333335  #-16.0/75.0
      #ct(5)=-0.05                 #-1.0/20.0
      #ct(6)=0.24                  #6.0/25.0
      #avepo=0.0

        #forall (mat=1:nTotalMat) sstor(mat)=(poroIni(mat)**3)/(ss(mat)*(1-poroIni(mat))**2)
        #mate=0.0
        #do node=1,nNode
          #column=0.0
          #oldColumn=0.0

          #do j=1,jti-1
              #do k=j+1,jti
                #do mat=1,nTotalMat
                  #column(j)=column(j)+((1-porosity(node,k,mat)*grid(node,k,mat)/&
                            #density(mat))+(porosity(node,k,mat)/density(nTotalMat+1)))
                #enddo
              #enddo
              #column(j)=column(j)+wDepth(node)/density(nTotalMat+1)
          #enddo

          #column(jti)=wDepth(node)/density(nTotalMat+1)

          #if (jti>1) then
            #do j=1,jti-1
              #do k=j+1,jti-1
                #do mat=1,nTotalMat
                  #oldColumn(j)=oldColumn(j)+((1-porosity(node,k,mat)*grid(node,k,mat)/ &
                              #density(mat))+(porosity(node,k,mat)/density(nTotalMat+1)))
                #enddo
              #enddo
              #do mat=1,nTotalMat
                    #oldColumn(j)=oldColumn(j)+(grid(node,jti,mat)*density(nTotalMat+1))
              #enddo
              #oldColumn(j)=oldColumn(j)+(wDepth(node)/density(nTotalMat+1))
            #enddo
          #else
            #oldColumn(jti)=wDepth(node)/density(nTotalMat+1)
          #endif
          #do j=1,jti
            #dp=column(j)-oldColumn(j)
            #sume=0.0
            #do mat=1,nTotalMat
              #mate(mat)=(1-porosity(node,j,mat))*grid(node,j,mat)
              #sume=sume+grid(node,j,mat)
            #enddo
            #avedummy=0.0
            #if (sume>0.0) then
              #do mat=1, nTotalMat
                #avedummy=avedummy+(dlog10(sstor(mat)))*grid(node,j,mat)/sume
              #enddo
            #endif
            #consum=10**avedummy
            #avepo=0.0
            #do mat=1,nTotalMat
              #avepo=avepo+(porosity(node,j,mat)*grid(node,j,mat))
              #if (isnan(avepo)) then
                #write(*,*) "avepo is nan"
                #write(*,*) mate(mat)
                #write(*,*) sume
                #write(*,*) avedummy
                #write(*,*)consum
                #write(*,*) porosity(node,j,mat)
                #write(*,*) grid(node,j,mat)
                #write(*,*)porosity(node,j,mat)*grid(node,j,mat)
                #write(*,*)"program stop"
                #stop
              #endif
            #enddo
            #avepor=avepo
##here start the compaction calculation
            #if (dp>0.0 .and.avepo>0.0.and.sume>0.0)then
                #x=0.0
                #h=0.001
                #tol=0.0000001
                #hold=h
                #do while (x<dp)
                  #xtemp=x
                  #z=-(1-avepo)*(avepo**3/(consum*(1-avepo)**2))
                  #f(1)=z
                  #ytemp=avepo
                  #temax=1000.0
                  #do while(temax>tol)
                        #do k=2,6
                            #x=xtemp+a(k)*h
                            #avepo=ytemp
                            #k1=k-1
                            #do l=1,k1
                                #avepo=avepo+h*b(k,l)*f(l)
                            #enddo
                          #z=-(1-avepo)*(avepo**3/(consum*(1-avepo)**2))
                            #f(k)=avepo
                        #enddo
                        #te=0.0
                        #avepo=ytemp
                        #do k=1,6
                            #avepo=avepo+h*ch(k)*f(k)
                            #te=te+h*ct(k)*f(k)
                        #enddo
                        #te=abs(te)
                        #temax=0.0
                        #if(temax<te) temax=te
                        #hold=h
                        #h=0.9*h*(tol/temax)**0.2
                        ##if(h>dp*0.5) h=0.01
                        #if(h>10000.0) h=10000.0
                        #if((dp-x)<h) h=dp-x
                        #if(h<0.001) then
                          #h=0.001
                          #teMax=tol
                          #hold=h
                        #end if
                    #enddo
                    #x=xtemp+hold
                #enddo
            #endif
            #dpdt=avepo-avepor
            #if (dpdt<0.0) dpdt=0.0
            #avepo=avepo-dpdt
            #if (dpdt>0.0)then
              #do mat=1,nTotalMat
                #dPor(mat)=(porosity(node,j,mat)*grid(node,j,mat))*avepo
                #dPor(mat)=dPor(mat)/grid(node,j,mat)
                #porosity(node,j,mat)=porosity(node,j,mat)-dPor(mat)
                #if (porosity(node,j,mat)<poroMin(mat))then
                  #porosity(node,j,mat)=poroMin(mat)
                #endif
              #enddo
              #do mat=1,nTotalMat
                #grid(node,j,mat)=mate(mat)/(1-porosity(node,j,mat))
                #if (node==221) write(198,*) j, totalTime,porosity(node,j,mat)
              #enddo
            #endif
          #enddo
        #enddo

#end subroutine compac