# -*- coding: utf-8 -*-

import modules.var as v


def compac():
    # It calculates sediment compaction
    for i in range(v.nTotalMat):
          sstor[mat] = (poroIni[mat]**3)/(ss[mat]*(1-poroIni[mat])**2)
          mate = 0.0
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