subroutine matsolFlow (dt,nNode,q,bw,typBC,head,nElem,in,tt,cn,nl,xCo,yCo)

!       solution of system of equations from fe mesh
!       adapted from Kinzelbach 1985 p. 290
      integer bw,nl,nNode,nElem
      real*8 xCo(nNode),yCo(nNode)
      integer in(nElem,3), typBC(nNode)
      real*8 head(nNode),bb(nNode), q(nNode)
      double precision tt(nElem),a(nNode,bw+1)
      double precision b(3),c(3),re(3,3)
      double precision pe(3,3)
      real*8 dt,cn
!       local variables
      real*8 rd,de,pi,su,h
      integer i,l,ii,jj,kd,ni,nj,ib,nk,ll,jb
      integer mb, kb, lb, ko

     
      a=0.0
      rd=1/dt
      bb=q
      forall(i=1:nNode, typBC(i)==0)
        a(i,nl+1)=1.0
        bb(i)=head(i)
      end forall
!       matrix assembly loop over elements
      do l=1,nElem
        b(1) = yCo(in(l,2))-yCo(in(l,3))
        b(2) = yCo(in(l,3))-yCo(in(l,1))
        b(3) = yCo(in(l,1))-yCo(in(l,2))
        c(1) = xCo(in(l,3))-xCo(in(l,2))
        c(2) = xCo(in(l,1))-xCo(in(l,3))
        c(3) = xCo(in(l,2))-xCo(in(l,1))
        de=(b(1)*c(2)-b(2)*c(1))*0.5
        pe(1,1)=(b(1)*b(1)*tt(l)+c(1)*c(1)*tt(l)+(b(1)*c(1)+c(1)*b(1))*0.0)/de*0.25
        pe(1,2)=(b(1)*b(2)*tt(l)+c(1)*c(2)*tt(l)+(b(1)*c(2)+c(1)*b(2))*0.0)/de*0.25
        pe(1,3)=(b(1)*b(3)*tt(l)+c(1)*c(3)*tt(l)+(b(1)*c(3)+c(1)*b(3))*0.0)/de*0.25
        pe(2,1)=(b(2)*b(1)*tt(l)+c(2)*c(1)*tt(l)+(b(2)*c(1)+c(2)*b(1))*0.0)/de*0.25
        pe(2,2)=(b(2)*b(2)*tt(l)+c(2)*c(2)*tt(l)+(b(2)*c(2)+c(2)*b(2))*0.0)/de*0.25
        pe(2,3)=(b(2)*b(3)*tt(l)+c(2)*c(3)*tt(l)+(b(2)*c(3)+c(2)*b(3))*0.0)/de*0.25
        pe(3,1)=(b(3)*b(1)*tt(l)+c(3)*c(1)*tt(l)+(b(3)*c(1)+c(3)*b(1))*0.0)/de*0.25
        pe(3,2)=(b(3)*b(2)*tt(l)+c(3)*c(2)*tt(l)+(b(3)*c(2)+c(3)*b(2))*0.0)/de*0.25
        pe(3,3)=(b(3)*b(3)*tt(l)+c(3)*c(3)*tt(l)+(b(3)*c(3)+c(3)*b(3))*0.0)/de*0.25
        re=de/12.0
        re(1,1)=de/6.0
        re(2,2)=de/6.0
        re(3,3)=de/6.0
! global matrix assembly
        ii=in(l,1)
        if (typBC(ii)/=0)then
          bb(ii)=bb(ii)-(pe(1,1))*head(in(l,1))
          if (typBC(in(l,1))/=0) then
            jj=in(l,1)-ii+nl+1
            if (jj>0) then
              a(ii,jj)=a(ii,jj)+(pe(1,1))*cn+ re(1,1)*rd
            endif
          endif
          bb(ii)=bb(ii)-(pe(1,2))*head(in(l,2))
          if (typBC(in(l,2))/=0) then
            jj=in(l,2)-ii+nl+1
            if (jj>0) then
              a(ii,jj)=a(ii,jj)+(pe(1,2))*cn+ re(1,2)*rd
            endif
          endif
          bb(ii)=bb(ii)-(pe(1,3))*head(in(l,3))
          if (typBC(in(l,3))/=0) then
            jj=in(l,3)-ii+nl+1
            if (jj>0) then
              a(ii,jj)=a(ii,jj)+(pe(1,3))*cn+ re(1,3)*rd
            endif
          endif
        endif

        ii=in(l,2)
        if (typBC(ii)/=0)then
          bb(ii)=bb(ii)-(pe(2,1))*head(in(l,1))
          if (typBC(in(l,1))/=0) then
            jj=in(l,1)-ii+nl+1
            if (jj>0) then
              a(ii,jj)=a(ii,jj)+(pe(2,1))*cn+ re(2,1)*rd
            endif
          endif
          bb(ii)=bb(ii)-(pe(2,2))*head(in(l,2))
          if (typBC(in(l,2))/=0) then
            jj=in(l,2)-ii+nl+1
            if (jj>0) then
              a(ii,jj)=a(ii,jj)+(pe(2,2))*cn+ re(2,2)*rd
            endif
          endif
          bb(ii)=bb(ii)-(pe(2,3))*head(in(l,3))
          if (typBC(in(l,3))/=0) then
            jj=in(l,3)-ii+nl+1
            if (jj>0) then
              a(ii,jj)=a(ii,jj)+(pe(2,3))*cn+ re(2,3)*rd
            endif
          endif
        endif
        ii=in(l,3)
        if (typBC(ii)/=0)then
          bb(ii)=bb(ii)-(pe(3,1))*head(in(l,1))
          if (typBC(in(l,1))/=0) then
            jj=in(l,1)-ii+nl+1
            if (jj>0) then
              a(ii,jj)=a(ii,jj)+(pe(3,1))*cn+ re(3,1)*rd
            endif
          endif
          bb(ii)=bb(ii)-(pe(3,2))*head(in(l,2))
          if (typBC(in(l,2))/=0) then
            jj=in(l,2)-ii+nl+1
            if (jj>0) then
              a(ii,jj)=a(ii,jj)+(pe(3,2))*cn+ re(3,2)*rd
            endif
          endif
          bb(ii)=bb(ii)-(pe(3,3))*head(in(l,3))
          if (typBC(in(l,3))/=0) then
            jj=in(l,3)-ii+nl+1
            if (jj>0) then
              a(ii,jj)=a(ii,jj)+(pe(3,3))*cn+ re(3,3)*rd
            endif
          endif
        endif
      enddo
    do i=1,nNode
    do j=1, bw+1
      write(963,*)a(i,j)   
    enddo
    write(964,*)bb(i)
  enddo    
      
! Solution of equation system
      kd=nl+1
! Factor matrix
      do ni=1,(nNode-1)
        pi=-a(ni,kd)
        if (pi/=0.0) then
          nj=ni+1
          ib=kd
          nk=ni+nl
          if (nk>nNode) nk=nNode
          do ll=nj,nk
            ib=ib-1
            h=a(ll,ib)
            if (h/=0.0) then
              h=h/pi
              a(ll,ib)=h
              jb=ib+1
              kb=ib+nk-ni
              lb=kd-ib
              do mb=jb,kb
                a(ll,mb)=a(ll,mb)+h*a(ni,lb+mb)
              enddo
            endif
          enddo
        else
          write(*,*) "error in matsolflow"
        endif
      enddo

!       solution using factored matrix

      do ni=2,nNode
        ko=ni-kd
        ib=1
        if (ko<=0) ib=1-ko
        nj=ib+ko
        su=0.0
        do jb=ib,nl
          su=su+a(ni,jb)*bb(nj)
          nj=nj+1
        enddo
        bb(ni)=bb(ni)+su
      enddo

!       back substitution
  write(7,*)bb
      bb(nNode)=bb(nNode)/a(nNode,kd)
      ll=kd+1
      ni=nNode

      do ib=2,nNode
        ni=ni-1
        nj=ni
        mb=bw
        if (ib<=nl) mb=nl+ib
        su=0.0
        do jb=ll,mb
          nj=nj+1
          su=su+a(ni,jb)*bb(nj)
        enddo
        bb(ni)=(bb(ni)-su)/a(ni,kd)
      enddo

! Equation solver finished
! Update concentrations

      forall(i=1:nNode, typBC(i)/=1) head(i)=head(i)+bb(i)

    end

