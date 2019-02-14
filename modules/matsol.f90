!_________________________________________________________________matsol
    subroutine matsol (dt,nNode,pm,bw,itn,cc,nElem,ve3t,velDi3,&
                     in,xd,yd,xy,th,nl,x,y)

!       solution of system of equations from fe mesh
!       adapted from Kinzelbach 1985 p. 290
      integer bw,nl,nNode,nElem
      real(KIND=8) ve3t(nElem),velDi3(nElem),x(nNode),y(nNode)
      integer in(nElem,3),itn(nNode)
      real(KIND=8) cc(nNode),bb(nNode),pm(nNode)
      double precision xd(nElem),yd(nElem),xy(nElem),a(nNode,bw+1)
      double precision b(3),c(3),re(3,3)
      double precision pe(3,3),ue(3,3)
      real(KIND=8) dt
!       local variables
      real(KIND=8) rd,ux,uy,de,th,pi,su,h
      integer i,j,l,ii,jj,fl,kj,kd,nu,ni,nj,ib,nk,ll,jb
      integer mb, kb, lb, ko

      fl=-1
      a=0
      rd=1/dt

!!OMP PARALLEL DO DEFAULT(SHARED) private(i,j) firstprivate(nNode,pm,fl,bw,itn)
      do i=1,nNode
        bb(i)=pm(i)
        if (fl<=0)then
          do j=1,bw
            a(i,j)=0.0
          enddo
        endif
        if (itn(i)==1)then
          a(i,nl+1)=1.0
          bb(i)=cc(i)
        endif
      enddo
!!OMP END PARALLEL DO
!       matrix assembly loop over elements
!!OMP PARALLEL DO DEFAULT(SHARED) private(l,ux,uy,b,c,de,pe,ue,re,i,ii,j,kj) firstprivate(nElem,velDi3,x,y,xd,yd,in,itn,fl,nl,th,rd)
      do l=1,nElem
!       flow components in each element in m/year

        ux=ve3t(l)*sin(velDi3(l)*0.017453292519943295)
        uy=ve3t(l)*cos(velDi3(l)*0.017453292519943295)

!       element matrix assembly

        b(1)=y(in(l,2))-y(in(l,3))
        b(2)=y(in(l,3))-y(in(l,1))
        b(3)=y(in(l,1))-y(in(l,2))
        c(1)=x(in(l,3))-x(in(l,2))
        c(2)=x(in(l,1))-x(in(l,3))
        c(3)=x(in(l,2))-x(in(l,1))

        de=(b(1)*c(2)-b(2)*c(1))*0.5

        pe(1,1)=(b(1)*b(1)*xd(l)+c(1)*c(1)*yd(l)+(b(1)*c(1)+c(1)*b(1))*xy(l))/de*0.25
        ue(1,1)=(b(1)*ux+c(1)*uy)/6.0
        pe(1,2)=(b(1)*b(2)*xd(l)+c(1)*c(2)*yd(l)+(b(1)*c(2)+c(1)*b(2))*xy(l))/de*0.25
        ue(1,2)=(b(2)*ux+c(2)*uy)/6.0
        pe(1,3)=(b(1)*b(3)*xd(l)+c(1)*c(3)*yd(l)+(b(1)*c(3)+c(1)*b(3))*xy(l))/de*0.25
        ue(1,3)=(b(3)*ux+c(3)*uy)/6.0

        pe(2,1)=(b(2)*b(1)*xd(l)+c(2)*c(1)*yd(l)+(b(2)*c(1)+c(2)*b(1))*xy(l))/de*0.25
        ue(2,1)=(b(1)*ux+c(1)*uy)/6.0
        pe(2,2)=(b(2)*b(2)*xd(l)+c(2)*c(2)*yd(l)+(b(2)*c(2)+c(2)*b(2))*xy(l))/de*0.25
        ue(2,2)=(b(2)*ux+c(2)*uy)/6.0
        pe(2,3)=(b(2)*b(3)*xd(l)+c(2)*c(3)*yd(l)+(b(2)*c(3)+c(2)*b(3))*xy(l))/de*0.25
        ue(2,3)=(b(3)*ux+c(3)*uy)/6.0

        pe(3,1)=(b(3)*b(1)*xd(l)+c(3)*c(1)*yd(l)+(b(3)*c(1)+c(3)*b(1))*xy(l))/de*0.25
        ue(3,1)=(b(1)*ux+c(1)*uy)/6.0
        pe(3,2)=(b(3)*b(2)*xd(l)+c(3)*c(2)*yd(l)+(b(3)*c(2)+c(3)*b(2))*xy(l))/de*0.25
        ue(3,2)=(b(2)*ux+c(2)*uy)/6.0
        pe(3,3)=(b(3)*b(3)*xd(l)+c(3)*c(3)*yd(l)+(b(3)*c(3)+c(3)*b(3))*xy(l))/de*0.25
        ue(3,3)=(b(3)*ux+c(3)*uy)/6.0
        re=de/12.0
        re(1,1)=de/6.0
        re(2,2)=de/6.0
        re(3,3)=de/6.0

! global matrix assembly

        do i=1,3
          ii=in(l,i)
          if (itn(ii)/=1)then
            do j=1,3
              kj=in(l,j)
              bb(ii)=bb(ii)-(pe(i,j)+ue(i,j))*cc(kj)
              if (fl<=0) then
                if (itn(kj)/=1) then
                  jj=kj-ii+nl+1
                  if (jj>0) then
                    a(ii,jj)=a(ii,jj)+(pe(i,j)+ue(i,j))*th+ re(i,j)*rd
                  endif
                endif
              endif
            enddo
          endif
        enddo
      enddo
!!OMP END PARALLEL DO
! Solution of equation system
      if (fl<=0) then
        kd=nl+1
! Factor matrix
        nu=nNode-1
        do ni=1,nu
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
        call error(6,ni)
          endif
        enddo
      endif

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
!$OMP PARALLEL WORKSHARE private(i) firstprivate(nNode)
      forall (i=1:nNode,itn(i)/=1)cc(i)=cc(i)+bb(i)
!$OMP END PARALLEL WORKSHARE
      if (fl<=0) fl=1

    end
