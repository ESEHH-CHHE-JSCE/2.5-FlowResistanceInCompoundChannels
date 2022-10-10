!***********************************************************************
!     水理公式集例題集（令和4年度版）　例題2.5　複断面河道の抵抗則
!***********************************************************************
      implicit none
      real(8),parameter::rho=1000.d0 !水の密度(kg/m3)
      real(8),parameter::g=9.8d0 !重力加速度(m/s2)
      real(8),parameter::eps0=1.d-8!打ち切り判定(次ステップの計算値-現時点の計算値）
      real(8),parameter::phi=0.1d0!収束安定化パラメータ（0～１）
      integer(4),parameter::kmax=10000 !収束計算の打ち切り回数
      integer(4)::i,k,ck,imax
      real(8)::zs,ib,f,fw,qs,as,us,beta,cr,cl,eps
      real(8),dimension(999)::zbi,zb0i,bi,hi,ai,ui,uni,qi,
     1                         ri,si,ni,tri,hwi,taus,tauds
      real(8),dimension(999)::swj,swdj
 

!断面情報の読み込み
      open(10,file='input-L3.dat',status='old')
      read(10,*)
      read(10,*) imax,zs,ib,f,fw
      read(10,*)
      do i=1,imax
      read(10,*) zb0i(i),bi(i),ni(i),hwi(i)
      end do

      close(10)

!分割断面の実質河床高(河床高＋樹木群の高さ)の定義
      do i=1,imax
      zbi(i)=zb0i(i)+hwi(i)
      end do

!分割断面境界の潤辺長(swdj)の計算
      do i=1,imax-1
      swdj(i)=zs-max(zbi (i),zbi (i+1))
      swj (i)=zs-max(zb0i(i),zb0i(i+1))
      swj (i)=swj(i)-swdj(i)
      end do

!分割断面の水深，断面積，潤辺，径深の計算,重力の流下方向成分の計算
      as=0.d0
      do i=1,imax
      hi(i)=max(zs-zbi(i),0.d0)
      ai(i)=hi(i)*bi(i)
      as=as+ai(i)
      if(i.eq.1) then
      si(i)=bi(i)+hi(i)+max(zb0i(i+1)-zb0i(i),0.d0)
      else if(i.eq.imax) then
      si(i)=bi(i)+hi(i)+max(zb0i(i-1)-zb0i(i),0.d0)
      else
      si(i)=bi(i)+max(zb0i(i-1)-zb0i(i),0.d0)
     1           +max(zb0i(i+1)-zb0i(i),0.d0)
      end if
      ri (i)=ai(i)/si(i)
      tri(i)=rho*g*ai(i)*ib
      end do

!準二次元解析***********************************************************

!初期条件
      qs =0.d0
      do i=1,imax
      ui (i)=0.d0
      uni(i)=0.d0
      end do

!ニュートン法による繰り返し計算開始
      k=0

 10   continue
      eps=0.d0

      do i=1,imax
      ck=0; cl=0.d0; cr=0.d0

!ゼロ割り防止
      if(ni (i).le.0.d0) ck=1
      if(tri(i).le.0.d0) ck=1
      if(ck.eq.1) then
      ui (i)=0.d0
      uni(i)=0.d0
      cycle
      end if

!分割断面の流速を計算
      if(i.eq.1) then
       if(zb0i(i+1).ge.zb0i(i)) cr=1.d0
       tauds(i)= rho*f *(ui(i)-ui(i+1))*abs(ui(i)-ui(i+1))*swdj(i  )
       taus (i)= rho*fw*ui(i)**2.d0*swj(i  )*cr
      else if(i.eq.imax) then
       if(zb0i(i-1).ge.zb0i(i)) cl=1.d0
       tauds(i)= rho*f *(ui(i)-ui(i-1))*abs(ui(i)-ui(i-1))*swdj(i-1)
       taus (i)= rho*fw*ui(i)**2.d0*swj(i-1)*cl
      else
       if(zb0i(i-1).ge.zb0i(i)) cl=1.d0
       if(zb0i(i+1).ge.zb0i(i)) cr=1.d0
       tauds(i)= rho*f *(ui(i)-ui(i-1))*abs(ui(i)-ui(i-1))*swdj(i-1)
     1          +rho*f *(ui(i)-ui(i+1))*abs(ui(i)-ui(i+1))*swdj(i  )
       taus (i)= rho*fw*ui(i)**2.d0*swj(i-1)*cl
     1          +rho*fw*ui(i)**2.d0*swj(i  )*cr
      end if

      if(tauds(i)+taus(i).ge.tri(i)) then
      uni(i)=0.d0
      else
      uni(i)=(ri(i)**(2.d0/3.d0)*ib**0.5d0/ni(i))
     1     *(1.d0-(tauds(i)+taus(i))/tri(i))**0.5d0
      end if

      eps=max(eps,abs(uni(i)-ui(i)))

      end do

      do i=1,imax
      ui(i)=phi*uni(i)+(1.d0-phi)*ui(i)
      end do

      k=k+1

      write(*,*) k,(ui(i),i=1,imax)
      if(k.ge.kmax) then
      write(*,*) "erro"
      stop
      end if

      if(eps.ge.eps0) goto 10

      do i=1,imax
      qs=qs+ai(i)*ui(i)
      end do
      us=qs/as

      beta=0.d0
      do i=1,imax
      beta=beta+1.1d0*ui(i)**2.d0*ai(i)
      end do
      beta=beta/(us**2.d0*as)
!***********************************************************************

!断面情報の読み込み
      open(10,file='ouput.dat',status='replace')
      write(10,*) "U",",","Q",",","β"
      write(10,100) us,qs,beta
      write(10,*) 
      write(10,*) "i",",","Ui",",","Ai",",","hi",",","Si",",","Ri"
      do i=1,imax
      write(10,101) i,ui(i),ai(i),hi(i),si(i),ri(i)
      end do
      write(10,*) 
      write(10,*) "j",",","Sw'",",","Sw"

      do i=1,imax-1
      write(10,102) i,swdj(i),swj(i)
      end do
      close(10)

 100  format(3f12.2)
 101  format(i8,5f12.2)
 102  format(i8,3f12.2)

      write(*,*) 'flowend'
      stop
      end program
