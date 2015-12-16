c------------------------------------------------------------------------------
c--- NBO2MOLDEN: a utility to transform NBO graphical data to MOLDEN format.
c--- Written in FORTRAN 77 and a bit of FORTRAN 90.
c---
c--- Tested for NBO 3.0, 5.x, and 6.x.
c---
c--- Updated:
c--- Ver.1.0.3, 12/07/2015, Bug fix for d-, f-, and g-functions.
c--- Ver.1.0.4, 12/16/2015, Plot file 41 by NBO6.
c---
c--- E-mail: qcband@gmail.com
c------------------------------------------------------------------------------

      program NBO2MOLDEN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*27 fmolden
      character*22 title
      character*2 fext

      character*10 dt
      character*5 ver

      ifbs=31
      ifob=32
      ifmo=40

c------------------------------------------------------------------------------
c--- head
c------------------------------------------------------------------------------
      ver="1.0.4"
      dt="12/16/2015"
      call head(ver,dt)

c------------------------------------------------------------------------------
c--- define orbital type
c------------------------------------------------------------------------------
      call obrtyp(fext,title)

c------------------------------------------------------------------------------
c--- define NBO graphical data
c------------------------------------------------------------------------------
      call nbodata(fext,fmolden)

c------------------------------------------------------------------------------
c--- read *.31
c------------------------------------------------------------------------------
      call rdf31(ist)
      if(ist.ne.0)goto 9999

c------------------------------------------------------------------------------
c--- read *.3x MO
c------------------------------------------------------------------------------
      call rdfmo(ist)
      if(ist.ne.0)goto 9999

c------------------------------------------------------------------------------
c--- write molden
c------------------------------------------------------------------------------
      OPEN(ifmo,FILE=trim(fmolden))
      rewind(ifmo)
      call wrcoord
      call wrgto
      call wrmo

c------------------------------------------------------------------------------
c--- finished
c------------------------------------------------------------------------------
      write(*,"(//,1x,55('-'),/,2x
     *'The MOLDEN file has been generated successfully!')")
      write(*,"(2x,'File Name = ',a,/,1x,55('-'))")trim(fmolden)
      close (ifmo)

9999  continue
      close (ifbs)
      close (ifob)

      call estop

      end

c------------------------------------------------------------------------------
c--- read *.3x MO
c------------------------------------------------------------------------------
      subroutine rdfmo(ist)
      implicit double precision (a-h,o-z)
      parameter(maxexp=5000)
      common/BASIC/NATOMS,NSHELL,NEXP,NBAS,labsph,maxlq
      common/MOORB/ititle,iab,iocc,occ(2,maxexp)
      character*23 tmp,tag00,tag37,tag39

      ifob=32
      ist=0
      tag00=' ----------------------'
      tag37=' NBOs in the AO basis: '
      tag39=' NLMOs in the AO basis:'

c--- test
      rewind (ifob)
      read(ifob,"(//,a23)")tmp
      if(tmp.ne.tag00)then
c--- USCF or ROSCF
c--- .32, and .33: no title (NBO-3 only), no occupations, no alpha/beta
c--- .34-.36, .38, and .40-.41: no occupations
        ititle=0
        iab=0
        iocc=0
        goto 100
      else
        ititle=1
      end if

      read(ifob,"(a23)")tmp
      if(tmp(1:11).eq.' ALPHA SPIN')then
        iab=1    ! open
      else
        iab=0
      end if
      rewind (ifob)
      read(ifob,"(/,a23)")tmp
      if(tmp.eq.tag37.or.tmp.eq.tag39)then
        iocc=1
      else
        iocc=0
      end if

100   continue
c      write(*,"(3i3)")ititle,iab,iocc

c--- occ.
      if(iocc.eq.0)then
        do i=1,NBAS
          occ(1,i)=0.d0
        end do
      else
        call lines(NBAS,5,nline,last)
        if(iab.eq.0)then
          rewind (ifob)
          read(ifob,"(//)",err=1000,END=1000)
          do i=1,NBAS
            do j=1,nline
              read(ifob,"()",err=1100,END=1100)
            end do
          end do
          read(ifob,"(1x,5f15.9)",err=1200,END=1200)(occ(1,i),i=1,NBAS)
        else
c--- alpha
          rewind (ifob)
200       read(ifob,"(a23)",err=2000,END=2000)tmp
          if(tmp(1:11).ne.' ALPHA SPIN')goto 200
          do i=1,NBAS
            do j=1,nline
              read(ifob,"()",err=2100,END=2100)
            end do
          end do
          read(ifob,"(1x,5f15.9)",err=2200,END=2200)(occ(1,i),i=1,NBAS)
c--- beta
          rewind (ifob)
400       read(ifob,"(a23)",err=3000,END=3000)tmp
          if(tmp(1:11).ne.' BETA  SPIN')goto 400
          do i=1,NBAS
            do j=1,nline
              read(ifob,"()",err=3100,END=3100)
            end do
          end do
          read(ifob,"(1x,5f15.9)",err=3200,END=3200)(occ(2,i),i=1,NBAS)
        end if
      end if
c--- DEBUG
c      write(9,"(1x,5f15.9)")(occ(1,i),i=1,NBAS)
c      write(9,"(1x,5f15.9)")(occ(2,i),i=1,NBAS)
c---

      return

1000  ist=1
      write(*,"(' Error = #NBO.100')")
      return
1100  ist=1
      write(*,"(' Error = #NBO.110')")
      return
1200  ist=1
      write(*,"(' Error = #NBO.200')")
      return
2000  ist=1
      write(*,"(' Error = #NBO.200')")
      return
2100  ist=1
      write(*,"(' Error = #NBO.210')")
      return
2200  ist=1
      write(*,"(' Error = #NBO.220')")
      return
3000  ist=1
      write(*,"(' Error = #NBO.300')")
      return
3100  ist=1
      write(*,"(' Error = #NBO.310')")
      return
3200  ist=1
      write(*,"(' Error = #NBO.320')")
      return

      end

c------------------------------------------------------------------------------
c--- write MO
c------------------------------------------------------------------------------
      subroutine wrmo
      implicit double precision (a-h,o-z)
      parameter(maxshell=800,maxexp=5000)
      common/BASIC/NATOMS,NSHELL,NEXP,NBAS,labsph,maxlq
      common/MOORB/ititle,iab,iocc,occ(2,maxexp)
      common/BASIS/iatm(maxshell),lq(maxshell),ncomp(maxshell),
     *nptr(maxshell),nprim(maxshell),label(maxexp),map(maxexp)
      dimension fmo(maxexp)
      character*11 tmp

      ifob=32
      ifmo=40
      write(ifmo,"('[MO]')")

c---  total-E
      tote=0.d0
      do i=1,NBAS
        tote=tote+occ(1,i)
      end do
      if(iab.eq.1)then
        do i=1,NBAS
          tote=tote+occ(1,i)
        end do
      end if
      if(iocc.eq.1)then
        write(*,"('  The total charge is ',f7.1,'.')")tote
      else
        write(*,"('  Writing......')")
      end if

      if(iab.eq.0)then
        rewind (ifob)
        if(ititle.eq.1)read(ifob,"(//)")
        do i=1,NBAS
          read(ifob,"(1x,5f15.9)")(fmo(j),j=1,NBAS)
          write(ifmo,"(' Ene= 0.0')")
          write(ifmo,"(' Spin= Alpha')")
          write(ifmo,"(' Occup= ',f11.6)")occ(1,i)
          do j=1,NBAS
            write(ifmo,"(i4,e18.9)")j,fmo(map(j))
          end do
        end do
      else
        rewind (ifob)
c--- alpha
100     read(ifob,"(a11)")tmp
        if(tmp.ne.' ALPHA SPIN')goto 100
        do i=1,NBAS
          read(ifob,"(1x,5f15.9)")(fmo(j),j=1,NBAS)
          write(ifmo,"(' Ene= 0.0')")
          write(ifmo,"(' Spin= Alpha')")
          write(ifmo,"(' Occup= ',f11.6)")occ(1,i)
          do j=1,NBAS
          end do
        end do
c--- beta
200     read(ifob,"(a11)")tmp
        if(tmp.ne.' BETA  SPIN')goto 200
        do i=1,NBAS
          read(ifob,"(1x,5f15.9)")(fmo(j),j=1,NBAS)
          write(ifmo,"(' Ene= 0.0')")
          write(ifmo,"(' Spin= Beta')")
          write(ifmo,"(' Occup= ',f11.6)")occ(2,i)
          do j=1,NBAS
            write(ifmo,"(i4,e18.9)")j,fmo(map(j))
          end do
        end do
      end if
c

      return
      end

c------------------------------------------------------------------------------
c--- count the number of lines
c------------------------------------------------------------------------------
      subroutine lines(n1,n2,nline,last)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      nline=n1/n2
      last=n2
      if(mod(n1,n2).ne.0)then
        last=mod(n1,n2)
        nline=nline+1
      end if

      return
      end

c------------------------------------------------------------------------------
c--- read *.31
c------------------------------------------------------------------------------
      subroutine rdf31(ist)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(maxnatm=100,maxshell=800,maxexp=5000)
      common/COORD/xyz(3,maxnatm),iz(maxnatm)
      common/BASIC/NATOMS,NSHELL,NEXP,NBAS,labsph,maxlq
      common/BASIS/iatm(maxshell),lq(maxshell),ncomp(maxshell),
     *nptr(maxshell),nprim(maxshell),label(maxexp),map(maxexp),
     *EX(maxexp),CS(maxexp),CP(maxexp),CD(maxexp),CF(maxexp),CG(maxexp)

      ifbs=31
      ist=0

      rewind (ifbs)
      read(ifbs,"(///,1x,3i6,/)",err=100)NATOMS,NSHELL,NEXP
c--- check
      if(NATOMS.gt.maxnatm)goto 200
      if(NSHELL.gt.maxshell)goto 210

c--- the default unit is Ang.
      do i=1,NATOMS
        read(ifbs,"(i5,3f14.9)",err=300)iz(i),xyz(1:3,i)
      end do
      read(ifbs,"()")

c--- read basis set
      NBAS=0
      ne=0
      maxlq=0
      
      do i=1,NSHELL
        read(ifbs,"(1x,4i6)",err=400)iatm(i),ncomp(i),nptr(i),nprim(i)
        read(ifbs,"(1x,10i6)",err=500)label(NBAS+1:NBAS+ncomp(i))
        call lquant(lq(i),ncomp(i),labsph)
        if(lq(i).eq.-2)goto 600
c--- check.
c--- There is a BUG in Gaussian-NBO 5 for g-functions: g is treated as f.
        lq2=max(label(NBAS+1),label(NBAS+ncomp(i)))/100
        if(lq2.ne.abs(lq(i)))goto 700
c--- reorder label and save to map
        call moldenorder(label,map,lq(i),NBAS,ncomp(i),labsph,ist)
        if(ist.ne.0)goto 710
        maxlq=max(maxlq,abs(lq(i)))
        NBAS=NBAS+ncomp(i)
        ne=ne+nprim(i)
      end do
c--- check
      if(NBAS.gt.maxexp)goto 800
      if(ne.ne.NEXP)goto 900
c--- DEBUG
c      write(9,"(10i4)")label(1:NBAS)
c      write(9,"(10i4)")label(map(1:NBAS))
c---

c--- read basis functions
c--- EX
      read(ifbs,"()",err=1000,END=1000)
      read(ifbs,"(2x,4e18.9)",err=1000,END=1000)(EX(i),i=1,NEXP)
c--- DEBUG
c      write(9,"(2x,4e18.9)")(EX(i),i=1,NEXP)
c---
c--- CS
      if(maxlq.ge.0)then
        read(ifbs,"()",err=1100,END=1100)
        read(ifbs,"(2x,4e18.9)",err=1100,END=1100)(CS(i),i=1,NEXP)
        call denorm(EX,CS,0,NEXP)
c--- DEBUG
c        write(9,"(2x,4e18.9)")(CS(i),i=1,NEXP)
c---
      end if
c--- CP
      if(maxlq.ge.1)then
        read(ifbs,"()",err=1200,END=1200)
        read(ifbs,"(2x,4e18.9)",err=1200,END=1200)(CP(i),i=1,NEXP)
        call denorm(EX,CP,1,NEXP)
c--- DEBUG
c        write(9,"(2x,4e18.9)")(CP(i),i=1,NEXP)
c---
      end if
c--- CD
      if(maxlq.ge.2)then
        read(ifbs,"()",err=1300,END=1300)
        read(ifbs,"(2x,4e18.9)",err=1300,END=1300)(CD(i),i=1,NEXP)
        call denorm(EX,CD,2,NEXP)
c--- DEBUG
c        write(9,"(2x,4e18.9)")(CD(i),i=1,NEXP)
c---
      end if
c--- CF
      if(maxlq.ge.3)then
        read(ifbs,"()",err=1400,END=1400)
        read(ifbs,"(2x,4e18.9)",err=1400,END=1400)(CF(i),i=1,NEXP)
        call denorm(EX,CF,3,NEXP)
c--- DEBUG
c        write(9,"(2x,4e18.9)")(CF(i),i=1,NEXP)
c---
      end if
c--- CG
      if(maxlq.ge.4)then
        read(ifbs,"()",err=1500,END=1500)
        read(ifbs,"(2x,4e18.9)",err=1500,END=1500)(CG(i),i=1,NEXP)
        sumg=0.d0
        do i=1,NEXP
          sumg=sumg+abs(CG(i))
        end do
        if(sumg.lt.1.d-8)goto 1600
        call denorm(EX,CG,4,NEXP)
c--- DEBUG
c        write(9,"(2x,4e18.9)")(CG(i),i=1,NEXP)
c---
      end if

      return

100   ist=1
      write(*,"(' Error = #31.01')")
      return
200   ist=1
      write(*,"(' Error = #31.020')")
      write(*,"(' NATOMS=',i4,' > MAXNATM=',i4)")NATOMS,maxnatm
      return
210   ist=1
      write(*,"(' Error = #31.021')")
      write(*,"(' NSHELL=',i4,' > MAXSHELL=',i4)")NSHELL,maxshell
      return
300   ist=1
      write(*,"(' Error = #31.03')")
      return
400   ist=1
      write(*,"(' Error = #31.04')")
      return
500   ist=1
      write(*,"(' Error = #31.05')")
      return
600   ist=1
      write(*,"(' Error = #31.06')")
      write(*,"(' Unknown basis type. NCOMP=',i2)")ncomp(i)
      return
700   ist=1
      write(*,"(' Error = #31.070')")
      write(*,"(
     *' The basis types do not match. It seems that the .31 file is',/,
     *' generated by Gaussian-NBO5. This is a BUG of Gaussian.')")
      return
710   ist=1
      write(*,"(' Error = #31.071')")
      return
800   ist=1
      write(*,"(' Error = #31.08')")
      write(*,"(' NBAS=',i4,' > MAXEXP=',i4)")NBAS,maxexp
      return
900   ist=1
      write(*,"(' Error = #31.09')")
      write(*,"(' NE=',i4,' .ne. NEXP=',i4)")ne,NEXP
      return
1000  ist=1
      write(*,"(' Error = #31.10',/,' No EXP part.')")
      return
1100  ist=1
      write(*,"(' Error = #31.11',/,' No CS part.')")
      return
1200  ist=1
      write(*,"(' Error = #31.12',/,' No CP part.')")
      return
1300  ist=1
      write(*,"(' Error = #31.13',/,' No CD part.')")
      return
1400  ist=1
      write(*,"(' Error = #31.14',/,' No CF part.')")
      return
1500  ist=1
      write(*,"(' Error = #31.15',/,' No CG part.')")
      write(*,"(' It seems that the .31 file is generated by NBO3.')")
      return
1600  ist=1
      write(*,"(' Error = #31.16')")
      write(*,"(
     *' The factors in the CG part are all zero. It seems that the',/,
     *' .31 file is generated by by Gaussian-NBO5. This is a BUG of',/,
     *' Gaussian.')")
      return
      end

c------------------------------------------------------------------------------
c--- de-normalization
c------------------------------------------------------------------------------
      subroutine denorm(E,C,lq,n)
      implicit double precision (a-h,o-z)
      dimension E(*),C(*)

      do i=1,n
        if(abs(c(i)).lt.1.d-10)cycle
        cf=1.d0
        call fnorm(lq,e(i),cf)
        c(i)=c(i)/cf
      end do

      return
      end

c------------------------------------------------------------------------------
c--- calculate the normalization factor for GTO(L)
c------------------------------------------------------------------------------
      subroutine fnorm(lq,ex,cf)
      implicit double precision (a-h,o-z)

      pi=acos(-1.d0)
      pi3=pi**3.d0

c--- unnormalize primitives
c--- Normal^4 = 2^n1 * a^n2 / (pi^3 * nf^2)
c---   n1=3+4*L; n2=3+2*L, nf=(2L-1)!!
      call power(lq,n1,n2,nf)
      f = (2.d0**dble(n1)) * (ex**dble(n2)) / (pi3 * dble(nf))
      cf=cf*sqrt(sqrt(f))

      return
      end

c------------------------------------------------------------------------------
c--- get power(n1,n2,nf) for GTO(L) normalization
c--- n1=3+4*L; n2=3+2*L, nf=(2L-1)
c------------------------------------------------------------------------------
      subroutine power(lq,n1,n2,nf)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      n1=0
      n2=0
      nf=0
      select case(lq)
        case(0)  ! s
          n1=3
          n2=3
          nf=1
        case(1)  ! p
          n1=7
          n2=5
          nf=1
        case(2)  ! d
          n1=11
          n2=7
          nf=9
        case(3)  ! f
          n1=15
          n2=9
          nf=225
        case(4)  ! g
          n1=19
          n2=11
          nf=11025
      end select

      return
      end

c------------------------------------------------------------------------------
c--- reorder labels according to MOLDEN format
c
c  SP: S Px Py Pz
c
c  5D: D 0, D+1, D-1, D+2, D-2
c  6D: xx, yy, zz, xy, xz, yz
c
c  7F: F 0, F+1, F-1, F+2, F-2, F+3, F-3
c 10F: xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
c
c  9G: G 0, G+1, G-1, G+2, G-2, G+3, G-3, G+4, G-4
c 15G: xxxx yyyy zzzz xxxy xxxz yyyx yyyz zzzx zzzy xxyy xxzz yyzz xxyz yyxz zzxy
c------------------------------------------------------------------------------
      subroutine moldenorder(label,map,lq,idx,nc,lsc,ist)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension label(*),map(*)
      dimension libp(3),libl(4),libd0(5),libd1(6),libf0(7),libf1(10),
     *libg0(9),libg1(15)
      data libp/101,102,103/,libl/1,101,102,103/,
     *libd0/255,252,253,254,251/,libd1/201,204,206,202,203,205/,
     *libf0/351,352,353,354,355,356,357/,
     *libf1/301,307,310,304,302,303,306,309,308,305/,
     *libg0/451,452,453,454,455,456,457,458,459/,
     *libg1/401,411,415,402,403,407,412,410,414,404,406,413,405,408,409/

      do i=1,nc
        map(idx+i)=0
      end do
      ist=0

c--- spherical
      if(lsc.eq.0)then
        if(lq.eq.0)then           ! s
          map(idx+1)=idx+1
        else if(lq.eq.1)then      ! p
          do k=1,nc
            do j=idx+1,idx+nc
              if(label(j).eq.libp(k).or.label(j).eq.(libp(k)+50))then
                map(idx+k)=j
                exit
              end if
            end do
          end do
        else if(lq.eq.-1)then     ! sp
          do k=1,nc
            do j=idx+1,idx+nc
              if(label(j).eq.libl(k).or.label(j).eq.(libl(k)+50))then
                map(idx+k)=j
                exit
              end if
            end do
          end do
        else if(lq.eq.2)then     ! d
          do k=1,nc
            do j=idx+1,idx+nc
              if(label(j).eq.libd0(k))then
                map(idx+k)=j
                exit
              end if
            end do
          end do
        else if(lq.eq.3)then     ! f
          do k=1,nc
            do j=idx+1,idx+nc
              if(label(j).eq.libf0(k))then
                map(idx+k)=j
                exit
              end if
            end do
          end do
        else if(lq.eq.4)then     ! g
          do k=1,nc
            do j=idx+1,idx+nc
              if(label(j).eq.libg0(k))then
                map(idx+k)=j
                exit
              end if
            end do
          end do
        end if

c--- cartesian
      else
        if(lq.eq.0)then           ! s
          map(idx+1)=idx+1
        else if(lq.eq.1)then      ! p
          do k=1,nc
            do j=idx+1,idx+nc
              if(label(j).eq.libp(k).or.label(j).eq.(libp(k)+50))then
                map(idx+k)=j
                exit
              end if
            end do
          end do
        else if(lq.eq.-1)then     ! sp
          do k=1,nc
            do j=idx+1,idx+nc
              if(label(j).eq.libl(k).or.label(j).eq.(libl(k)+50))then
                map(idx+k)=j
                exit
              end if
            end do
          end do
        else if(lq.eq.2)then     ! d
          do k=1,nc
            do j=idx+1,idx+nc
              if(label(j).eq.libd1(k))then
                map(idx+k)=j
                exit
              end if
            end do
          end do
        else if(lq.eq.3)then     ! f
          do k=1,nc
            do j=idx+1,idx+nc
              if(label(j).eq.libf1(k))then
                map(idx+k)=j
                exit
              end if
            end do
          end do
        else if(lq.eq.4)then     ! g
          do k=1,nc
            do j=idx+1,idx+nc
              if(label(j).eq.libg1(k))then
                map(idx+k)=j
                exit
              end if
            end do
          end do
        end if
      end if

c--- check map
      do i=1,nc
        if(map(idx+i).eq.0)then
          ist=1
          goto 9999
        end if
      end do

9999  continue

      return
      end

c------------------------------------------------------------------------------
c--- L-quantum number
c--- s=0, p=1, sp=-1, d=2, f=3, g=4, unknown=-2
c--- lab= 0 (spherical) / 1 (cartesian)
c------------------------------------------------------------------------------
      subroutine lquant(lq,nc,lab)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      lab=0    ! spherical
      select case(nc)
        case(1)
          lq=0
        case(3)
          lq=1
        case(4)
          lq=-1
        case(5)
          lq=2
        case(6)
          lq=2
          lab=1    ! cartesian
        case(7)
          lq=3
        case(10)
          lq=3
          lab=1    ! cartesian
        case(9)
          lq=4
        case(15)
          lq=4
          lab=1    ! cartesian
        case default
          lq=-2
      end select

      return
      end

c------------------------------------------------------------------------------
c--- write coordinates
c------------------------------------------------------------------------------
      subroutine wrcoord
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(maxnatm=100)
      common/COORD/xyz(3,maxnatm),iz(maxnatm)
      common/BASIC/NATOMS,NSHELL,NEXP,NBAS,labsph,maxlq
      parameter (max_za=103)
      character*2 ATOMLIB(max_za)
      data (atomlib(i),i=1,50) /
     1'H ','HE','LI','BE','B ','C ','N ','O ','F ','NE',
     2'NA','MG','AL','SI','P ','S ','CL','AR','K ','CA',
     3'SC','TI','V ','CR','MN','FE','CO','NI','CU','ZN',
     4'GA','GE','AS','SE','BR','KR','RB','SR','Y ','ZR',
     5'NB','MO','TC','RU','RH','PD','AG','CD','IN','SN'/
      data (atomlib(i),i=51,100) /
     1'SB','TE','I ','XE','CS','BA','LA','CE','PR','ND',
     2'PM','SM','EU','GD','TB','DY','HO','ER','TM','YB',
     3'LU','HF','TA','W ','RE','OS','IR','PT','AU','HG',
     4'TL','PB','BI','PO','AT','RN','FR','RA','AC','TH',
     5'PA','U ','NP','PU','AM','CM','BK','CF','ES','FM'/
      data (atomlib(i),i=101,103) /
     1'MD','NO','LR'/

      ifmo=40

      write(ifmo,"('[MOLDEN FORMAT]')")
      write(ifmo,"('[Atoms] Angs')")

      do i=1,NATOMS
        write(ifmo,"(a2,2i5,1x,3f20.10)")
     *  atomlib(iz(i)),i,iz(i),xyz(1:3,i)
      end do

      return
      end

c------------------------------------------------------------------------------
c--- write gto basis set
c------------------------------------------------------------------------------
      subroutine wrgto
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(maxnatm=100,maxshell=800,maxexp=5000)
      common/BASIC/NATOMS,NSHELL,NEXP,NBAS,labsph,maxlq
      common/BASIS/iatm(maxshell),lq(maxshell),ncomp(maxshell),
     *nptr(maxshell),nprim(maxshell),label(maxexp),map(maxexp),
     *EX(maxexp),CS(maxexp),CP(maxexp),CD(maxexp),CF(maxexp),CG(maxexp)
      character*2 bstyp(3),clq
      data bstyp/'5D','7F','9G'/

      ifmo=40

c--- spherical basis set
      if(labsph.eq.0)then
        do i=1,3
          if((i+1).le.maxlq)write(ifmo,"('[',a2,']')")bstyp(i)
        end do
      end if

c--- gto
      write(ifmo,"('[GTO]')")
      do i=1,NATOMS
        write(ifmo,"(i3,' 0')")i
        do j=1,NSHELL
          if(iatm(j).ne.i)cycle
          write(ifmo,"(1x,a2,1x,i3,3x,'1.0')")clq(lq(j)),nprim(j)
          do k=nptr(j),nptr(j)+nprim(j)-1
            select case (lq(j))
              case(-1)
                write(ifmo,"(e18.10,2e14.6)")EX(k),CS(k),CP(k)
              case(0)
                write(ifmo,"(e18.10,e14.6)")EX(k),CS(k)
              case(1)
                write(ifmo,"(e18.10,e14.6)")EX(k),CP(k)
              case(2)
                write(ifmo,"(e18.10,e14.6)")EX(k),CD(k)
              case(3)
                write(ifmo,"(e18.10,e14.6)")EX(k),CF(k)
              case(4)
                write(ifmo,"(e18.10,e14.6)")EX(k),CG(k)
            end select
          end do
        end do
        write(ifmo,"()")
      end do

      return
      end

c------------------------------------------------------------------------------
c--- int LQ --> char LQ
c------------------------------------------------------------------------------
      character*2 function clq(l)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      select case(l)
        case(-1)
          clq='sp'
        case(0)
          clq='s '
        case(1)
          clq='p '
        case(2)
          clq='d '
        case(3)
          clq='f '
        case(4)
          clq='g '
      end select

      return
      end

c------------------------------------------------------------------------------
c--- define orbital type
c------------------------------------------------------------------------------
      subroutine obrtyp(fext,title)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(ntp=10)
      character*2 fext
      character*22 nbo(ntp),title
      data nbo/
     *'PNAOs in the AO basis ',
     *'NAOs in the AO basis  ',
     *'PNHOs in the AO basis ',
     *'NHOs in the AO basis  ',
     *'PNBOs in the AO basis ',
     *'NBOs in the AO basis  ',
     *'PNLMOs in the AO basis',
     *'NMLOs in the AO basis ',
     *'MOs in the AO basis   ',
     *'NOs in the AO basis   '/

      write(*,"(/,
     *' Which NBO file do you want to convert (1~',i2,'):')")ntp
      do inbo=1,ntp
        if(inbo.eq.6)then
          write(*,800)inbo,nbo(inbo),inbo+31
        else
          write(*,820)inbo,nbo(inbo),inbo+31
        end if
      end do
100   write(*,"(' > ',$)")
      read(*,"(a2)")fext
c--- default
      if(len_trim(fext).eq.0)fext=' 6'
      read(fext,*,err=200)i
      if(i.lt.1 .or. i.gt.ntp)then
        write(*,900)
        goto 100
      end if

      inbo=i+31
      write(fext,"(i2)")inbo
      write(*,"(/,' You select',/,2x,i1,') ',a22,
     *' ...... (*.',i2,')',/,1x,55('-'),//)")i,nbo(i),inbo

      return

200   write(*,900)
      goto 100

800   format(2x,i2,') ',a22,' ...... (*.',i2,') <--- (Default)')
820   format(2x,i2,') ',a22,' ...... (*.',i2,')')
900   format(/,' Unknown selection. Please try again.',/)
      return
      end

c------------------------------------------------------------------------------
c--- define NBO graphical data
c--- fext='32'~'41'; fmd: MOLDEN file name (length=leng<27)
c------------------------------------------------------------------------------
      subroutine nbodata(fext,fmd)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*23 fbs,fob
      character*20 ftmp
      character*2 fext
      character*27 fmd

      ifbs=31
      ifob=32

100   write(*,"(
     *' Input the base name of NBO graphical data (no extension)',/,
     *' within 20 characters:')")
      write(*,"(' > ',$)")
      read(*,"(a20)")ftmp
      lstr=nonspace(ftmp)
      lend=LEN_TRIM(ftmp)
      if(lend.eq.0)then                 ! use default file name
        lstr=1
        lend=4
        ftmp(lstr:lend)='FILE'
      end if
      fbs=ftmp(lstr:lend)//'.31'
      fob=ftmp(lstr:lend)//'.'//fext
      fmd=ftmp(lstr:lend)//'.molden'
      leng=lend-lstr+4

      open(ifbs,file=fbs(1:leng),status='old',err=110)
      goto 200
110   write(*,800)trim(fbs)
      goto 100

200   open(ifob,file=fob(1:leng),status='old',err=210)
      goto 300
210   write(*,820)trim(fob)
      goto 100

300   write(*,"(/,' The following files have been found.',/,
     *2(2x,a,/),1x,55('-'),//)")trim(fbs),trim(fob)

      return
800   format(/,' Can not find the basis set file ',a,/,
     *' Please try again.',/)
820   format(/,' Can not find the NBO graphical data file ',a,/,
     *' Please try again.',/)
      end

c------------------------------------------------------------------------------
c--- read an <ENTER> and stop
c------------------------------------------------------------------------------
      subroutine estop
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      write(*,"(//,' Press <ENTER> to exit',/)")
      read(*,*)

      stop

      return
      end

c------------------------------------------------------------------------------
c--- length of a string without the first and last spaces.
c------------------------------------------------------------------------------
      function nonspace(string)
      implicit double precision (a-h,o-z)
      character*(*) string
      character*1 space

      space=' '
      length=LEN_TRIM(string)
      if(length.eq.0) then
       i=1
      else
       do i=1,length
         if(string(i:i).ne.space) goto 20
       end do
      endif

20    nonspace=i

      return
      end

c------------------------------------------------------------------------------
c--- head
c------------------------------------------------------------------------------
      subroutine head(ver,dt)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*10 dt
      character*5 ver

      write(*,"(//,1x,46('*'),/
     *'           * * *   NBO2MOLDEN   * * *',/,
     *'          Version ',a5,',',4x,a10,//,
     *' Transform NBO graphical data to MOLDEN format',/,
     *1x,46('*'),/)")ver,dt

      return
      end