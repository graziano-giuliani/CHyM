        module mod_hourdb

        use mod_time
        use mod_museo

        contains

	subroutine validazione(lat,lon,source,p,n)
	integer source(1)
	real lat(1),lon(1),p(1)
	common /timefl/ ora,giorno,mese,anno,idate,logf
	integer ora,giorno,mese,anno,idate ; logical logf 
	integer scart(100) ; data scart /100*0/ ; save scart
	m=n
	do i=1,m 
	   if (source(i).eq.1) then
	      if (abs(lat(i)-41.9611).lt.1.0e-10.and. &
		  abs(lon(i)-13.6517).lt.1.0e-10.and.anno.eq.2000) then
	          if (scart(1).eq.0) then
	             if (logf)write(6,'(10x,a,2f10.5,i5,f10.1,a)') 'Discarded ', &
		            lat(i),lon(i),source(i),p(i) ,' for the whole year'
		     scart(1)=1
		  endif
	          lat(i)=lat(n)
	          lon(i)=lon(n)
	          source(i)=source(n)
	          p(i)=p(n)
	          n=n-1
	      else if (abs(lat(i)-42.0325).lt.1.0e-10.and.   &
		       abs(lon(i)-13.8806).lt.1.0e-10.and.anno.eq.2001) then
	          if (scart(2).eq.0) then
	             if (logf)write(6,'(10x,a,2f10.5,i5,f10.1,a)') 'Discarded ', &
		            lat(i),lon(i),source(i),p(i) ,' for the whole year'
		     scart(2)=1
		  endif
	          lat(i)=lat(n)
	          lon(i)=lon(n)
	          source(i)=source(n)
	          p(i)=p(n)
	          n=n-1
	      else if (abs(lat(i)-42.4436).lt.1.0e-10.and. &
		       abs(lon(i)-13.9033).lt.1.0e-10.and.anno.eq.2003) then
	          if (scart(3).eq.0) then
	             if (logf)write(6,'(10x,a,2f10.5,i5,f10.1,a)') 'Discarded ', &
		            lat(i),lon(i),source(i),p(i) ,' for the whole year'
		     scart(3)=1
		  endif
	          lat(i)=lat(n)
	          lon(i)=lon(n)
	          source(i)=source(n)
	          p(i)=p(n)
	          n=n-1
	      else if (abs(lat(i)-42.0394).lt.1.0e-10.and. &
		       abs(lon(i)-13.4425).lt.1.0e-10.and.anno.eq.2003) then
	          if (scart(4).eq.0) then
	             if (logf)write(6,'(10x,a,2f10.5,i5,f10.1,a)') 'Discarded ', &
		            lat(i),lon(i),source(i),p(i) ,' for the whole year'
		     scart(4)=1
		  endif
	          lat(i)=lat(n)
	          lon(i)=lon(n)
	          source(i)=source(n)
	          p(i)=p(n)
	          n=n-1
	      else if (abs(lat(i)-41.9611).lt.1.0e-10.and. &
		       abs(lon(i)-13.6517).lt.1.0e-10.and.anno.eq.2005) then
	          if (scart(5).eq.0) then
	             if (logf)write(6,'(10x,a,2f10.5,i5,f10.1,a)') 'Discarded ', &
		            lat(i),lon(i),source(i),p(i) ,' for the whole year'
		     scart(5)=1
		  endif
	          lat(i)=lat(n)
	          lon(i)=lon(n)
	          source(i)=source(n)
	          p(i)=p(n)
	          n=n-1
	      else if (abs(lat(i)-42.0394).lt.1.0e-10.and. &
		       abs(lon(i)-13.4425).lt.1.0e-10.and.anno.eq.2006) then
	          if (scart(6).eq.0) then
	             if (logf)write(6,'(10x,a,2f10.5,i5,f10.1,a)') 'Discarded ', &
		            lat(i),lon(i),source(i),p(i) ,' for the whole year'
		     scart(6)=1
		  endif
	          lat(i)=lat(n)
	          lon(i)=lon(n)
	          source(i)=source(n)
	          p(i)=p(n)
	          n=n-1
	      else if (abs(lat(i)-41.9611).lt.1.0e-10.and. &
		       abs(lon(i)-13.6517).lt.1.0e-10.and.anno.eq.2006) then
	          if (scart(7).eq.0) then
	             if (logf)write(6,'(10x,a,2f10.5,i5,f10.1,a)') 'Discarded ', &
		            lat(i),lon(i),source(i),p(i) ,' for the whole year'
		     scart(7)=1
		  endif
	          lat(i)=lat(n)
	          lon(i)=lon(n)
	          source(i)=source(n)
	          p(i)=p(n)
	          n=n-1
	      else if (p(i).gt.50.0) then
	          if (logf)write(6,'(10x,a,2f10.5,i5,i10,f10.1)') 'Discarded ', &
		   lat(i),lon(i),source(i),idate,p(i) 
	          lat(i)=lat(n)
	          lon(i)=lon(n)
	          source(i)=source(n)
	          p(i)=p(n)
	          n=n-1
	      endif
!----------------------------------------------------------------------IIRA-----
	   else if (source(i).eq.2) then
	      if (abs(lat(i)-42.4436).lt.1.0e-10.and. &
		  abs(lon(i)-13.9033).lt.1.0e-10.and.anno.eq.2003) then 
	          if (scart(8).eq.0) then
	             if (logf)write(6,'(10x,a,2f10.5,i5,f10.1,a)') 'Discarded ', &
		            lat(i),lon(i),source(i),p(i) ,' for the whole year'
		     scart(8)=1
		  endif
	          lat(i)=lat(n)
	          lon(i)=lon(n)
	          source(i)=source(n)
	          p(i)=p(n)
	          n=n-1
	      else if (abs(lat(i)-42.1208).lt.1.0e-10.and. &
		       abs(lon(i)-15.5072).lt.1.0e-10.and. &
	              (anno.eq.2006.or.anno.eq.2007)) then
	          if (scart(9).eq.0) then
	             if (logf)write(6,'(10x,a,2f10.5,i5,f10.1,a)') 'Discarded ', &
		            lat(i),lon(i),source(i),p(i) ,' for the whole year'
		     scart(9)=1
		  endif
	          lat(i)=lat(n)
	          lon(i)=lon(n)
	          source(i)=source(n)
	          p(i)=p(n)
	          n=n-1
	      else if (abs(lat(i)-42.5822).lt.1.0e-10.and. &
		       abs(lon(i)-14.0556).lt.1.0e-10.and.anno.eq.2006) then
	          if (scart(10).eq.0) then
	             if (logf)write(6,'(10x,a,2f10.5,i5,f10.1,a)') 'Discarded ', &
		            lat(i),lon(i),source(i),p(i) ,' for the whole year'
		     scart(10)=1
		  endif
	          lat(i)=lat(n)
	          lon(i)=lon(n)
	          source(i)=source(n)
	          p(i)=p(n)
	          n=n-1
	      else if (abs(lat(i)-42.3494).lt.1.0e-10.and. &
		       abs(lon(i)-14.1667).lt.1.0e-10.and.anno.eq.2006) then
	          if (scart(11).eq.0) then
	             if (logf)write(6,'(10x,a,2f10.5,i5,f10.1,a)') 'Discarded ', &
		            lat(i),lon(i),source(i),p(i) ,' for the whole year'
		     scart(11)=1
		  endif
	          lat(i)=lat(n)
	          lon(i)=lon(n)
	          source(i)=source(n)
	          p(i)=p(n)
	          n=n-1
	      else if (abs(lat(i)-42.5678).lt.1.0e-10.and. &
		       abs(lon(i)-14.0139).lt.1.0e-10.and. &
		       (anno.eq.2006.or.anno.eq.2007)) then
	          if (scart(12).eq.0) then
	             if (logf)write(6,'(10x,a,2f10.5,i5,f10.1,a)') 'Discarded ', &
		            lat(i),lon(i),source(i),p(i) ,' for the whole year'
		     scart(12)=1
		  endif
	          lat(i)=lat(n)
	          lon(i)=lon(n)
	          source(i)=source(n)
	          p(i)=p(n)
	          n=n-1
	      else if (abs(lat(i)-42.3250).lt.1.0e-10.and. &
		       abs(lon(i)-13.5831).lt.1.0e-10.and.anno.eq.2006) then
	          if (scart(13).eq.0) then
	             if (logf)write(6,'(10x,a,2f10.5,i5,f10.1,a)') 'Discarded ', &
		            lat(i),lon(i),source(i),p(i) ,' for the whole year'
		     scart(13)=1
		  endif
	          lat(i)=lat(n)
	          lon(i)=lon(n)
	          source(i)=source(n)
	          p(i)=p(n)
	          n=n-1
	      else if (abs(lat(i)-41.7422).lt.1.0e-10.and. &
		       abs(lon(i)-14.3831).lt.1.0e-10.and.anno.eq.2006) then
	          if (scart(14).eq.0) then
	             if (logf)write(6,'(10x,a,2f10.5,i5,f10.1,a)') 'Discarded ', &
		            lat(i),lon(i),source(i),p(i) ,' for the whole year'
		     scart(14)=1
		  endif
	          lat(i)=lat(n)
	          lon(i)=lon(n)
	          source(i)=source(n)
	          p(i)=p(n)
	          n=n-1
	      else if (abs(lat(i)-41.9611).lt.1.0e-10.and. &
		       abs(lon(i)-13.6517).lt.1.0e-10.and.anno.eq.2006) then
	          if (scart(15).eq.0) then
	             if (logf)write(6,'(10x,a,2f10.5,i5,f10.1,a)') 'Discarded ', &
		            lat(i),lon(i),source(i),p(i) ,' for the whole year'
		     scart(15)=1
		  endif
	          lat(i)=lat(n)
	          lon(i)=lon(n)
	          source(i)=source(n)
	          p(i)=p(n)
	          n=n-1
	      else if (abs(lat(i)-41.8858).lt.1.0e-10.and. &
		       abs(lon(i)-14.0706).lt.1.0e-10.and. &
		       (anno.eq.2006)) then
	          if (scart(16).eq.0) then
	             if (logf)write(6,'(10x,a,2f10.5,i5,f10.1,a)') 'Discarded ', &
		            lat(i),lon(i),source(i),p(i) ,' for the whole year'
		     scart(16)=1
		  endif
	          lat(i)=lat(n)
	          lon(i)=lon(n)
	          source(i)=source(n)
	          p(i)=p(n)
	          n=n-1
	      else if (p(i).gt.40.0) then
	          if (logf)write(6,'(10x,a,2f10.5,i5,i10,f10.1)') 'Discarded ', &
		   lat(i),lon(i),source(i),idate,p(i) 
	          lat(i)=lat(n)
	          lon(i)=lon(n)
	          source(i)=source(n)
	          p(i)=p(n)
	          n=n-1
	      endif
	   endif
	enddo
	return
	end subroutine validazione

	subroutine arssandb(iflag,anno,irec,lat,lon,source,p,nn)
	integer anno,irec,source(1),nn
	real lat(1),lon(1),p(1)
	character ofile*50 ; data ofile /'pippuzzo'/ ; save ofile
	character nfile*50 ; data nfile /'pippuzzo'/ ; save nfile
	integer orec ; data orec /-1/ ; save orec
	parameter (nsens=60)
	real mat(nsens,366*24+1),xlat(nsens),xlon(nsens) ; save mat,xlat,xlon
	integer lun1,lun2
	character strum*40
	write(nfile,'(a,i4,a)') 'arssa',anno,'.dat'
	if (nfile.ne.ofile) then
!	   call mvsetflags('Log Level',10.0)
	   call getlun(lun1)
	   call openmuseofiles(lun1,nfile,0)
	   call getlun(lun2)
	   call openmuseofiles(lun2,'comuni.arssa',1)
	   if (iflag.eq.2) then			! Skippa i 60 pluviometri
	      do i=1,nsens 
		 read(lun1) n 
	         read(lun2,'(i3)') n
	         read(lun2,'(i10)') n
	      enddo
	   endif
	   mat=-9999
	   do i=1,nsens
	      read(lun1) n,(mat(i,j),j=1,n)
	      read(lun2,'(i3)') n
	      read(lun2,'(i10,2x,a40,2f8.4)') n,strum,xlat(i),xlon(i)
	   enddo
	   close(lun1) ; close(lun2)
!	   call mvsetflags('Log Level',0.0)
	   ofile=nfile
	endif
	nstart=nn
	do i=1,nsens
	   if (mat(i,irec).gt.-1000.0) then
	      nn=nn+1
	      lat(nn)=xlat(i)
	      lon(nn)=xlon(i)
	      p(nn)=mat(i,irec)
	      source(nn)=1
	   endif
	enddo
	if (iflag.eq.2) then
!	   write (6,'(15x,i3,a)') nn-nstart,' temperature from new ARSSA db.'
	else
!	   write (6,'(15x,i3,a)') nn-nstart,' rain gauge from new ARSSA db.'
	endif
	return
	end subroutine arssandb

	integer function isensorchk(isource,icode)
	isensorchk=1

!			Il sensore di ortucchio da i numeri da gennaio 2006
	if (isource.eq.1.and.icode.eq.51) isensorchk=0
	return
	end function isensorchk

	subroutine takel1comp(irec,iday,ihou,l1)
	integer irec,iday,ihou,l1
	if (mod(irec,24).eq.0) then
	   iday=irec/24
	   ihou=24
	else
	   iday=irec/24+1
	   ihou=mod(irec,24)
	endif
	l1=(iday-1)*96+(ihou-1)*4+1
	return
	end subroutine takel1comp

	subroutine dewetra(iflag,anno,irec,lat,lon,source,p,n)
	integer anno,irec,source(1),n
	real lat(1),lon(1),p(1)
	parameter (maxcom=4000)
	real v(maxcom,367*24*4),xlat(maxcom),xlon(maxcom)
	logical cistanno /.true./
	logical ifirst /.true./
	save ifirst,cistanno,ncom,v,xlat,xlon
	character*45 com,sens
	integer code,lun
	maxc=mvlibmagicnum(1)
	if (ifirst) then
	   do i=1,maxc
	      call museoanagrafica('dewetra',i,com,sens,xlat(i),xlon(i))
	   enddo
	   ifirst=.false.
	   call getlun(lun)
	   call mvsetflags('Openmuseofiles Behaviour',1.0)
	   if (iflag.eq.1) then
	      call openmuseodb(lun,'dewrain',anno)
	   else
	      call openmuseodb(lun,'dewtemp',anno)
	   endif
	   call mvsetflags('Openmuseofiles Behaviour',0.0)
	   if (igetmvlibint(1).eq.1) then
	      cistanno=.false.
	      return
	   endif
	   v=-9999
	   do i=1,maxc
	      read(lun) nn,(v(i,j),j=1,nn)
	   enddo
!	   do i=1,maxc ; do j=1,367*24*4	! Perchè lo avevo fatto?
!	      if (v(i,j).gt.20) v(i,j)=-9999
!	   enddo ; enddo
	   close(lun)
	endif
	n=0
	if (.not. cistanno) return
	jj1=(irec-1)*4+1 ; jj2=(irec-1)*4+4
	do ii=1,maxc
	   nn=0 ; vv=0
	   do jj=jj1,jj2
	      if (nint(v(ii,jj)).ne.-9999) then
	         nn=nn+1 ; vv=vv+v(ii,jj)
	      endif
	   enddo
	   if (nn.gt.0) then
	      n=n+1 ; lat(n)=xlat(ii) ; lon(n)=xlon(ii) ; source(n)=5
	      if (iflag.eq.1) then
		 if (vv.lt.0) vv=0
		 p(n)=vv
	      else
		 p(n)=vv/nn
	      endif
	   endif
	enddo
	return
	end subroutine dewetra

	subroutine acqwa(iflag,anno,irec,lat,lon,source,p,n)
	integer anno,irec,source(1),n,nd
	real lat(1),lon(1),p(1)
	parameter (maxcom=260,nrec=501)
	real v(maxcom,367*24),xlat(maxcom),xlon(maxcom)
	character*45 com,sens
	integer code,lun
	logical ifirst /.true./ 
	save ifirst,ncom,v,xlat,xlon,ng
	if (ifirst) then
	   ng=0
	   call getlun(lun) ; call openmuseodb(lun,'acqwapo',anno) ; v=-9999
	   do ii=1,nrec
	      call museoanagrafica('acqwapo',ii,com,sens,xlat(ng+1),xlon(ng+1))
	      if (sens(1:len_trim(sens)).eq.'mm'.and.iflag.eq.1) then
	         ng=ng+1
	         read (lun) nd,(v(ng,i),i=1,nd)
	      else if (sens(1:len_trim(sens)).eq.'oC'.and.iflag.eq.2) then
	         ng=ng+1
	         read (lun) nd,(v(ng,i),i=1,nd)
	      else
		 read (lun) nd
	      endif
	   enddo
	   close(lun)
	   ifirst=.false.
	endif
	do ii=1,ng
	   if (nint(v(ii,irec)).ne.-9999) then
	      n=n+1
	      lat(n)=xlat(ii)
	      lon(n)=xlon(ii)
	      p(n)=v(ii,irec)
	      source(n)=4
	   endif
	enddo
	return
	end subroutine acqwa

	subroutine lazio(iflag,anno,irec,lat,lon,source,p,n)
	integer anno,irec,source(1),n,nd
	real lat(1),lon(1),p(1)
	parameter (maxcom=200)
	real v(maxcom,367*24*4),xlat(maxcom),xlon(maxcom)
	character*35 com
	integer code,lun
	logical ifirst /.true./ 
	save ifirst,ncom,v,xlat,xlon
	if (iflag.eq.2) return
	if (ifirst) then
	   call getlun(lun)
	   call openmuseofiles(lun,'comuni.lazio',1)
	   do i=1,maxcom
	      read(lun,'(a35,2x,i9,2f10.5)',end=10) com,code,xlat(i),xlon(i)
	      ncom=i
	   enddo
  10	   close(lun)
	   call openmuseodb(lun,'lazio',anno)
	   v=-9999
	   do ii=1,ncom
	      read (lun) nd,(v(ii,i),i=1,nd)
	   enddo
	   close(lun)
	   ifirst=.false.
	endif
	call takel1comp(irec,iday,ihou,l1)
	do ii=1,ncom
	   pio=0.0
	   nm=0
	   do i=l1,l1+3
	      if (nint(v(ii,i)).ne.-9999) then
	         nm=nm+1
	         pio=pio+v(ii,i)
	      endif
	   enddo
	   if (nm.gt.0) then
	      n=n+1
	      lat(n)=xlat(ii)
	      lon(n)=xlon(ii)
	      p(n)=pio
	      source(n)=3
	   endif
	enddo
	close(41)
	return
	end subroutine lazio

	real function xmindist(xlat,xlon,lat,lon,n)
	real lat(1),lon(1),xlat,xlon
	xmindist=1.e30
	do i=1,n
	   xx=distance(xlat,xlon,lat(i),lon(i))
	   if (xx.lt.xmindist) xmindist=xx
	enddo
	return
	end function xmindist

	subroutine iira(iflag,anno,irec,lat,lon,source,p,n)
	integer anno,irec,source(1),n,nd
	real lat(1),lon(1),p(1)
	integer ir,icd,good(500)
	real v(500,367*24*4),xlat(500),xlon(500) ; save v,xlat,xlon
	character mis*50,sens*50,loc*50,com*50,prov*2,str*60
	character strm*30
	logical ifirst /.true./ 
	save ifirst,good,ncom
	if (iflag.eq.2) then
	   strm='termometro aria'
	else
	   strm='pluviometro'
	endif
	nstrm=len_trim(strm)
	if (ifirst) then
	   call openmuseofiles(21,'comuni.iira',1)
	   do i=1,500
	      read(21,'(i3,1x,a45,a20,a2)',end=10) ir,loc,com,prov
	      read(21,'(4x,i6,2x,a30,a10,2f8.4)') icd,sens,mis,xlat(i),xlon(i)
	      ncom=i
	      call cv2lower(sens,sens)
	      if (sens(1:len_trim(sens)).eq.strm(1:nstrm)) then
	         good(i)=1
	      else
	         good(i)=0
	      endif
	   enddo
  10	   close(21)
	   write (str,'(a,i4,a)') 'iira',anno,'.dat'
	   call openmuseofiles(31,str,0)
	   v=-9999
	   do ii=1,ncom
	      read (31) nd,(v(ii,i),i=1,nd)
	   enddo
	   close(31)
	   ifirst=.false.
	endif
	nstart=n
	rewind (31)
	if (mod(irec,24).eq.0) then
	   iday=irec/24
	   ihou=24
	else
	   iday=irec/24+1
	   ihou=mod(irec,24)
	endif
	l1=(iday-1)*96+(ihou-1)*4+1
	do ii=1,ncom
	   if (good(ii).eq.1) then
	      pio=0.0
	      nm=0
	      do i=l1,l1+3
	         if (nint(v(ii,i)).ne.-9999) then
	            nm=nm+1
	            pio=pio+v(ii,i)
		 endif
	      enddo
	      if (nm.gt.0) then
		 n=n+1
	         lat(n)=xlat(ii)
	         lon(n)=xlon(ii)
	         if (iflag.eq.2) then
	            p(n)=pio/nm
		 else
	            p(n)=pio
		 endif
	         source(n)=2
	      endif
	   endif
	enddo
	return
        end subroutine iira
!   Questa utility aggiorna il database orario sull'Italia delle precipitazioni
!   e temperature. Per ogni anno vengono collezionati 8760 record in ognuno
!   dei quali ci si scrivono tutte le precipitazioni (temperature) osservate
!   ad ogni intervallo di un'ora. Si usa con
!
!    ${IO} [FLAG] [ANNO] [LOG]
!
!   La flag può essere 01, 02, 11 o 12
!
!   La seconda cifra seleziona le piogge (1) o le temperature (2). La prima cifra
!   seleziona se si vuole utilizzare il solo database dewetra che copr-e(irebbe)
!   tutta l'Italia (primo digit=0), ovvero le varie serie temporali che sono 
!   state collezionate nel corso degli anni: iira, arssa, idrografico 
!   lazio-umbria, acqwa, ecc. (primo digit=1)
!
!   Se un qualunque terzo argomento viene passato si attiva il log. Per creare
!   i db dai dati di acqwa.
!
!   hourlydb 11 2002 a
!   hourlydb 12 2002 a
end module mod_hourdb
