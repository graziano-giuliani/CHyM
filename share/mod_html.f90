!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP CHyM.
!
!    ICTP CHyM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP CHyM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP CHyM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
module mod_html

  use mod_libmv
  use mod_strman
  use mod_time

  character(len=20) , dimension(31) :: colors
  data colors/'cyan' , 'darkcyan' , 'cornflowerblue' , 'crimson' ,     &
    'darkgray' , 'coral' , 'darkorange' , 'darkturquoise' ,            &
    'deeppink' , 'deepskyblue' , 'forestgreen' , 'fuchsia' , 'gold' ,  &
    'green' , 'lightblue' , 'lightgreen' , 'lightseagreen' , 'lime' ,  &
    'magenta' , 'maroon' , 'navy' , 'olive' , 'orangered' ,            &
    'palevioletred' , 'purple' , 'red' , 'saddlebrown' , 'skyblue' ,   &
    'springgreen' , 'tomato' , 'yellow'/

  contains

  subroutine htmlheader(titolo,mese,lun)
    implicit none
    integer :: lun
    character(len=*) :: mese , titolo
    intent (in) lun
    write (lun,'(a)') '<HTML><HEAD><TITLE>'//trim(titolo)
    write (lun,'(a)') '</TITLE>'
    if ( len_trim(mese)>0 ) then
      write (lun,'(a)') '<SCRIPT>'
      write (lun,'(a)') 'function '//trim(mese)//'(str) {'
      write (lun,'(a,i3,a,i3,a)') 'searchWin = window.open(str'//','''//    &
              trim(mese)//''',''scrollbars=no'//                  &
              ',resizable=no,width=' , iflg(47) ,',height=' , iflg(47) ,    &
              ',status=no,location=no,toolbar=no'');'
      write (lun,'(a)') '//        searchWin.refer = self;'
      write (lun,'(/,a,/)') '}'
      write (lun,'(/,a,/)') '</SCRIPT>'
    end if
    write (lun,'(a)') '</HEAD>'
    if ( iflg(9)<1 .or. iflg(9)>3 ) iflg(9) = 1
    if ( iflg(9)==1 ) then
      write (lun,'(a)') '<BODY bgcolor="white" ><BR><CENTER>'
      write (lun,'(a)') '<FONT face="Lucida" color="red" size="5">'
      write (lun,'(a)') titolo
      write (lun,'(a)') '</FONT></CENTER><P>'
      write (lun,'(a)') '<font size="1" face="Verdana, Arial, '//           &
                        'Helvetica, sans-serif">'
    else if ( iflg(9)==2 ) then
      write (lun,'(a)') '<body bgcolor="white"><CENTER>'
    else if ( iflg(9)==3 ) then
      write (lun,'(a)') '<BODY TEXT="YELLOW" BGCOLOR="BLACK" '//            &
                        'VLINK="00ff00" LINK="00ff00"> <CENTER> '
    end if
  end subroutine htmlheader
 
  subroutine javascript1(link,scritta,window,icon,lun)
    implicit none
    character(len=*) :: icon , link , scritta , window
    integer :: lun
    intent (in) lun , scritta
    integer :: ii
    character(len=40) :: ltmp
    character(len=60) :: tbcl
    call mvgetiflags(66,ii)
    if ( ii>=1 .and. ii<=31 ) then
      tbcl = 'BGCOLOR="'//trim(colors(ii))//'"'
    else
      tbcl = ' '
    end if
    ltmp = scritta
    call strsub('''',' ',ltmp)
    call strsub('''',' ',ltmp)
    if ( len_trim(icon)>0 ) then
      write (lun,'(a)') '<TD ALIGN="CENTER"'//trim(tbcl)   &
                //'> <a href="javascript:'//trim(window)   &
                //'('''//trim(link)                        &
                //''');"onmouseover="window.status='''//   &
                trim(ltmp)//''';return true"><IMG SRC="'// &
                trim(icon)//'" alt="'//                    &
                trim(ltmp)//'"></a><BR>'//                 &
                '<font size="2" face="Verdana, Arial, '//  &
                'Helvetica, sans-serif">'//trim(ltmp)      &
                //' </TD>'
    else
      write (lun,'(a)') '<TD ALIGN="CENTER"> '//                &
             '<font size="1" face="Verdana, Arial, Helvetica, sans-serif">'
      write (lun,'(a)') '<a href="javascript:'//trim(window)    &
             //'('''//trim(link)                                &
             //''');"onmouseover="window.status='''//           &
             trim(ltmp)//''';return true">'//trim(ltmp)//'</a></TD>'
    end if
  end subroutine javascript1
     
  subroutine htmlsignature(programma,versione,autore,lun)
    implicit none
    character(len=*) :: autore , programma , versione
    integer :: lun
    intent (in) lun
    integer :: anno , giorno , mese , minuti , ora
    character(len=60) :: date
    call whattimeisit(minuti,ora,giorno,mese,anno)
    call mvsetflags('Data Style',2.0)
    call datafrommin(minuti,ora,giorno,mese,anno,date)
    write (lun,'(a)') '<CENTER> <TABLE BORDER="0" WIDTH="90%">'//  &
                     '<TR><TD ALIGN="LEFT">'//                     &
             '<font size="1" face="Verdana, Arial, Helvetica, sans-serif">'
    write (lun,'(a)') '<B>HTML code automatically generated using'
    write (lun,'(a)') '"'//trim(programma)//'" software '
    if ( len_trim(versione)>0 ) write (lun,'(a)') 'version '//versione
    if ( len_trim(autore)>0 ) write (lun,'(a)') 'developped by '//autore
    write (lun,'(a)') '</B></TD><TD ALIGN="RIGHT">'//                       &
             '<font size="1" face="Verdana, Arial, Helvetica, sans-serif">' &
            //'<B>This document has been updated on '//date//               &
             '</B></TD></TR></TABLE>'
    write (lun,'(a)') '</CENTER></FONT>'
  end subroutine htmlsignature
     
  subroutine htmltable(ele,nele,nxrow,startl,endl,iunit)
    implicit none
    character(len=*) :: endl , startl
    integer :: iunit , nele , nxrow
    character(len=*) , dimension(nele) :: ele
    intent (in) nele , nxrow
    integer :: i , nn
    character(len=132) :: str
    write (iunit,'(a)') '<TABLE border=0 align=center ><TR>'
    str = ' '
    nn = len_trim(startl)
    do i = 1 , nele
      if ( nn>0 ) call makefilename(startl,i,endl,str)
      call htmltabel(ele(i),str,0,iunit)
      if ( mod(i,nxrow)==0 ) write (iunit,'(a)') '</TR><TR>'
    end do
    write (iunit,'(a)') '</TR></TABLE>'
  end subroutine htmltable
     
  subroutine htmltabel(element,link,align,lunit)
    implicit none
    integer :: align , lunit
    character(len=*) :: element , link
    intent (in) align , element , lunit
    integer :: nn
    character(len=60) :: tbcl
    if ( iflg(66)>=1 .and. iflg(66)<=31 ) then
      tbcl = 'BGCOLOR="'//trim(colors(iflg(66)))//'"'
    else
      tbcl = ' '
    end if
    if ( align<0 ) then
      write (lunit,'(a)') '<TD ALIGN="LEFT" '//trim(tbcl)//'>'
    else if ( align>0 ) then
      write (lunit,'(a)') '<TD ALIGN="RIGHT" '//trim(tbcl)//'>'
    else
      write (lunit,'(a)') '<TD ALIGN="CENTER" '//trim(tbcl)//'>'
    end if
    nn = len_trim(link)
    write (lunit,'(a)') '<font size="2" face="Verdana, Arial, Helvetica, '//&
                       'sans-serif">'
    if ( nn>0 ) write (lunit,'(a)') '<A HREF="'//link(1:nn)//'">'
    write (lunit,'(a)') element
    if ( nn>0 ) write (lunit,'(a)') '</A>'
  end subroutine htmltabel
     
  subroutine mkhtmlpres(files,htmlname,title)
    implicit none
    character(len=*) :: files , htmlname , title
    intent (in) files
    character(len=132) :: cfile , html , lsfile
    integer :: i , n , nn
    call system('ls -1 '//files//' >/tmp/elencofiles.dat')
    open (10,file='/tmp/elencofiles.dat',status='old')
    n = 0
    do while ( n>=0 .or. n<=0 )
      read (10,'(a)',end=100) cfile
      lsfile = cfile
      n = n + 1
      call makefilename(htmlname,n,'.html',html)
      open (20,file=html,status='unknown')
      call htmlheader(title,' ',20)
      write (20,'(a)') '<CENTER>'
      if ( n>1 ) then
        call makefilename(htmlname,n-1,'.html',html)
        nn = 1
        do i = 1 , len_trim(html)
          if ( html(i:i)=='/' ) nn = i + 1
        end do
        write (20,'(a)') '<A HREF="'//trim(html)//'">'
        write (20,'(a,i10,a)') 'Slice # ' , n - 1 , '</A>'
      else
        write (6,'(7x,a)') 'MvLib: mkhtmlpres, first file is: '//trim(html)
      end if
      write (20,'(a,i10)') 'Slice # ' , n
      call makefilename(htmlname,n+1,'.html',html)
      nn = 1
      do i = 1 , len_trim(html)
        if ( html(i:i)=='/' ) nn = i + 1
      end do
      write (20,'(a)') '<A HREF="'//trim(html)//'">'
      write (20,'(a,i10,a)') 'Slice # ' , n + 1 , '</A><BR>'
      nn = 1
      do i = 1 , len_trim(cfile)
        if ( cfile(i:i)=='/' ) nn = i + 1
      end do
      write (20,'(a)') '<IMG SRC="'//trim(cfile)//'"><BR>'
      close (20)
    end do
  100  close (10)
    call system('/bin/rm /tmp/elencofiles.dat')
    call makefilename(htmlname,n,'.html',html)
    open (20,file=html,status='unknown')
    call htmlheader(title,' ',20)
    write (20,'(a)') '<CENTER>'
    call makefilename(htmlname,n-1,'.html',html)
    nn = 1
    do i = 1 , len_trim(html)
      if ( html(i:i)=='/' ) nn = i + 1
    end do
    write (20,'(a)') '<A HREF="'//trim(html)//'">'
    write (20,'(a,i10,a)') 'Slice # ' , n - 1 , '</A>'
    write (20,'(a,i10,a)') 'Slice # ' , n , '<BR>'
    nn = 1
    do i = 1 , len_trim(lsfile)
      if ( lsfile(i:i)=='/' ) nn = i + 1
    end do
    write (20,'(a)') '<IMG SRC="'//trim(lsfile)//'"><BR>'
    close (20)
  end subroutine mkhtmlpres

end module mod_html
