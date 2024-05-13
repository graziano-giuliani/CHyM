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
module mod_libmv
  public
  !
  ! PARAMETER definitions
  !
  integer , parameter :: nintpar = 10
  integer , parameter :: nchymp = 60
  integer , parameter :: maxdum = 150000
  integer , parameter :: mm5d = 120
  integer , parameter :: ngks = 2000
  !
  ! COMMON /CELLULA/
  !
  integer , dimension(512) :: irule
  !
  ! COMMON /CHYMCB/
  !
  logical :: chymdom , chymzoom
  integer :: chymlai1 , chymlai2 , chymloi1 , chymloi2 , chymnow
  real :: chymlat1 , chymlat2 , chymlon1 , chymlon2
  integer , dimension(nchymp) :: mchym
  real , dimension(nchymp) :: rchym
  character(150) , dimension(nchymp) :: schym
  !
  ! COMMON /CNVUTM/
  !
  real(8) :: lambda0
  integer , dimension(2) :: utm_grid_zone
  !
  ! COMMON /GKSDUM/
  !
  integer , dimension(ngks) :: ids1 , ids2
  real , dimension(ngks) :: rq , xra , yra
  !
  ! COMMON /LOGUNT/
  !
  integer :: lun00 , lun51 , lun53 , lun54 , lun55 , lun56
  !
  ! COMMON /MM5CB1/
  !
  integer , dimension(1000,20) :: mif
  character(80) , dimension(1000,20) :: mifc , mrfc
  real , dimension(1000,20) :: mrf
  !
  ! COMMON /MM5CB2/
  !
  integer , parameter :: maxfl = 50
  integer , parameter :: mxdf = 100
  character(80) , dimension(mxdf) :: ch1 , ch2
  integer , dimension(mxdf) :: ic1 , ic2 , ic3 , ic4 , id1 , ir1 , ir2 ,  &
                               ir3 , ir4
  integer , dimension(4) :: nd
  real , dimension(mxdf) :: rd1
  !
  ! COMMON /MM5CB4/
  !
  integer , dimension(100) :: lindex
  !
  ! COMMON /MM5V34/
  !
  integer , dimension(50,20) :: bhi
  character(80) , dimension(50,20) :: bhic
  real , dimension(20,20) :: bhr
  character(80) , dimension(20,20) :: bhrc
  character(24) :: current_date
  character(46) :: description
  integer , dimension(4) :: end_index , start_index
  integer :: mm5flag , mm5ndim
  character(len=9) :: mm5id , mm5_name
  character(len=4) :: ordering , staggering
  character(len=25) :: mm5_unit
  real :: xtime
  !
  ! COMMON /MM5ZOOM/
  !
  integer :: latreal , lonreal , latz1 , latz2 , lonrea , lonz1 , lonz2
  logical :: mm5zooming
  !
  ! COMMON /MVFLAG/
  !
  integer , dimension(100) :: iflg
  integer , dimension(nintpar) :: isub
  real , dimension(100) :: rflg
  real , dimension(nintpar) :: rsub
  !
  ! COMMON /MVGRFL/
  !
  integer :: i53 , i54 , i55 , i56 , iopen , lmeth , nhistplot , nplotta ,&
             nplottb
  real , dimension(100) :: rlevels
  !
  ! COMMON /MVINTER/
  !
  integer :: iplatf
  !
  ! COMMON /MVPATH/
  !
  character(len=80) :: filename , pathb , pathm , filenamerst
  integer :: npathb , npathm
  !
  ! COMMON /MVTIME/
  !
  integer , dimension(13) :: dfc
  integer :: di33 , oldyear
  character(9) , dimension(12) :: emonths , months
  character(11) , dimension(7) :: eweekd , weekd
  integer , dimension(12) :: mesi
  !
  ! COMMON /NEWDUM/
  !
  integer , dimension(maxdum) :: i1dum
  real , dimension(mm5d,mm5d) :: mm51 , mm52
  real , dimension(maxdum) :: x1dum , x2dum , y1dum , y2dum
  !
  ! COMMON /PINOUFLE/
  !
  integer , dimension(100) :: icolori
  integer :: numero
  real , dimension(100) :: valori
  !
  ! COMMON /HANDLETIME/
  !
  integer , public :: time,hour,day,month,year

  integer , private :: i

  data iflg(1)/0/            ! log level
  data iflg(2)/700/          ! X size grafici
  data iflg(3)/1/            ! boxplot should or should not draw the frame
  data iflg(4)/0/            ! Plot 1D
  data iflg(5)/0/            ! Dominio 0-3
  data iflg(6)/1/            ! Mese Corrente
  data iflg(7)/1/            ! Check Landuse
  data iflg(8)/0/            ! Logo sui plot
  data iflg(9)/0/            ! Stile HTML
  data iflg(10)/2/           ! Stile con cui produce le date
  data iflg(11)/-1/          ! confini politici (no longer used)
  data iflg(12)/1/           ! colore delle scritte dei grafici
  data iflg(13)/1/           ! colore degli assi nei plot 1-d
  data iflg(14)/4/           ! colore delle freccette del vento
  data iflg(15)/0/           ! Se > 0 nella sub. meteogrammi setta l'asse x
  data iflg(16)/0/           ! colore dei confini politici
  data iflg(17)/-1/          ! se >= 0 disegna i confini di default gks ncar
  data iflg(18)/1/           ! unita' di misura u-v (1=m/s 2=Km/h)
  data iflg(19)/1/           ! colore delle isolinee (routine mm5contorni)
  data iflg(20)/0/           ! Numero di curve di livello
  data iflg(21)/2/           ! Dimensioni delle scritte (1,2,3 della wtstr)
  data iflg(22)/1/           ! coord.grafiche (1=x,y 2=lat-lon 3=correnti)
  data iflg(23)/0/           ! primo indice delle matrici (0=lon altro=lat)
  data iflg(24)/0/           ! step tra i livelli degli automi cellulari
  data iflg(25)/1/           ! colore sfondo dei grafici: 0=nero altro=bianco
  data iflg(26)/18/          ! numero di colori
  data iflg(27)/1/           ! palette di colori
  data iflg(28)/0/           ! contorno ai pallocchi
  data iflg(29)/1/           ! La label bar di mm5colormap
  data iflg(30)/0/           ! se = 0 viacolvento cancella i singoli frames
  data iflg(31)/0/           ! label asse x ( 1 = ore 2 = mesi 3 = nulla )
  data iflg(32)/0/           ! Numero di colori nei plot 2d
  data iflg(33)/1/           ! Ogni quanti punti griglia fa le freccette
  data iflg(34)/875949887/   ! Random seed for gauss and flat distribution
  data iflg(35)/0/           ! Qnorm - number of bin to be used
  data iflg(36)/0/           ! label asse y ( 3 = nulla )
  data iflg(37)/1/           ! Primo colore per disegnare le palette
  data iflg(38)/1/           ! Ogni quanti colori per disegnare le palette
  data iflg(39)/1/           ! Scrive il massimo e il minimo - subr. contorni
  data iflg(40)/-1/          ! Colore Plotta
  data iflg(41)/-13/         ! Tipo di carattere vedi man gstxfp
  data iflg(42)/20/          ! Colore delle freccette (subr. freccetelle)
  data iflg(43)/0/           ! Colore sfondo dei plot 1-d
  data iflg(44)/0/           ! Used in plotta package - internal only
  data iflg(45)/1/           ! boundary routine also write main town if > 0
  data iflg(46)/3/           ! Visulaizzaplot format (0=bmp, 1=gif, 2=tiff)
  data iflg(47)/620/         ! Size delle finestre negli script html
  data iflg(48)/0/           ! nei plot mm5 se = 1 disegna l'italia
                             !              se = 2 disegna le regioni
                             !              se = 3 disegna le provincie
  data iflg(49)/0/           ! Se = 1 disegna regioni Se = 2 provincie
  data iflg(50)/0/           ! Internal use
  data iflg(51)/0/           ! Se = 1 label orizzont  Se = 2 label vertic
  data iflg(52)/0/           ! Stile delle freccette (routine boxarrow)
  data iflg(53)/-1/          ! Histogram color (routine mvbookplot)
  data iflg(54)/20/          ! Numero di colori della palette corrente
  data iflg(55)/0/           ! Internal use only - plottaframe, plottavw
  data iflg(56)/-1/          ! Size delle label - plotpalette, boxlegend
  data iflg(57)/-1/          ! secolo - subroutine datafrommm5index ecc.
  data iflg(58)/1/           ! controllo video. <>0 = ON
  data iflg(59)/0/           ! Visulaizzaplot behaviour 0=vis ; 1=salva
  data iflg(60)/1/           ! Used by somepoch to define neural initializ.
  data iflg(61)/1/           ! Used by somepoch to define neural distance
  data iflg(62)/0/           ! Used by openmuseofile
  data iflg(63)/0/           ! Used by openmuseofile (internal)
  data iflg(64)/0/           ! Used by boxplot
  data iflg(65)/0/           ! Used by chymdownscaling
  data iflg(66)/0/           ! Background color of html table
  data iflg(67)/3/           ! Histogram frame (subroutine mvbookplot)
  data iflg(68)/1/           ! Histogram statistic (subroutine mvbookplot)
  data iflg(69)/0/           ! Plotta fill
  data iflg(70)/6/           ! Logical unit used by CHyM Model for log
  data iflg(71)/0/           ! CHyM Standard domain
  data iflg(72)/0/           ! Calendar
  data (rflg(i),i=1,4)/0.05 , 0.95 , 0.05 , 0.95/   ! limite grafici
  data (rflg(i),i=5,7)/0.0 , 0.0 , 0.0/             ! curve di livello
  data rflg(8)/1.0/                                 ! size delle linee
  data (rflg(i),i=9,16)/8*1.0/                      ! cellular weight
  data rflg(17)/1.0/            ! Size del carattere, vedi man gschxp
  data rflg(18)/0.5/            !
  data rflg(19)/0.5/            ! rflg(18-20) are defined in jnplot and used in
  data rflg(20)/1.0/            ! jnlegend
  data (rflg(i),i=21,24)/4*-100.0/  ! Lon-Lat boundaries of current plot.
 
  data isub/nintpar* - 9999/
  data rsub/nintpar* - 9999.0/
 
  data mesi/31 , 28 , 31 , 30 , 31 , 30 , 31 , 31 , 30 , 31 , 30 , 31/
  data weekd/'Lunedi'' ' , 'Martedi'' ' , 'Mercoledi'' ' , 'Giovedi'' ' , &
             'Venerdi'' ' , 'Sabato' , 'Domenica'/
  data months/'Gennaio' , 'Febbraio' , 'Marzo' , 'Aprile' , 'Maggio' ,   &
              'Giugno' , 'Luglio' , 'Agosto' , 'Settembre' , 'Ottobre' , &
              'Novembre' , 'Dicembre'/
  data eweekd/'Monday' , 'Tuesday' , 'Wednesday' , 'Thursday' , 'Friday' ,&
              'Saturday' , 'Sunday'/
  data emonths/'January' , 'February' , 'March' , 'April' , 'May' ,   &
               'June' , 'July' , 'August' , 'September' , 'October' , &
               'November' , 'December'/
 
  data mm5zooming/.false./
  data chymzoom/.false./           ! Zoom over the chym plot?
  data chymdom/.false./            ! Chym domain has been defined?
  data mchym(8)/0/                 ! To avoid problems inside chymsavedfield
  data chymnow/-1/                 ! Last time slide read
  data nplotta , nplottb/0 , 0/
 
  data utm_grid_zone/-1 , -1/      ! Conversion utm to lat-lon e viceversa
  data lmeth/0/                    ! Used inside mm5colormap routines
 
  data nhistplot/0/                ! Used inside mvbookplot

  integer :: chymcrec, hourstep
  data chymcrec/-1/
 
end module mod_libmv
