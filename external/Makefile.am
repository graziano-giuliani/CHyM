#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#    This file is part of ICTP CHyM.
#
#    ICTP CHyM is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    ICTP CHyM is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with ICTP CHyM.  If not, see <http://www.gnu.org/licenses/>.
#
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
include $(CHYM_ROOT)/makeinc

vpath %.F90 $(CHYM_ROOT)/external/mpi-serial
vpath %.F90 $(CHYM_ROOT)/external/mpi-serial
vpath %.c $(CHYM_ROOT)/external/mpi-serial
vpath %.h $(CHYM_ROOT)/external/mpi-serial


if CHYM_PARALLEL
noinst_LIBRARIES =
else
noinst_LIBRARIES = libmpiser.a
endif
libmpiser_a_SOURCES = mpi.c send.c recv.c collective.c req.c list.c \
                      handles.c comm.c group.c time.c pack.c fort.F90 \
                      mpistub.F90

MPIFH= mpif.$(FORT_SIZE).h

mpif.h:
	cp -f $(CHYM_ROOT)/external/mpi-serial/$(MPIFH) mpif.h

%.o: %.c
	$(CC) $(CFLAGS) $(MPEUCFLAGS) -c $<

%.o: %.F90
	  $(FC) $(CPPFLAGS) $(FCFLAGS) $(MPEUFCFLAGS) -c $<

%.o: %.F90
	  $(FC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) $(MPEUFCFLAGS) -c $<

distclean-local:
	rm -f *.mod mpif.h

clean-local:
	rm -f *.mod

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

mpi.o: mpif.h
fort.o: fort.F90
mpistub.o: mpistub.F90 
