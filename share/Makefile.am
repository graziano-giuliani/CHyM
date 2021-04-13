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
if CHYM_PARALLEL
rcsuffix = MPI
else
rcsuffix = Serial
endif

##vpath %.inc .
##vpath %.o $(CHYM_ROOT)/external:.

noinst_LIBRARIES = librchlib.a
librchlib_a_SOURCES = mod_global.f90 mod_internal.f90 mod_libmv.f90 mod_runparams.f90 mod_param.f90 mod_phys.f90 mod_time.f90 mod_strman.f90 mod_statparams.f90 mod_mssg.f90 mod_stats.f90 mod_vector.f90 mod_html.f90 mod_crtstatic.f90 mod_museo.f90 mod_ncio.f90 mod_cellaut.f90 mod_chymlib.f90 mod_interp.f90 mod_varandtypes.f90 mod_mpimess.f90 mod_mm5.f90

## FCFLAGS += -I. `nf-config --fflags`
FCFLAGS += -I. $(NC_INCLUDES)
FCFLAGS += -I. -I/$(CHYM_ROOT)/external
FCFLAGS += -L. -L/$(CHYM_ROOT)/external -lmpiser

FCLD = $(MPIFC)

%.o: %.f90
	$(FCLD) $(CPPFLAGS) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod *__genmod.*

clean-local:
	rm -f *.mod *__genmod.*

mod_global.o: mod_global.f90
	 $(MPIFC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<
mod_interp.o: mod_interp.f90 mod_libmv.o mod_internal.o mod_vector.o mod_phys.o mod_cellaut.o mod_museo.o
mod_html.o: mod_html.f90 mod_strman.o mod_libmv.o mod_time.o
mod_vector.o: mod_vector.f90
mod_mssg.o: mod_mssg.f90 mod_strman.o mod_html.o
mod_statparams.o: mod_statparams.f90
mod_strman.o: mod_strman.f90
mod_phys.o: mod_phys.f90
mod_runparams.o: mod_runparams.f90
mod_libmv.o: mod_libmv.f90
mod_chymlib.o: mod_chymlib.f90 mod_statparams.o mod_libmv.o mod_runparams.o mod_museo.o mod_phys.o
mod_ncio.o: mod_ncio.f90 mod_internal.o mod_statparams.o mod_libmv.o mod_runparams.o mod_time.o mod_mssg.o mod_museo.o mod_varandtypes.o mod_mpimess.o
	 $(MPIFC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<
mod_cellaut.o: mod_cellaut.f90 mod_libmv.o mod_strman.o mod_vector.o
mod_stats.o: mod_stats.f90 mod_libmv.o mod_vector.o
mod_time.o: mod_time.f90 mod_libmv.o mod_strman.o mod_phys.o mod_internal.o
mod_museo.o: mod_museo.f90 mod_time.o mod_strman.o mod_vector.o mod_internal.o mod_statparams.o
mod_mm5.o: mod_mm5.f90 mod_strman.o mod_museo.o mod_libmv.o
mod_internal.o: mod_internal.f90 mod_libmv.o mod_stats.o mod_strman.o
mod_crtstatic.o: mod_crtstatic.f90 mod_libmv.o mod_internal.o mod_museo.o mod_statparams.o mod_vector.o mod_runparams.o mod_ncio.o mod_strman.o mod_cellaut.o mod_param.o mod_chymlib.o mod_phys.o
mod_param.o: mod_param.f90 mod_runparams.o mod_internal.o mod_phys.o mod_strman.o mod_statparams.o mod_time.o mod_mssg.o mod_libmv.o mod_varandtypes.o mod_mpimess.o
mod_varandtypes.o: mod_varandtypes.f90
mod_mpimess.o: mod_mpimess.f90 mod_varandtypes.o
