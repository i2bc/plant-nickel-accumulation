# Makefile.in generated by automake 1.11 from Makefile.am.
# src/shcl/Makefile.  Generated from Makefile.in by configure.

# Copyright (C) 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002,
# 2003, 2004, 2005, 2006, 2007, 2008, 2009  Free Software Foundation,
# Inc.
# This Makefile.in is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.





pkgdatadir = $(datadir)/mcl
pkgincludedir = $(includedir)/mcl
pkglibdir = $(libdir)/mcl
pkglibexecdir = $(libexecdir)/mcl
am__cd = CDPATH="$${ZSH_VERSION+.}$(PATH_SEPARATOR)" && cd
install_sh_DATA = $(install_sh) -c -m 644
install_sh_PROGRAM = $(install_sh) -c
install_sh_SCRIPT = $(install_sh) -c
INSTALL_HEADER = $(INSTALL_DATA)
transform = $(program_transform_name)
NORMAL_INSTALL = :
PRE_INSTALL = :
POST_INSTALL = :
NORMAL_UNINSTALL = :
PRE_UNINSTALL = :
POST_UNINSTALL = :
build_triplet = x86_64-unknown-linux-gnu
host_triplet = x86_64-unknown-linux-gnu
DIST_COMMON = $(srcdir)/Makefile.am $(srcdir)/Makefile.in \
	$(top_srcdir)/include/include.am
bin_PROGRAMS = clmformat$(EXEEXT) clm$(EXEEXT)
noinst_PROGRAMS = clmdag$(EXEEXT) clmframe$(EXEEXT)
subdir = src/shcl
ACLOCAL_M4 = $(top_srcdir)/aclocal.m4
am__aclocal_m4_deps = $(top_srcdir)/acinclude.m4 \
	$(top_srcdir)/configure.ac
am__configure_deps = $(am__aclocal_m4_deps) $(CONFIGURE_DEPENDENCIES) \
	$(ACLOCAL_M4)
mkinstalldirs = $(install_sh) -d
CONFIG_HEADER = $(top_builddir)/config.h
CONFIG_CLEAN_FILES =
CONFIG_CLEAN_VPATH_FILES =
am__installdirs = "$(DESTDIR)$(bindir)" "$(DESTDIR)$(bindir)"
PROGRAMS = $(bin_PROGRAMS) $(noinst_PROGRAMS)
am_clm_OBJECTS = clm.$(OBJEXT) clminfo.$(OBJEXT) clmimac.$(OBJEXT) \
	clmadjust.$(OBJEXT) clmresidue.$(OBJEXT) clmclose.$(OBJEXT) \
	clmmate.$(OBJEXT) clmmeet.$(OBJEXT) clmps.$(OBJEXT) \
	clmdist.$(OBJEXT) clmorder.$(OBJEXT) report.$(OBJEXT)
clm_OBJECTS = $(am_clm_OBJECTS)
clm_LDADD = $(LDADD)
clm_DEPENDENCIES = ../mcl/libmcl.a ../clew/libclew.a \
	../gryphon/libgryphon.a ../impala/libimpala.a \
	../../util/libutil.a
am_clmdag_OBJECTS = clmdag.$(OBJEXT) report.$(OBJEXT)
clmdag_OBJECTS = $(am_clmdag_OBJECTS)
clmdag_LDADD = $(LDADD)
clmdag_DEPENDENCIES = ../mcl/libmcl.a ../clew/libclew.a \
	../gryphon/libgryphon.a ../impala/libimpala.a \
	../../util/libutil.a
am_clmformat_OBJECTS = clmformat.$(OBJEXT) report.$(OBJEXT)
clmformat_OBJECTS = $(am_clmformat_OBJECTS)
clmformat_LDADD = $(LDADD)
clmformat_DEPENDENCIES = ../mcl/libmcl.a ../clew/libclew.a \
	../gryphon/libgryphon.a ../impala/libimpala.a \
	../../util/libutil.a
am_clmframe_OBJECTS = clmframe.$(OBJEXT) report.$(OBJEXT)
clmframe_OBJECTS = $(am_clmframe_OBJECTS)
clmframe_LDADD = $(LDADD)
clmframe_DEPENDENCIES = ../mcl/libmcl.a ../clew/libclew.a \
	../gryphon/libgryphon.a ../impala/libimpala.a \
	../../util/libutil.a
am__vpath_adj_setup = srcdirstrip=`echo "$(srcdir)" | sed 's|.|.|g'`;
am__vpath_adj = case $$p in \
    $(srcdir)/*) f=`echo "$$p" | sed "s|^$$srcdirstrip/||"`;; \
    *) f=$$p;; \
  esac;
am__strip_dir = f=`echo $$p | sed -e 's|^.*/||'`;
am__install_max = 40
am__nobase_strip_setup = \
  srcdirstrip=`echo "$(srcdir)" | sed 's/[].[^$$\\*|]/\\\\&/g'`
am__nobase_strip = \
  for p in $$list; do echo "$$p"; done | sed -e "s|$$srcdirstrip/||"
am__nobase_list = $(am__nobase_strip_setup); \
  for p in $$list; do echo "$$p $$p"; done | \
  sed "s| $$srcdirstrip/| |;"' / .*\//!s/ .*/ ./; s,\( .*\)/[^/]*$$,\1,' | \
  $(AWK) 'BEGIN { files["."] = "" } { files[$$2] = files[$$2] " " $$1; \
    if (++n[$$2] == $(am__install_max)) \
      { print $$2, files[$$2]; n[$$2] = 0; files[$$2] = "" } } \
    END { for (dir in files) print dir, files[dir] }'
am__base_list = \
  sed '$$!N;$$!N;$$!N;$$!N;$$!N;$$!N;$$!N;s/\n/ /g' | \
  sed '$$!N;$$!N;$$!N;$$!N;s/\n/ /g'
SCRIPTS = $(bin_SCRIPTS) $(noinst_SCRIPTS)
DEFAULT_INCLUDES = -I. -I$(top_builddir)
depcomp = $(SHELL) $(top_srcdir)/autofoo/depcomp
am__depfiles_maybe = depfiles
am__mv = mv -f
COMPILE = $(CC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) \
	$(CPPFLAGS) $(AM_CFLAGS) $(CFLAGS)
CCLD = $(CC)
LINK = $(CCLD) $(AM_CFLAGS) $(CFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@
SOURCES = $(clm_SOURCES) $(clmdag_SOURCES) $(clmformat_SOURCES) \
	$(clmframe_SOURCES)
DIST_SOURCES = $(clm_SOURCES) $(clmdag_SOURCES) $(clmformat_SOURCES) \
	$(clmframe_SOURCES)
ETAGS = etags
CTAGS = ctags
DISTFILES = $(DIST_COMMON) $(DIST_SOURCES) $(TEXINFOS) $(EXTRA_DIST)
ACLOCAL = ${SHELL} /home/cpereira/Phylogeny/mcl-12-068/autofoo/missing --run aclocal-1.11
AMTAR = ${SHELL} /home/cpereira/Phylogeny/mcl-12-068/autofoo/missing --run tar
AUTOCONF = ${SHELL} /home/cpereira/Phylogeny/mcl-12-068/autofoo/missing --run autoconf
AUTOHEADER = ${SHELL} /home/cpereira/Phylogeny/mcl-12-068/autofoo/missing --run autoheader
AUTOMAKE = ${SHELL} /home/cpereira/Phylogeny/mcl-12-068/autofoo/missing --run automake-1.11
AWK = gawk
CC = gcc
CCDEPMODE = depmode=gcc3
CFLAGS = -g -O2 -pthread
CPP = gcc -E
CPPFLAGS = 
CYGPATH_W = echo
DEFS = -DHAVE_CONFIG_H
DEPDIR = .deps
ECHO_C = 
ECHO_N = -n
ECHO_T = 
EGREP = /bin/grep -E
EXEEXT = 
GREP = /bin/grep
INSTALL = /usr/bin/install -c
INSTALL_DATA = ${INSTALL} -m 644
INSTALL_PROGRAM = ${INSTALL}
INSTALL_SCRIPT = ${INSTALL}
INSTALL_STRIP_PROGRAM = $(install_sh) -c -s
LDFLAGS = 
LIBOBJS = 
LIBS =  -lm 
LN_S = ln -s
LTLIBOBJS = 
MAINT = #
MAKEINFO = ${SHELL} /home/cpereira/Phylogeny/mcl-12-068/autofoo/missing --run makeinfo
MKDIR_P = /bin/mkdir -p
OBJEXT = o
PACKAGE = mcl
PACKAGE_BUGREPORT = mcl-devel@lists.micans.org
PACKAGE_NAME = mcl
PACKAGE_STRING = mcl 12-068
PACKAGE_TARNAME = mcl
PACKAGE_URL = 
PACKAGE_VERSION = 12-068
PATH_SEPARATOR = :
PTHREAD_CC = gcc
PTHREAD_CFLAGS = -pthread
PTHREAD_LIBS = 
RANLIB = ranlib
SET_MAKE = 
SHELL = /bin/sh
STRIP = 
VERSION = 12-068
abs_builddir = /home/cpereira/Phylogeny/mcl-12-068/src/shcl
abs_srcdir = /home/cpereira/Phylogeny/mcl-12-068/src/shcl
abs_top_builddir = /home/cpereira/Phylogeny/mcl-12-068
abs_top_srcdir = /home/cpereira/Phylogeny/mcl-12-068
ac_ct_CC = gcc
acx_pthread_config = 
am__include = include
am__leading_dot = .
am__quote = 
am__tar = ${AMTAR} chof - "$$tardir"
am__untar = ${AMTAR} xf -
bindir = ${exec_prefix}/bin
build = x86_64-unknown-linux-gnu
build_alias = 
build_cpu = x86_64
build_os = linux-gnu
build_vendor = unknown
builddir = .
datadir = ${datarootdir}
datarootdir = ${prefix}/share
docdir = ${datarootdir}/doc/${PACKAGE_TARNAME}
dvidir = ${docdir}
exec_prefix = ${prefix}
host = x86_64-unknown-linux-gnu
host_alias = 
host_cpu = x86_64
host_os = linux-gnu
host_vendor = unknown
htmldir = ${docdir}
includedir = ${prefix}/include
infodir = ${datarootdir}/info
install_sh = ${SHELL} /home/cpereira/Phylogeny/mcl-12-068/autofoo/install-sh
libdir = ${exec_prefix}/lib
libexecdir = ${exec_prefix}/libexec
localedir = ${datarootdir}/locale
localstatedir = ${prefix}/var
mandir = ${datarootdir}/man
mkdir_p = /bin/mkdir -p
oldincludedir = /usr/include
pdfdir = ${docdir}
prefix = /home/cpereira/Phylogeny
program_transform_name = s,x,x,
psdir = ${docdir}
sbindir = ${exec_prefix}/sbin
sharedstatedir = ${prefix}/com
srcdir = .
sysconfdir = ${prefix}/etc
target_alias = 
top_build_prefix = ../../
top_builddir = ../..
top_srcdir = ../..
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)
LDADD = ../mcl/libmcl.a ../clew/libclew.a  ../gryphon/libgryphon.a ../impala/libimpala.a ../../util/libutil.a
clmformat_SOURCES = clmformat.c report.c
clmframe_SOURCES = clmframe.c report.c
clmdag_SOURCES = clmdag.c report.c
clm_SOURCES = clm.c clminfo.c clmimac.c clmadjust.c clmresidue.c clmclose.c clmmate.c \
					clmmeet.c clmps.c clmdist.c clmorder.c report.c

bin_SCRIPTS = mclpipeline
noinst_SCRIPTS = mclgraga mclgrep
EXTRA_DIST = clmorder.h clmdist.h clminfo.h clmclose.h clmmate.h clmmeet.h clm.h \
					clmadjust.h clmresidue.h clmimac.h clmps.h clmps.defs.h\
					report.h report.c \
					$(bin_SCRIPTS) \
               mclgraga mclgrep clmformat.zmm clmformat.zmm.h \
					clmps.defs clmps.defs.h

all: all-am

.SUFFIXES:
.SUFFIXES: .c .o .obj
$(srcdir)/Makefile.in: # $(srcdir)/Makefile.am $(top_srcdir)/include/include.am $(am__configure_deps)
	@for dep in $?; do \
	  case '$(am__configure_deps)' in \
	    *$$dep*) \
	      ( cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh ) \
	        && { if test -f $@; then exit 0; else break; fi; }; \
	      exit 1;; \
	  esac; \
	done; \
	echo ' cd $(top_srcdir) && $(AUTOMAKE) --gnu src/shcl/Makefile'; \
	$(am__cd) $(top_srcdir) && \
	  $(AUTOMAKE) --gnu src/shcl/Makefile
.PRECIOUS: Makefile
Makefile: $(srcdir)/Makefile.in $(top_builddir)/config.status
	@case '$?' in \
	  *config.status*) \
	    cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh;; \
	  *) \
	    echo ' cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@ $(am__depfiles_maybe)'; \
	    cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@ $(am__depfiles_maybe);; \
	esac;

$(top_builddir)/config.status: $(top_srcdir)/configure $(CONFIG_STATUS_DEPENDENCIES)
	cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh

$(top_srcdir)/configure: # $(am__configure_deps)
	cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh
$(ACLOCAL_M4): # $(am__aclocal_m4_deps)
	cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh
$(am__aclocal_m4_deps):
install-binPROGRAMS: $(bin_PROGRAMS)
	@$(NORMAL_INSTALL)
	test -z "$(bindir)" || $(MKDIR_P) "$(DESTDIR)$(bindir)"
	@list='$(bin_PROGRAMS)'; test -n "$(bindir)" || list=; \
	for p in $$list; do echo "$$p $$p"; done | \
	sed 's/$(EXEEXT)$$//' | \
	while read p p1; do if test -f $$p; \
	  then echo "$$p"; echo "$$p"; else :; fi; \
	done | \
	sed -e 'p;s,.*/,,;n;h' -e 's|.*|.|' \
	    -e 'p;x;s,.*/,,;s/$(EXEEXT)$$//;$(transform);s/$$/$(EXEEXT)/' | \
	sed 'N;N;N;s,\n, ,g' | \
	$(AWK) 'BEGIN { files["."] = ""; dirs["."] = 1 } \
	  { d=$$3; if (dirs[d] != 1) { print "d", d; dirs[d] = 1 } \
	    if ($$2 == $$4) files[d] = files[d] " " $$1; \
	    else { print "f", $$3 "/" $$4, $$1; } } \
	  END { for (d in files) print "f", d, files[d] }' | \
	while read type dir files; do \
	    if test "$$dir" = .; then dir=; else dir=/$$dir; fi; \
	    test -z "$$files" || { \
	      echo " $(INSTALL_PROGRAM_ENV) $(INSTALL_PROGRAM) $$files '$(DESTDIR)$(bindir)$$dir'"; \
	      $(INSTALL_PROGRAM_ENV) $(INSTALL_PROGRAM) $$files "$(DESTDIR)$(bindir)$$dir" || exit $$?; \
	    } \
	; done

uninstall-binPROGRAMS:
	@$(NORMAL_UNINSTALL)
	@list='$(bin_PROGRAMS)'; test -n "$(bindir)" || list=; \
	files=`for p in $$list; do echo "$$p"; done | \
	  sed -e 'h;s,^.*/,,;s/$(EXEEXT)$$//;$(transform)' \
	      -e 's/$$/$(EXEEXT)/' `; \
	test -n "$$list" || exit 0; \
	echo " ( cd '$(DESTDIR)$(bindir)' && rm -f" $$files ")"; \
	cd "$(DESTDIR)$(bindir)" && rm -f $$files

clean-binPROGRAMS:
	-test -z "$(bin_PROGRAMS)" || rm -f $(bin_PROGRAMS)

clean-noinstPROGRAMS:
	-test -z "$(noinst_PROGRAMS)" || rm -f $(noinst_PROGRAMS)
clm$(EXEEXT): $(clm_OBJECTS) $(clm_DEPENDENCIES) 
	@rm -f clm$(EXEEXT)
	$(LINK) $(clm_OBJECTS) $(clm_LDADD) $(LIBS)
clmdag$(EXEEXT): $(clmdag_OBJECTS) $(clmdag_DEPENDENCIES) 
	@rm -f clmdag$(EXEEXT)
	$(LINK) $(clmdag_OBJECTS) $(clmdag_LDADD) $(LIBS)
clmformat$(EXEEXT): $(clmformat_OBJECTS) $(clmformat_DEPENDENCIES) 
	@rm -f clmformat$(EXEEXT)
	$(LINK) $(clmformat_OBJECTS) $(clmformat_LDADD) $(LIBS)
clmframe$(EXEEXT): $(clmframe_OBJECTS) $(clmframe_DEPENDENCIES) 
	@rm -f clmframe$(EXEEXT)
	$(LINK) $(clmframe_OBJECTS) $(clmframe_LDADD) $(LIBS)
install-binSCRIPTS: $(bin_SCRIPTS)
	@$(NORMAL_INSTALL)
	test -z "$(bindir)" || $(MKDIR_P) "$(DESTDIR)$(bindir)"
	@list='$(bin_SCRIPTS)'; test -n "$(bindir)" || list=; \
	for p in $$list; do \
	  if test -f "$$p"; then d=; else d="$(srcdir)/"; fi; \
	  if test -f "$$d$$p"; then echo "$$d$$p"; echo "$$p"; else :; fi; \
	done | \
	sed -e 'p;s,.*/,,;n' \
	    -e 'h;s|.*|.|' \
	    -e 'p;x;s,.*/,,;$(transform)' | sed 'N;N;N;s,\n, ,g' | \
	$(AWK) 'BEGIN { files["."] = ""; dirs["."] = 1; } \
	  { d=$$3; if (dirs[d] != 1) { print "d", d; dirs[d] = 1 } \
	    if ($$2 == $$4) { files[d] = files[d] " " $$1; \
	      if (++n[d] == $(am__install_max)) { \
		print "f", d, files[d]; n[d] = 0; files[d] = "" } } \
	    else { print "f", d "/" $$4, $$1 } } \
	  END { for (d in files) print "f", d, files[d] }' | \
	while read type dir files; do \
	     if test "$$dir" = .; then dir=; else dir=/$$dir; fi; \
	     test -z "$$files" || { \
	       echo " $(INSTALL_SCRIPT) $$files '$(DESTDIR)$(bindir)$$dir'"; \
	       $(INSTALL_SCRIPT) $$files "$(DESTDIR)$(bindir)$$dir" || exit $$?; \
	     } \
	; done

uninstall-binSCRIPTS:
	@$(NORMAL_UNINSTALL)
	@list='$(bin_SCRIPTS)'; test -n "$(bindir)" || exit 0; \
	files=`for p in $$list; do echo "$$p"; done | \
	       sed -e 's,.*/,,;$(transform)'`; \
	test -n "$$list" || exit 0; \
	echo " ( cd '$(DESTDIR)$(bindir)' && rm -f" $$files ")"; \
	cd "$(DESTDIR)$(bindir)" && rm -f $$files

mostlyclean-compile:
	-rm -f *.$(OBJEXT)

distclean-compile:
	-rm -f *.tab.c

include ./$(DEPDIR)/clm.Po
include ./$(DEPDIR)/clmadjust.Po
include ./$(DEPDIR)/clmclose.Po
include ./$(DEPDIR)/clmdag.Po
include ./$(DEPDIR)/clmdist.Po
include ./$(DEPDIR)/clmformat.Po
include ./$(DEPDIR)/clmframe.Po
include ./$(DEPDIR)/clmimac.Po
include ./$(DEPDIR)/clminfo.Po
include ./$(DEPDIR)/clmmate.Po
include ./$(DEPDIR)/clmmeet.Po
include ./$(DEPDIR)/clmorder.Po
include ./$(DEPDIR)/clmps.Po
include ./$(DEPDIR)/clmresidue.Po
include ./$(DEPDIR)/report.Po

.c.o:
	$(COMPILE) -MT $@ -MD -MP -MF $(DEPDIR)/$*.Tpo -c -o $@ $<
	$(am__mv) $(DEPDIR)/$*.Tpo $(DEPDIR)/$*.Po
#	source='$<' object='$@' libtool=no \
#	DEPDIR=$(DEPDIR) $(CCDEPMODE) $(depcomp) \
#	$(COMPILE) -c $<

.c.obj:
	$(COMPILE) -MT $@ -MD -MP -MF $(DEPDIR)/$*.Tpo -c -o $@ `$(CYGPATH_W) '$<'`
	$(am__mv) $(DEPDIR)/$*.Tpo $(DEPDIR)/$*.Po
#	source='$<' object='$@' libtool=no \
#	DEPDIR=$(DEPDIR) $(CCDEPMODE) $(depcomp) \
#	$(COMPILE) -c `$(CYGPATH_W) '$<'`

ID: $(HEADERS) $(SOURCES) $(LISP) $(TAGS_FILES)
	list='$(SOURCES) $(HEADERS) $(LISP) $(TAGS_FILES)'; \
	unique=`for i in $$list; do \
	    if test -f "$$i"; then echo $$i; else echo $(srcdir)/$$i; fi; \
	  done | \
	  $(AWK) '{ files[$$0] = 1; nonempty = 1; } \
	      END { if (nonempty) { for (i in files) print i; }; }'`; \
	mkid -fID $$unique
tags: TAGS

TAGS:  $(HEADERS) $(SOURCES)  $(TAGS_DEPENDENCIES) \
		$(TAGS_FILES) $(LISP)
	set x; \
	here=`pwd`; \
	list='$(SOURCES) $(HEADERS)  $(LISP) $(TAGS_FILES)'; \
	unique=`for i in $$list; do \
	    if test -f "$$i"; then echo $$i; else echo $(srcdir)/$$i; fi; \
	  done | \
	  $(AWK) '{ files[$$0] = 1; nonempty = 1; } \
	      END { if (nonempty) { for (i in files) print i; }; }'`; \
	shift; \
	if test -z "$(ETAGS_ARGS)$$*$$unique"; then :; else \
	  test -n "$$unique" || unique=$$empty_fix; \
	  if test $$# -gt 0; then \
	    $(ETAGS) $(ETAGSFLAGS) $(AM_ETAGSFLAGS) $(ETAGS_ARGS) \
	      "$$@" $$unique; \
	  else \
	    $(ETAGS) $(ETAGSFLAGS) $(AM_ETAGSFLAGS) $(ETAGS_ARGS) \
	      $$unique; \
	  fi; \
	fi
ctags: CTAGS
CTAGS:  $(HEADERS) $(SOURCES)  $(TAGS_DEPENDENCIES) \
		$(TAGS_FILES) $(LISP)
	list='$(SOURCES) $(HEADERS)  $(LISP) $(TAGS_FILES)'; \
	unique=`for i in $$list; do \
	    if test -f "$$i"; then echo $$i; else echo $(srcdir)/$$i; fi; \
	  done | \
	  $(AWK) '{ files[$$0] = 1; nonempty = 1; } \
	      END { if (nonempty) { for (i in files) print i; }; }'`; \
	test -z "$(CTAGS_ARGS)$$unique" \
	  || $(CTAGS) $(CTAGSFLAGS) $(AM_CTAGSFLAGS) $(CTAGS_ARGS) \
	     $$unique

GTAGS:
	here=`$(am__cd) $(top_builddir) && pwd` \
	  && $(am__cd) $(top_srcdir) \
	  && gtags -i $(GTAGS_ARGS) "$$here"

distclean-tags:
	-rm -f TAGS ID GTAGS GRTAGS GSYMS GPATH tags

distdir: $(DISTFILES)
	@srcdirstrip=`echo "$(srcdir)" | sed 's/[].[^$$\\*]/\\\\&/g'`; \
	topsrcdirstrip=`echo "$(top_srcdir)" | sed 's/[].[^$$\\*]/\\\\&/g'`; \
	list='$(DISTFILES)'; \
	  dist_files=`for file in $$list; do echo $$file; done | \
	  sed -e "s|^$$srcdirstrip/||;t" \
	      -e "s|^$$topsrcdirstrip/|$(top_builddir)/|;t"`; \
	case $$dist_files in \
	  */*) $(MKDIR_P) `echo "$$dist_files" | \
			   sed '/\//!d;s|^|$(distdir)/|;s,/[^/]*$$,,' | \
			   sort -u` ;; \
	esac; \
	for file in $$dist_files; do \
	  if test -f $$file || test -d $$file; then d=.; else d=$(srcdir); fi; \
	  if test -d $$d/$$file; then \
	    dir=`echo "/$$file" | sed -e 's,/[^/]*$$,,'`; \
	    if test -d "$(distdir)/$$file"; then \
	      find "$(distdir)/$$file" -type d ! -perm -700 -exec chmod u+rwx {} \;; \
	    fi; \
	    if test -d $(srcdir)/$$file && test $$d != $(srcdir); then \
	      cp -fpR $(srcdir)/$$file "$(distdir)$$dir" || exit 1; \
	      find "$(distdir)/$$file" -type d ! -perm -700 -exec chmod u+rwx {} \;; \
	    fi; \
	    cp -fpR $$d/$$file "$(distdir)$$dir" || exit 1; \
	  else \
	    test -f "$(distdir)/$$file" \
	    || cp -p $$d/$$file "$(distdir)/$$file" \
	    || exit 1; \
	  fi; \
	done
check-am: all-am
check: check-am
all-am: Makefile $(PROGRAMS) $(SCRIPTS)
installdirs:
	for dir in "$(DESTDIR)$(bindir)" "$(DESTDIR)$(bindir)"; do \
	  test -z "$$dir" || $(MKDIR_P) "$$dir"; \
	done
install: install-am
install-exec: install-exec-am
install-data: install-data-am
uninstall: uninstall-am

install-am: all-am
	@$(MAKE) $(AM_MAKEFLAGS) install-exec-am install-data-am

installcheck: installcheck-am
install-strip:
	$(MAKE) $(AM_MAKEFLAGS) INSTALL_PROGRAM="$(INSTALL_STRIP_PROGRAM)" \
	  install_sh_PROGRAM="$(INSTALL_STRIP_PROGRAM)" INSTALL_STRIP_FLAG=-s \
	  `test -z '$(STRIP)' || \
	    echo "INSTALL_PROGRAM_ENV=STRIPPROG='$(STRIP)'"` install
mostlyclean-generic:

clean-generic:

distclean-generic:
	-test -z "$(CONFIG_CLEAN_FILES)" || rm -f $(CONFIG_CLEAN_FILES)
	-test . = "$(srcdir)" || test -z "$(CONFIG_CLEAN_VPATH_FILES)" || rm -f $(CONFIG_CLEAN_VPATH_FILES)

maintainer-clean-generic:
	@echo "This command is intended for maintainers to use"
	@echo "it deletes files that may require special tools to rebuild."
clean: clean-am

clean-am: clean-binPROGRAMS clean-generic clean-noinstPROGRAMS \
	mostlyclean-am

distclean: distclean-am
	-rm -rf ./$(DEPDIR)
	-rm -f Makefile
distclean-am: clean-am distclean-compile distclean-generic \
	distclean-tags

dvi: dvi-am

dvi-am:

html: html-am

html-am:

info: info-am

info-am:

install-data-am:

install-dvi: install-dvi-am

install-dvi-am:

install-exec-am: install-binPROGRAMS install-binSCRIPTS

install-html: install-html-am

install-html-am:

install-info: install-info-am

install-info-am:

install-man:

install-pdf: install-pdf-am

install-pdf-am:

install-ps: install-ps-am

install-ps-am:

installcheck-am:

maintainer-clean: maintainer-clean-am
	-rm -rf ./$(DEPDIR)
	-rm -f Makefile
maintainer-clean-am: distclean-am maintainer-clean-generic

mostlyclean: mostlyclean-am

mostlyclean-am: mostlyclean-compile mostlyclean-generic

pdf: pdf-am

pdf-am:

ps: ps-am

ps-am:

uninstall-am: uninstall-binPROGRAMS uninstall-binSCRIPTS

.MAKE: install-am install-strip

.PHONY: CTAGS GTAGS all all-am check check-am clean clean-binPROGRAMS \
	clean-generic clean-noinstPROGRAMS ctags distclean \
	distclean-compile distclean-generic distclean-tags distdir dvi \
	dvi-am html html-am info info-am install install-am \
	install-binPROGRAMS install-binSCRIPTS install-data \
	install-data-am install-dvi install-dvi-am install-exec \
	install-exec-am install-html install-html-am install-info \
	install-info-am install-man install-pdf install-pdf-am \
	install-ps install-ps-am install-strip installcheck \
	installcheck-am installdirs maintainer-clean \
	maintainer-clean-generic mostlyclean mostlyclean-compile \
	mostlyclean-generic pdf pdf-am ps ps-am tags uninstall \
	uninstall-am uninstall-binPROGRAMS uninstall-binSCRIPTS


clmformat.zmm.h: clmformat.zmm
	../../scripts/cfy clmformat.zmm clmformat.zmm.h defs_zmm

clmps.defs.h: clmps.defs
	../../scripts/cfy clmps.defs clmps.defs.h defs_ps

clmformat.c: clmformat.zmm.h

# Tell versions [3.59,3.63) of GNU make to not export all variables.
# Otherwise a system limit (for SysV at least) may be exceeded.
.NOEXPORT:
