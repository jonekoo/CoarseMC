dist_doc_DATA = README LICENSE DoxygenLayout.xml

include $(top_srcdir)/aminclude.am

SUBDIRS = claf90
SUBDIRS += cephes Examples
SUBDIRS += mt_stream_f90-1.11 src tests

if DX_COND_html

nodist_doc_DATA = @DX_DOCDIR@/html/index.html

$(nodist_doc_DATA): doxygen-doc

DX_CLEANFILES=@DX_DOCDIR@/html/* @DX_DOCDIR@/ptgbcyl.tag

endif

ACLOCAL_AMFLAGS = -I m4

EXTRA_DIST = $(DX_CONFIG)

MOSTLYCLEANFILES = $(DX_CLEANFILES)


