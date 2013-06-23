# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#

install-data-local: docbuild

all: docbuild
	
docbuild: common/builder.py build_options.py
	export SAGE_DOC_SRC="$$PWD"; \
	export SAGE_DOC="@datadir@/sage/doc"; \
	PYTHONPATH=common:$$PYTHONPATH $(PYTHON) $< --no-pdf-links all html

docenvconfdir=$(sysconfdir)/env.d
docenvconf_DATA = doc.sh

.PHONY: docbuild
