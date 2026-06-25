SAGE_PATH=$(shell dirname $(shell readlink -f $(shell which sage)))
$(info Sage Path: $(SAGE_PATH))

SAGE=$(SAGE_PATH)/sage

ifeq (,$(wildcard $(SAGE_PATH)/local/bin/python3))
PYTHON_PATH=python3
PYTHON_SYSTEM=true
else
PYTHON_PATH=$(SAGE_PATH)/local/bin/python3
PYTHON_SYSTEM=false
endif

#PYTHON_VERSION=3.9
$(info Python path: $(PYTHON_PATH))
PYTHON_VERSION=$(shell $(PYTHON_PATH) -c 'import platform; major, minor, patch = platform.python_version_tuple(); print(major + "." + minor)')
$(info Python version: $(PYTHON_VERSION))

PYTHON_SITE_PACKAGES_INCLUDE=$(shell $(PYTHON_PATH) -c "import site; print('-I' + ' -I'.join(site.getsitepackages()))")
#PYTHON_SITE_PACKAGES_PATH=/usr/local/lib/python$(PYTHON_VERSION)/site-packages
$(info Python site packages includes: $(PYTHON_SITE_PACKAGES_INCLUDE))

PYTHON_INCLUDE=$(shell $(PYTHON_PATH) -c "from sysconfig import get_paths as gp; print(gp()['include'])")

#CYTHON=$(SAGE) -cython -I $(SAGE_PATH)/local/lib64/python$(PYTHON_VERSION)/site-packages
ifeq ($(PYTHON_SYSTEM), true)
ifeq (, $(shell which cython3))
CYTHON=cython $(PYTHON_SITE_PACKAGES_INCLUDE)
else
CYTHON=cython3 $(PYTHON_SITE_PACKAGES_INCLUDE)
endif
else
CYTHON=$(SAGE) -cython $(PYTHON_SITE_PACKAGES_INCLUDE)
endif
$(info Cython: $(CYTHON))

ifeq ($(PYTHON_SYSTEM), true)
SAGE_INCLUDES=$(shell $(PYTHON_PATH) -c "import site; print('-I' + ' -I'.join([path + '/sage/cpython' for path in site.getsitepackages()]))") \
			  $(shell $(PYTHON_PATH) -c "import site; print('-I' + ' -I'.join([path + '/sage/libs/ntl' for path in site.getsitepackages()]))")
else
PYTHON_INCLUDE=$(SAGE_PATH)/local/include/python$(PYTHON_VERSION)
SAGE_SRC=$(SAGE_PATH)/src
SAGE_INCLUDES=-I$(SAGE_SRC) -I$(SAGE_SRC)/sage/libs/ntl -I$(SAGE_SRC)/sage/cpython
endif
$(info Python include: $(PYTHON_INCLUDE))
$(info Sage include: $(SAGE_INCLUDES))

CC=gcc -shared -pthread -fPIC -fwrapv -g -O2 -Wall -fno-strict-aliasing \
	   -I$(PYTHON_INCLUDE) \
	   $(PYTHON_SITE_PACKAGES_INCLUDE) \
	   $(SAGE_INCLUDES)

all: \
LEE.so

LEE.so: LEE.c
	$(CC) -o LEE.so LEE.c

LEE.c: LEE.pyx
	$(CYTHON) -3 LEE.pyx

LEE.pyx: LEE.sage.py
	cp LEE.sage.py LEE.pyx

LEE.sage.py: LEE.sage
	$(SAGE) -preparse LEE.sage

clean:
	rm -f LEE.so LEE.c LEE.pyx
	rm -f *.sage.py
