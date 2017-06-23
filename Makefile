# A simple makefile for creating the oxyfuel (ASU optimization) distribution zip file
VERSION=$(shell git describe --tags --dirty)
PRODUCT=ASU Optimization (oxyfuel)
PROD_SNAME=Oxyfuel
LICENSE=LICENSE.md
PKG_DIR=CCSI_$(PROD_SNAME)_$(VERSION)
PACKAGE=$(PKG_DIR).zip

PAYLOAD=docs/*.pdf \
        ASUSpecFiles \
	DegeneracyExamples \
	InitFiles \
	ModelFiles \
	ProcsFiles \
	ResultInputs \
	ThrottleValveExampleSpecFiles \
	worker0 \
	*.gms \
	*.m \
	README.md \
	$(LICENSE)

# Get just the top part (not dirname) of each entry so cp -r does the right thing
PAYLOAD_TOPS=$(foreach v,$(PAYLOAD),$(shell echo $v | cut -d'/' -f1))
# And the payload with the PKG_DIR prepended
PKG_PAYLOAD=$(addprefix $(PKG_DIR)/, $(PAYLOAD))

# OS detection & changes
UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
  MD5BIN=md5sum
endif
ifeq ($(UNAME), Darwin)
  MD5BIN=md5
endif
ifeq ($(UNAME), FreeBSD)
  MD5BIN=md5
endif

.PHONY: all test clean

all: $(PACKAGE)

# Make zip file without extra file attribs (-X) so md5sum doesn't
# change if the payload hasn't
$(PACKAGE): $(PAYLOAD)
	@mkdir $(PKG_DIR)
	@cp -r $(PAYLOAD_TOPS) $(PKG_DIR)
	@zip -rqX $(PACKAGE) $(PKG_PAYLOAD)
	@$(MD5BIN) $(PACKAGE)
	@rm -rf $(PKG_DIR)

test:
	@$(MAKE) -sC test test

clean:
	@$(MAKE) -sC test clean
	@rm -rf $(PACKAGE) $(PKG_DIR) *.zip *~
