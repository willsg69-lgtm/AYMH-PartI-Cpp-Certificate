# AYMH Part I v6 CAPD certificate driver
#
# Proof target:
#   make mp
#
# The proof run uses CAPD's multiprecision MPFR interval backend.  If
# mpcapd-config is not on PATH, point MPCAPD_CONFIG to the installed binary:
#
#   MPCAPD_CONFIG=/path/to/capd/bin/mpcapd-config make mp

MPCAPD_CONFIG ?= mpcapd-config

# Apple's default make sets CXX=c++ (clang).  If CAPD was configured with
# Homebrew GCC, linking user programs with clang can fail because of C++ ABI
# differences.  Prefer the newest Homebrew g++ when available, while still
# allowing `make CXX=/path/to/compiler` to override this.
ifeq ($(origin CXX), default)
BREW_GXX := $(shell if command -v brew >/dev/null 2>&1; then ls "$$(brew --prefix gcc)"/bin/g++-* 2>/dev/null | tail -1; fi)
CXX := $(if $(BREW_GXX),$(BREW_GXX),g++)
endif

CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -pedantic
# Homebrew installs GMP/MPFR outside GCC's default header and library search
# paths on macOS.  CAPD's development mpcapd-config may list -lgmp/-lmpfr
# without those paths, so add them when the packages are available.
BREW_MP_FLAGS := $(shell if command -v brew >/dev/null 2>&1; then \
	gmp="$$(brew --prefix gmp 2>/dev/null)"; \
	mpfr="$$(brew --prefix mpfr 2>/dev/null)"; \
	if [ -n "$$gmp" ] && [ -n "$$mpfr" ]; then \
		echo "-I$$gmp/include -I$$mpfr/include -L$$gmp/lib -L$$mpfr/lib"; \
	fi; \
fi)
CXXFLAGS += $(BREW_MP_FLAGS)
MPCAPD_FLAGS := $(shell $(MPCAPD_CONFIG) --cflags 2>/dev/null)
MPCAPD_LIBS  := $(shell $(MPCAPD_CONFIG) --libs 2>/dev/null)

TARGET_V6_MP := lppssi_cert_v6_mp
SRC_V6 := src/lppssi_cert_v6.cpp

all: mp
mp: $(TARGET_V6_MP)

$(TARGET_V6_MP): $(SRC_V6) Makefile
	@if [ -z "$(MPCAPD_FLAGS)" ] || [ -z "$(MPCAPD_LIBS)" ]; then \
		echo "mpcapd-config was not found or did not return flags."; \
		echo "Install CAPD with MPFR support, then run:"; \
		echo "  MPCAPD_CONFIG=/path/to/install/bin/mpcapd-config make mp"; \
		exit 1; \
	fi
	$(CXX) $(CXXFLAGS) -DLPPSSI_USE_MP $(MPCAPD_FLAGS) -Iinclude $< $(MPCAPD_LIBS) -o $@

print-cxx:
	@printf '%s\n' "$(CXX)"

clean:
	rm -f $(TARGET_V6_MP) capd_certificate_v6_220bit_local.txt

.PHONY: all mp clean print-cxx
