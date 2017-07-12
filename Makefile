install:
	Rscript -e "rt::rmake()"

build:
	Rscript -e "rt::rbuild()"

check:
	Rscript -e "rt::rcheck()"

docs:
	Rscript -e "rt::rdoc()"

pkgdown:
	Rscript -e "rt::rpkgdown()"

winbuild:
	Rscript -e "rt::rwinbuild(devel=TRUE)"

all: install
