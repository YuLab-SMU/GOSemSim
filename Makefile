PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)

all: rd readme check clean

alldocs: rd readme site

rd:
	Rscript -e 'roxygen2::roxygenise(".")'

readme:
	Rscript -e 'rmarkdown::render("README.Rmd")'

build:
	cd ..;\
	R CMD build $(PKGSRC)

install:
	cd ..;\
	R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

check: build
	cd ..;\
	Rscript -e 'rcmdcheck::rcmdcheck("$(PKGNAME)_$(PKGVERS).tar.gz")'

check2: build
	cd ..;\
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz;\

bioccheck:
	cd ..;\
	Rscript -e 'BiocCheck::BiocCheck("$(PKGNAME)_$(PKGVERS).tar.gz")'

clean:
	cd ..;\
	$(RM) -r $(PKGNAME).Rcheck/

site:
	cd site_src;\
	ln -s ../../software/themes themes;\
	Rscript -e 'blogdown::build_site()';\
	rm themes;\
	cd ..

preview:
	cd site_src;\
	ln -s ../../software/themes themes;\
	Rscript -e 'blogdown::serve_site()';\
	rm themes;\
	cd ..


gitmaintain:
	git gc --auto;\
	git prune -v;\
	git fsck --full


update:
	git fetch --all;\
	git checkout master;\
	git merge upstream/master

release:
	git checkout RELEASE_3_6;\
	git fetch --all


push:
	git push upstream master;\
	git checkout github;\
	git merge -m 'merge from bioc repo' upstream/master;\
	git push -f origin HEAD:master;\
	git checkout master

