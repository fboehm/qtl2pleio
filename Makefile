WEB_PATH:= "~/Documents/website-mediumish/docs/static/software/qtl2pleio/."


all: build_site copy_site

.PHONY: build_site

build_site:
	Rscript -e 'pkgdown::build_site()'

copy_site: build_site
	cp -R docs/* $(WEB_PATH)
