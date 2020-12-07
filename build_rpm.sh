#!/usr/bin/bash

# clean all .pyc files, if any
find . -name "*.pyc" -exec rm -f {} \;

# create ~/rpmbuild tree
rpmdev-setuptree

# copy source
rsync -av --exclude=".*" ./ ~/rpmbuild/SOURCES

# copy spec file
cp ./SMHI-gridpp-lib.spec ~/rpmbuild/SPECS

# build rpm
rpmbuild --target=x86_64 -bb ~/rpmbuild/SPECS/SMHI-gridpp-lib.spec
