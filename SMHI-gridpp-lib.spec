Name: SMHI-gridpp-lib
Version: 0.1
Release: 1
Summary: Driver library for SMHI's GridPP version
License: SMHI
Group: Applications/System
Requires: swig, boost, gsl, netcdf, armadillo
Requires: python3, python3-scipy, python3-numpy, python3-six

%define AppOwner        gridpp
%define AppGroup        gridppg
%define AppOwnerTst     gridpp.t
%define AppGroupTst     gridppgt
%define AppOwnerUtv     gridpp.u
%define AppGroupUtv     gridppgu
%define INSTALLDIR %{_libdir}/python3.6/site-packages/smhi_gridpp_lib
%define _unpackaged_files_terminate_build 0

%description

Gridpp a is post-processing tool for gridded weather forecasts. It consists of a library of
commonly-used methods and a command-line tool that applies these methods to forecast fields in NetCDF files.
Gridpp is written in C++ but offers python bindings to the functions in the library. The tool is used at
MET Norway to produce operational weather forecasts for Yr (https://www.yr.no).
Gridpp is currently under active development and the current version is a prototype for testing. Feedback is
welcome, either by using the issue tracker in Github, or by contacting Thomas Nipen (thomasn@met.no).

This is an SMHI's adaptation of MetNO's library

%install

mkdir -p $RPM_BUILD_ROOT%{INSTALLDIR}

cp ${RPM_SOURCE_DIR}/build/extras/SWIG/python/gridpp.py $RPM_BUILD_ROOT%{INSTALLDIR}
cp ${RPM_SOURCE_DIR}/extras/python/__init__.py $RPM_BUILD_ROOT%{INSTALLDIR}
cp ${RPM_SOURCE_DIR}/build/extras/SWIG/python/_gridpp.so $RPM_BUILD_ROOT%{INSTALLDIR}

%postun

%clean

rm -rf ${RPM_BUILD_ROOT}
rm -rf ${RPM_SOURCE_DIR}

%files
%{INSTALLDIR}/*.py
%{INSTALLDIR}/_gridpp.so
%{INSTALLDIR}/__pycache__/*.pyc

%defattr(755,root,root,755)
%attr(644,root,root) %{INSTALLDIR}/*.py

%changelog

* Thu Oct 15 2020 Aliaksandr Rahachou <aliaksandr.rahachou@hiq.se> - 0.1-1
- First variant, GridPP 0.4.2b