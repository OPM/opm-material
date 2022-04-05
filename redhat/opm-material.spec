#
# spec file for package opm-material
#

%define tag final
%define rtype release
%define toolset devtoolset-9
%define build_openmpi 1
%define build_openmpi3 1
%define build_mpich 1

Name:           opm-material
Version:        2018.10
Release:        0
Summary:        Open Porous Media - thermodynamic framework library
License:        GPL-3.0
Group:          Development/Libraries/C and C++
Url:            http://www.opm-project.org/
Source0:        https://github.com/OPM/%{name}/archive/release/%{version}/%{tag}.tar.gz#/%{name}-%{version}.tar.gz
BuildRequires: blas-devel lapack-devel
BuildRequires: git suitesparse-devel doxygen bc openblas-devel
BuildRequires: tinyxml-devel
BuildRequires: cmake3
BuildRequires: %{toolset}-toolchain
BuildRequires: boost-devel tbb-devel python3-devel
BuildRequires: dune-common-devel
BuildRequires: dune-istl-devel
BuildRequires: opm-common-devel

%if %{build_openmpi}
BuildRequires: openmpi-devel
BuildRequires: dune-common-openmpi-devel
BuildRequires: dune-istl-openmpi-devel
BuildRequires: opm-common-openmpi-devel
%endif

%if %{build_openmpi3}
BuildRequires: openmpi3-devel
BuildRequires: dune-common-openmpi3-devel
BuildRequires: dune-istl-openmpi3-devel
BuildRequires: opm-common-openmpi3-devel
%endif

%if %{build_mpich}
BuildRequires: mpich-devel
BuildRequires: dune-common-mpich-devel
BuildRequires: dune-istl-mpich-devel
BuildRequires: opm-common-mpich-devel
%endif

BuildRoot:      %{_tmppath}/%{name}-%{version}-build

%description
The Open Porous Media (OPM) initiative provides a set of open-source tools centered around the simulation of flow and transport of fluids in porous media. The goal of the initiative is to establish a sustainable environment for the development of an efficient and well-maintained software suite.

%package devel
Summary:        Development and header files for opm-material
Group:          Development/Libraries/C and C++

%description devel
This package contains the development and header files for opm-material

%package doc
Summary:        Documentation files for opm-material
Group:          Documentation
BuildArch:	noarch

%description doc
This package contains the documentation files for opm-material

%if %{build_openmpi}
%package openmpi-devel
Summary:        Development and header files for opm-material with openmpi
Group:          Development/Libraries/C and C++

%description openmpi-devel
This package contains the development and header files for opm-material with
openMPI.
%endif

%if %{build_openmpi3}
%package openmpi3-devel
Summary:        Development and header files for opm-material with openmpi3
Group:          Development/Libraries/C and C++

%description openmpi3-devel
This package contains the development and header files for opm-material with
openMPI3.
%endif

%if %{build_mpich}
%package mpich-devel
Summary:        Development and header files for opm-material with mpich
Group:          Development/Libraries/C and C++

%description mpich-devel
This package contains the development and header files for opm-material with
mpich.
%endif

%global debug_package %{nil}

%prep
%setup -q -n %{name}-%{rtype}-%{version}-%{tag}

# consider using -DUSE_VERSIONED_DIR=ON if backporting
%build
mkdir serial
pushd serial
scl enable %{toolset} 'cmake3 -DENABLE_MPI=0 -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=RelWithDebInfo -DSTRIP_DEBUGGING_SYMBOLS=ON -DCMAKE_INSTALL_PREFIX=%{_prefix} -DCMAKE_INSTALL_DOCDIR=share/doc/%{name}-%{version} -DUSE_RUNPATH=OFF -DWITH_NATIVE=OFF  ..'
scl enable %{toolset} 'make %{?_smp_mflags}'
scl enable %{toolset} 'make test'
popd

%if %{build_openmpi}
mkdir openmpi
pushd openmpi
module load mpi/openmpi-x86_64
scl enable %{toolset} 'cmake3 -DUSE_MPI=1 -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=RelWithDebInfo -DSTRIP_DEBUGGING_SYMBOLS=ON -DCMAKE_INSTALL_PREFIX=%{_prefix}/lib64/openmpi -DCMAKE_INSTALL_LIBDIR=lib -DUSE_RUNPATH=OFF -DWITH_NATIVE=OFF  -DCMAKE_INSTALL_INCLUDE_DIR=%{_prefix}/include/openmpi-x86_64 ..'
scl enable %{toolset} 'make %{?_smp_mflags}'
scl enable %{toolset} 'make test'
module unload mpi/openmpi-x86_64
popd
%endif

%if %{build_openmpi3}
mkdir openmpi3
pushd openmpi3
module load mpi/openmpi3-x86_64
scl enable %{toolset} 'cmake3 -DUSE_MPI=1 -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=RelWithDebInfo -DSTRIP_DEBUGGING_SYMBOLS=ON -DCMAKE_INSTALL_PREFIX=%{_prefix}/lib64/openmpi3 -DCMAKE_INSTALL_LIBDIR=lib -DUSE_RUNPATH=OFF -DWITH_NATIVE=OFF  -DCMAKE_INSTALL_INCLUDE_DIR=%{_prefix}/include/openmpi3-x86_64 ..'
scl enable %{toolset} 'make %{?_smp_mflags}'
scl enable %{toolset} 'make test'
module unload mpi/openmpi3-x86_64
popd
%endif

%if %{build_mpich}
mkdir mpich
pushd mpich
module load mpi/mpich-x86_64
scl enable %{toolset} 'cmake3 -DUSE_MPI=1 -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=RelWithDebInfo -DSTRIP_DEBUGGING_SYMBOLS=ON -DCMAKE_INSTALL_PREFIX=%{_prefix}/lib64/mpich -DCMAKE_INSTALL_LIBDIR=lib -DUSE_RUNPATH=OFF -DWITH_NATIVE=OFF  -DCMAKE_INSTALL_INCLUDE_DIR=%{_prefix}/include/mpich-x86_64 ..'
scl enable %{toolset} 'make %{?_smp_mflags}'
scl enable %{toolset} 'make test'
%endif

%install
scl enable %{toolset} 'make install DESTDIR=${RPM_BUILD_ROOT} -C serial'
scl enable %{toolset} 'make install-html DESTDIR=${RPM_BUILD_ROOT} -C serial'

%if %{build_openmpi}
scl enable %{toolset} 'make install DESTDIR=${RPM_BUILD_ROOT} -C openmpi'
mv ${RPM_BUILD_ROOT}/%{_libdir}/openmpi/include/* ${RPM_BUILD_ROOT}/usr/include/openmpi-x86_64/
%endif

%if %{build_openmpi3}
scl enable %{toolset} 'make install DESTDIR=${RPM_BUILD_ROOT} -C openmpi3'
mv ${RPM_BUILD_ROOT}/%{_libdir}/openmpi3/include/* ${RPM_BUILD_ROOT}/usr/include/openmpi3-x86_64/
%endif

%if %{build_mpich}
scl enable %{toolset} 'make install DESTDIR=${RPM_BUILD_ROOT} -C mpich'
mv ${RPM_BUILD_ROOT}/%{_libdir}/mpich/include/* ${RPM_BUILD_ROOT}/usr/include/mpich-x86_64/
%endif

%clean
rm -rf %{buildroot}

%files doc
%{_docdir}/*

%files devel
%defattr(-,root,root,-)
/usr/lib/dunecontrol/*
/usr/lib/pkgconfig/*
%{_includedir}/*
%{_datadir}/cmake/*
%{_datadir}/opm/cmake/Modules/*

%if %{build_openmpi}
%files openmpi-devel
%defattr(-,root,root,-)
%{_libdir}/openmpi/lib/dunecontrol/*
%{_libdir}/openmpi/lib/pkgconfig/*
%{_includedir}/openmpi-x86_64/*
%{_libdir}/openmpi/share/cmake/*
%{_libdir}/openmpi/share/opm/cmake/Modules/*
%endif

%if %{build_openmpi3}
%files openmpi3-devel
%defattr(-,root,root,-)
%{_libdir}/openmpi3/lib/dunecontrol/*
%{_libdir}/openmpi3/lib/pkgconfig/*
%{_includedir}/openmpi3-x86_64/*
%{_libdir}/openmpi3/share/cmake/*
%{_libdir}/openmpi3/share/opm/cmake/Modules/*
%endif

%if %{build_mpich}
%files mpich-devel
%defattr(-,root,root,-)
%{_libdir}/mpich/lib/dunecontrol/*
%{_libdir}/mpich/lib/pkgconfig/*
%{_includedir}/mpich-x86_64/*
%{_libdir}/mpich/share/cmake/*
%{_libdir}/mpich/share/opm/cmake/Modules/*
%endif
