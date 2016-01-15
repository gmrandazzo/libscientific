Name:               libscientific
Version:            0.7.4
Release:            0%{?dist}
Summary:            Scientific Library for Dummies.
License:            GPLv3
Group:              System Environment/Libraries
Source:             %{name}-%{version}.tar.bz2
URL:                http://github.com/zeld/libscientific
BuildRequires:      cmake, gfortran

%description
 %{name} is a library born for computing statistics and numerical calculation.
 Libscientific provides tools for operation through matrix, vector and arrays,  make models (pca, pls, upca, upls),  polynomial interpolation
 and much more.

%package devel
Summary:        Development files for %{name}
Group:          Development/Libraries
Requires:       %{name}%{?_isa} = %{version}-%{release}

%description devel
The %{name}-devel package contains libraries, header files and documentation
for developing applications that use %{name}.

%prep
%setup -q -n %{name}-%{version}

%build
mkdir build
pushd build
%cmake -DCMAKE_BUILD_TYPE=Release \
       -DLIB_INSTALL_DIR=/usr/lib64/ \
       -DINCLUDE_INSTALL_DIR=/usr/include ..
make %{?_smp_mflags}

%install
rm -rf %{buildroot}
mkdir -p %{buildroot}
pushd build
make install DESTDIR=%{buildroot}

%post -p /sbin/ldconfig

%postun -p /sbin/ldconfig

%clean
rm -rf %{buildroot}

%files
%defattr(-,root,root,-)
#%doc COPYING README FORMAT
%{_libdir}/%{name}.so

%files devel
%defattr(-,root,root,-)
# %doc doc/html
%{_includedir}/scientific.h
%{_includedir}/scientific

%changelog
* Fri Sep  7 2012 Giuseppe Marco Randazzo <gmrandazzo@gmail.com> - 0.1.0
- New upstream release

