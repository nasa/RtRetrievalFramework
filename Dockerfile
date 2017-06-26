FROM centos:7

# Number of simultaneous build jobs for make
ENV BUILD_JOBS 4
ENV THIRDPARTY_SRC_DIR /src/third_party
ENV THIRDPARTY_INSTALL_DIR /rtr_framework/third_party
ENV ABSCO_DIR /absco
ENV MERRA_DIR /merra
ENV L2_SRC_DIR /src/l2_fp
ENV L2_BUILD_DIR /src/build
ENV L2_INSTALL_DIR /rtr_framework/install

VOLUME ["$ABSCO_DIR", "$MERRA_DIR", "$L2_SRC_DIR"]

# Get latest packages from Fedora Extra Packages for Enterprise Linux (EPEL) project
# https://fedoraproject.org/wiki/EPEL
RUN rpm -Uvh https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm

RUN yum install -y gcc gcc-gfortran gcc-c++ \
    make automake patch zlib-devel bzip2-devel \
    file which \
    ruby doxygen python2-future python-sphinx \
    python34 python34-devel python34-pip python34-numpy python34-nose

# For the initial build of the third party stuff copy over the source
ADD . $THIRDPARTY_SRC_DIR

# Build third party libraries to be included in image
# Don't use -j as an option to make as it screws up the build
RUN cd $THIRDPARTY_SRC_DIR && \
    ./configure THIRDPARTY=build PYTHON=python3 NOSETESTS=nosetests-3.4 PYTHON_VERSION=3.4 --prefix=$THIRDPARTY_INSTALL_DIR && \
    make thirdparty && \
    rm -rf $THIRDPARTY_SRC_DIR

# Compile software as default when running image
# Requires mounting of $L2_SRC_DIR
CMD mkdir -p $L2_BUILD_DIR && \
    cd $L2_BUILD_DIR && \
    $L2_SRC_DIR/configure THIRDPARTY=$THIRDPARTY_INSTALL_DIR PYTHON=python3 NOSETESTS=nosetests-3.4 PYTHON_VERSION=3.4 \
        --prefix=$L2_INSTALL_DIR --with-python-swig \
	--with-absco=$ABSCO_DIR --with-merra=$MERRA_DIR && \
    make all -j $BUILD_JOBS && make install
