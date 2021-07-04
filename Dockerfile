FROM ubuntu:18.04

# Install requirements
RUN apt-get update && apt-get install -y \
	python3.8 \
	python3.8-distutils \
	python3-pip \
	wget \
	cmake \
	tar

# Install C++ compiler
RUN apt-get install -y g++

# Install GROMACS
RUN apt-get install -y gromacs

# Install obabel
RUN wget https://github.com/openbabel/openbabel/archive/refs/tags/openbabel-3-1-1.tar.gz && \
	tar -xzvf openbabel-3-1-1.tar.gz && \
	cd openbabel-openbabel-3-1-1/ && \
	mkdir build && \
	cd build && \
	cmake -DRUN_SWIG=ON -DPYTHON_BINDINGS=ON .. && \
	make -j2 && \
	make install

# Install jcalc
RUN mkdir /home/jcalc/
COPY . /home/jcalc
RUN cd /home/jcalc && ./easy_install.py

# Turn scripts into executables
WORKDIR /home/data
ENTRYPOINT ["jcalc"]
CMD ["-h"]
