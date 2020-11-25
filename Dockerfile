FROM ubuntu:18.04

# Install requirements
RUN apt-get update && apt-get install -y \
	python3 \
	python3-distutils \
	wget \
	cmake \
	tar


# Install python libraries
RUN wget "https://bootstrap.pypa.io/get-pip.py" && \
	python3 get-pip.py
COPY requirements.txt /tmp/
RUN pip3 install --requirement /tmp/requirements.txt

# Install C++ compiler
RUN apt-get install -y g++

# Install GROMACS
RUN apt-get install -y gromacs

# Copy modules and test files
RUN mkdir /home/jcalc/
COPY . /home/jcalc

# Turn scripts into executables
RUN chmod +x /home/jcalc/JCalc.py
WORKDIR /home/data/
ENTRYPOINT ["/home/jcalc/JCalc.py"]
CMD ["-h"]
