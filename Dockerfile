FROM ubuntu:20.04

RUN apt-get update && apt-get install -y build-essential git zlib1g-dev wget \ 
    libbz2-dev liblzma-dev libcurl4-gnutls-dev autoconf libncurses5-dev \
    libncursesw5-dev libssl-dev libxml-xpath-perl libjson-perl bedtools \
    python3 python3-pip

RUN pip install pandas

RUN git clone -b 1.13 https://github.com/samtools/htslib.git && \
    cd htslib && git submodule update --init --recursive && make -j `nproc` && make install && \
    cd .. && rm -rf htslib

RUN git clone -b 1.13 https://github.com/samtools/bcftools.git && \
    cd bcftools && autoheader && autoconf && ./configure && make -j `nproc` && make install && \
    cd .. && rm -rf bcftools

RUN git clone -b 1.13 https://github.com/samtools/samtools.git && \
    cd samtools && autoheader && autoconf && ./configure && make -j `nproc` && make install && \
    cd .. && rm -rf samtools

RUN git clone https://github.com/gariem/SURVIVOR.git && cd SURVIVOR/Debug && make -j `nproc` && \
    ln /SURVIVOR/Debug/SURVIVOR /sbin/SURVIVOR

ENV LD_LIBRARY_PATH=/usr/local/lib

CMD ["bash"]