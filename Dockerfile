FROM debian:10

USER root

ENV LANG C.UTF-8
WORKDIR /root

RUN apt-get update \
    && apt-get install -y \
        python3-pip \
        python3-all \
        python3-all-dev \
        python3-tk \
        libpq-dev \
        vim \
        git \
        wget \
        unzip \
        build-essential \
        cmake \
        python3-dev \
        sqlite3 \
        libsqlite3-dev \
        libboost-dev \
        libboost-all-dev \
        libboost-system1.67-dev \
        libboost-thread1.67-dev \
        libboost-serialization1.67-dev \
        libboost-python1.67-dev \
        libboost-regex1.67-dev \
        libboost-iostreams1.67-dev \
        libxrender1 \
        libxext6 \ 
        libeigen3-dev \
        openjdk-11-jdk \
        openjdk-11-jre \
        sudo \
        ncbi-blast+ \
	    libigraph0v5 \
	    libigraph0-dev \
        zlib1g-dev \
        libcairo2-dev \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

ADD ./requirements.txt /root/
RUN pip3 install -U pip setuptools
RUN pip3 install -r requirements.txt
RUN pip3 install django

RUN git clone -b spin-setup https://github.com/satyarth934/django-mysite.git
WORKDIR /root/django-mysite

