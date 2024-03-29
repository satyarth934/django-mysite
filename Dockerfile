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
        vim \
        git \
        wget \
        unzip \
        python3-dev \
        mariadb-client \
        sudo \
        # libpq-dev \
        # build-essential \
        # cmake \
        # sqlite3 \
        # libsqlite3-dev \
        # libboost-dev \
        # libboost-all-dev \
        # libboost-system1.67-dev \
        # libboost-thread1.67-dev \
        # libboost-serialization1.67-dev \
        # libboost-python1.67-dev \
        # libboost-regex1.67-dev \
        # libboost-iostreams1.67-dev \
        # libxrender1 \
        # libxext6 \ 
        # libeigen3-dev \
        # openjdk-11-jdk \
        # openjdk-11-jre \
        # ncbi-blast+ \
	    # libigraph0v5 \
	    # libigraph0-dev \
        # zlib1g-dev \
        # libcairo2-dev \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

ADD ./requirements.txt /root/
RUN pip3 install -U pip setuptools
RUN pip3 install -r requirements.txt

# Install RetroTide API
COPY retrotide-master /root/retrotide
RUN cd /root/retrotide && \
    pip3 install . && \
    cd /root

# Clone Website code
RUN git clone -b spin-setup https://github.com/satyarth934/django-mysite.git
# RUN git clone -b db-integration https://github.com/satyarth934/django-mysite.git

# Adding MySQL config file to the container
# ADD ./mysql_config_spin.yaml /root/django-mysite/mysql_config.yaml

WORKDIR /root/django-mysite


# The '--insecure' flag is needed to render the 'static' and 'media' files if DEBUG=False.
# TODO: Find a better alternative to this for production website.
CMD python3 manage.py makemigrations; python3 manage.py migrate; python3 manage.py runserver 0.0.0.0:8000 --insecure
