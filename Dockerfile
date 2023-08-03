FROM python:3.11-slim

RUN apt update && \
    apt install -y --no-install-recommends r-base=4.2.2.20221110-2 r-base-dev=4.2.2.20221110-2 \
    build-essential software-properties-common dirmngr gnupg pandoc \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libfreetype6-dev libharfbuzz-dev libfribidi-dev \ 
    libxml2-dev libcurl4-gnutls-dev libssl-dev libgmp-dev libnlopt-dev cmake libcairo2-dev libxt-dev \
    plink1.9 sqlite3 wget vim ripgrep && \
    mkdir -p /home/r_scripts && \
    rm -rf /var/lib/apt/lists/*

COPY docker/ docker
RUN Rscript docker/requirements.r
RUN pip install -r docker/requirements.txt

COPY r_scripts /home/r_scripts
#RUN wget http://csg.sph.umich.edu/abecasis/metal/download/Linux-metal.tar.gz
#RUN install bcftools too?

CMD ["/bin/bash"]
