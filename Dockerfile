FROM python:3.11-slim

RUN apt update && \
    apt install -y --no-install-recommends build-essential software-properties-common dirmngr gnupg \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libfreetype6-dev libharfbuzz-dev libfribidi-dev \ 
    libxml2-dev libcurl4-gnutls-dev libssl-dev libgmp-dev libnlopt-dev cmake libcairo2-dev libxt-dev \
    plink1.9 sqlite3 && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys B8F25A8A73EACF41 && \
    echo "deb https://www.stats.bris.ac.uk/R/bin/linux/debian bullseye-cran40/" | tee -a /etc/apt/sources.list && \
    apt update && \
    apt install -y r-base=4.3.0-1~bullseyecran.0 r-base-dev=4.3.0-1~bullseyecran.0 && \ 
    rm -rf /var/lib/apt/lists/* #pin versions?

COPY docker/ docker
RUN Rscript docker/requirements.r && \
    pip install -r docker/requirements.txt

RUN mkdir -p /home/r_scripts
COPY r_scripts /home/r_scripts
#RUN wget http://csg.sph.umich.edu/abecasis/metal/download/Linux-metal.tar.gz
#RUN install bcftools too?

CMD ["/bin/bash"]
