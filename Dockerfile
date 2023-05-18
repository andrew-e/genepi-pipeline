FROM python:3.11-slim

RUN apt update && \
    apt install -y --no-install-recommends build-essential software-properties-common dirmngr gnupg libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libfreetype6-dev libharfbuzz-dev libfribidi-dev libxml2-dev libcurl4-gnutls-dev libssl-dev libgmp-dev libnlopt-dev cmake libcairo2-dev libxt-dev && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys B8F25A8A73EACF41 && \
    echo "deb https://www.stats.bris.ac.uk/R/bin/linux/debian bullseye-cran40/" | tee -a /etc/apt/sources.list && \
    apt update && \
    apt install -y r-base=4.3.0-1~bullseyecran.0 r-base-dev=4.3.0-1~bullseyecran.0 && \ 
    rm -rf /var/lib/apt/lists/* #pin versions?

COPY . .
RUN Rscript docker/requirements.r && \
    pip install -r docker/requirements.txt

CMD ["/bin/bash"]
