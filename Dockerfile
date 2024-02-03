R/gwas_formatting.rFROM andrewrrelmore/genepi_pipeline_base:latest

COPY scripts /home/scripts
COPY R /home/R
WORKDIR "/home/scripts"

CMD ["/bin/bash"]
