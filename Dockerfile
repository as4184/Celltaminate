FROM rocker/r-ver:4.3.3

LABEL org.opencontainers.image.title="Celltaminate"
LABEL org.opencontainers.image.description="Contamination-aware prioritization of microbial taxa from Kraken reports"

RUN apt-get update && apt-get install -y --no-install-recommends         python3         python3-pip         libcurl4-openssl-dev         libssl-dev         libxml2-dev         zlib1g-dev         && rm -rf /var/lib/apt/lists/*

RUN Rscript -e "install.packages(c('dplyr','ggplot2','ggrepel','fmsb','writexl','tibble','tidyr','shiny','shinyWidgets','DT','httr','jsonlite','zip','xml2'), repos='https://cloud.r-project.org')"

WORKDIR /opt/Celltaminate
COPY . /opt/Celltaminate
ENV PATH="/opt/Celltaminate/bin:${PATH}"

ENTRYPOINT ["celltaminate"]
CMD ["--help"]
