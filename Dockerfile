FROM rocker/r-ver:latest
RUN apt-get update && \
    apt-get -y upgrade

LABEL org.opencontainers.image.licenses="GPL-2.0-or-later"
LABEL org.opencontainers.image.source="https://github.com/Izar-de-villasante/SG0002"
LABEL org.opencontainers.image.vendor="IJC Bioinformatics Team"
LABEL org.opencontainers.image.authors="Izar de Villasante <idevillasante@carrerasresearch.org>"
LABEL org.opencontainers.image.description="Ready to use rstudio + quarto container to start your new projects. This image contains Rv.4.2 Python v.3.8+ rstudio v2.1.0.2 renv 0.17.3 shiny Bioconductor and quarto 1.3+ and the extensions shinylive and molstar."

ENV GCM_CREDENTIAL_STORE=gpg
ENV DEBIAN_FRONTEND noninteractive
ENV S6_VERSION=v2.1.0.2
ENV RSTUDIO_VERSION=latest
#2022.07.2+576
ENV DEFAULT_USER=rstudio
ENV PANDOC_VERSION=default
ENV QUARTO_VERSION=default

ENV RENV_VERSION 0.17.3
ENV RENV_PATHS_LIBRARY .cache/renv/


RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    git git-core \
    pass \
    make \
    libffi-dev \
    python3-dev python3-pip \
    build-essential \
    libmemcached-dev \
    curl \
    libxml2-dev \
    g++ \
    gcc \
    libxml2 \
    libjpeg-dev \
    libfreetype6-dev \
    webp \
    pcre++-dev \
    libpango1.0-dev \
    libxslt-dev \
    libcairo2-dev \
    libsqlite3-dev \
    libmariadbd-dev \
    libpq-dev \
    libev-dev \
    libssh2-1-dev \
    unixodbc-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    coinor-libcbc-dev coinor-libclp-dev libglpk-dev \
    libgtk2.0-dev libxt-dev xvfb xauth xfonts-base \
    ssh

#RUN "Y" | pip3 install shynilive --upgrade
RUN install2.r --error --skipinstalled --ncpus -1 \
    --repos https://ropensci.r-universe.dev --repos getOption \
    renv \
    BiocManager \
    devtools \
    data.table \
    shiny \
    Cairo \
    docopt \
    future \
    future.batchtools \
    && rm -rf /tmp/downloaded_packages \
    && strip /usr/local/lib/R/site-library/*/libs/*.so

RUN /rocker_scripts/install_jupyter.sh
RUN /rocker_scripts/install_rstudio.sh
RUN /rocker_scripts/install_pandoc.sh
RUN /rocker_scripts/install_quarto.sh 'prerelease'
RUN /rocker_scripts/install_python.sh
RUN /rocker_scripts/install_shiny_server.sh
RUN quarto install extension jmbuhr/quarto-molstar --no-prompt
RUN quarto install extension quarto-ext/shinylive --no-prompt

COPY renv.lock renv.lock
RUN mkdir -p renv
COPY .Rprofile .Rprofile
COPY renv/activate.R renv/activate.R
COPY renv/settings.json renv/settings.json

RUN R -e "renv::restore()"

EXPOSE 8787
RUN useradd -ms /bin/bash izar
RUN adduser izar sudo && echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

#RUN useradd --create-home --shell /bin/bash user1
#RUN echo 'user1:password' | chpasswd
#RUN apt -y install ssh
#RUN apt-get install ssh-client

EXPOSE 3838
EXPOSE 22
CMD ["/init"]
