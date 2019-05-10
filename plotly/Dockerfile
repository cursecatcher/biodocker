FROM r-base:3.5.3

RUN apt-get update
RUN apt-get install -y curl libssl-dev libgit2-dev
RUN apt-get install -y libcurl4-openssl-dev

RUN R -e "install.packages('curl')"
RUN R -e "install.packages('httr')"
RUN R -e "install.packages('gh')"
RUN R -e "install.packages('usethis')"

RUN R -e "install.packages('devtools')"
RUN R -e "install.packages(c('Rcpp','ggplot2','munsell','htmltools','DBI','assertthat','gridExtra','digest','fpc','TSP','registry','gclus','gplots','RColorBrewer','stringr','labeling','yaml'))"
RUN R -e "devtools::install_github('ropensci/plotly', ref='v4.8.0')"
RUN R -e "install.packages('heatmaply')"

COPY files/pandoc/* /usr/local/bin/

# Install our public gpg key to trusted store

# COPY files/conda.gpg /tmp/conda.gpg
# RUN install -o root -g root -m 644 /tmp/conda.gpg /etc/apt/trusted.gpg.d/
#
# # Add our debian repo
# RUN echo "deb [arch=amd64] https://repo.anaconda.com/pkgs/misc/debrepo/conda stable main" > /etc/apt/sources.list.d/conda.list
# RUN apt-get update && apt-get install -y conda
# RUN /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh"
# RUN /opt/conda/bin/conda install -c plotly plotly-orca psutil

RUN apt-get install -y python3

COPY files/heatmap.sh /bin/
COPY files/preprocessing.py /4seq/
COPY files/heatmap.R /4seq/


#COPY files/orca-1.2.1-x86_64.AppImage /tmp/
#
# RUN apt-get install -y libgtk2.0-0 libgconf-2-4
# RUN apt-get install -y fuse


# RUN ln -s /tmp/orca-1.2.1-x86_64.AppImage /bin/orca
