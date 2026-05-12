FROM python:3.12-slim

ARG USER="repro"

# System deps for R
RUN apt-get update && \
	apt-get install -y --no-install-recommends \
		latexmk \
		libcurl4-openssl-dev \
		libfontconfig1-dev \
		libfribidi-dev \
		libharfbuzz-dev \
		libssl-dev \
		libtiff5-dev \
		libuv1-dev \
		libxml2-dev \
		r-base \
		r-base-dev \
		sudo \
		texlive-base \
		texlive-latex-extra \
		texlive-luatex \
		texlive-pictures \
		texlive-publishers \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir ~/.R
RUN echo MAKEFLAGS = -j$(($(nproc)/4)) > ~/.R/Makevars

# Install R packages in system-wide context
RUN Rscript -e "install.packages(c('tidyverse', 'scales', 'patchwork', 'tikzDevice', 'ggnewscale'), Ncpus=parallel::detectCores(), repos='https://cloud.r-project.org', clean=TRUE)"


# Install python dependencies in system-wide context
WORKDIR /root
COPY ./reproduction/requirements.txt .
RUN pip install --no-cache-dir -r ./requirements.txt
RUN rm requirements.txt

# Create a local user
RUN useradd -m ${USER} && \
    echo "${USER} ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers
USER ${USER}
WORKDIR /home/${USER}

# Drop a bash instead of python
CMD [ "bash" ]
