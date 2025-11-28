# Use a multi-arch R image; r-ver provides base R across architectures
FROM rocker/r-ver:4.3.2

RUN apt-get update && apt-get install -y \
	build-essential \
	libcurl4-gnutls-dev \
	libxml2-dev \
	libssl-dev \
	libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \	
	libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
	cmake \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Set the working directory
WORKDIR /app

# Copy the R script to the container
# Copy full project structure needed for pipeline and report
COPY Code /app/Code
COPY Report /app/Report
COPY Results /app/Results
COPY tmp /app/tmp
COPY Data /app/Data
COPY run.sh /app/run.sh
COPY README.md /app/README.md

# Install packages
RUN Rscript Code/Install_packages.R

# Install TinyTeX for PDF rendering (disabled; HTML-only reports)
# RUN R -e "install.packages('tinytex', dependencies=TRUE); tinytex::install_tinytex(force=TRUE); tinytex::tlmgr_path()" && \
#     echo 'export PATH=/root/bin:/root/.TinyTeX/bin/aarch64-linux:/root/.TinyTeX/bin/x86_64-linux:$PATH' > /etc/profile.d/tinytex.sh
# ENV PATH="/root/bin:/root/.TinyTeX/bin/aarch64-linux:/root/.TinyTeX/bin/x86_64-linux:${PATH}"

# Optional: Install Quarto CLI (detect architecture and install matching package)
ARG QUARTO_VERSION=1.5.57
RUN arch=$(uname -m) && \
    if [ "$arch" = "x86_64" ] || [ "$arch" = "amd64" ]; then \
      wget -O /tmp/quarto.deb https://quarto.org/download/latest/quarto-linux-amd64.deb; \
    elif [ "$arch" = "aarch64" ] || [ "$arch" = "arm64" ]; then \
      wget -O /tmp/quarto.deb https://quarto.org/download/latest/quarto-linux-arm64.deb; \
    else \
      echo "Unsupported architecture: $arch" && exit 1; \
    fi && \
    apt-get update && apt-get install -y /tmp/quarto.deb && \
    rm /tmp/quarto.deb && rm -rf /var/lib/apt/lists/*

# Make the script executable
RUN chmod +x /app/run.sh

# Set the entrypoint to execute the script
ENTRYPOINT ["/app/run.sh"]
#ENTRYPOINT ["/bin/bash"]
