# Use the official R base image
FROM rocker/r-ubuntu:22.04

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

# Ensure Quarto CLI available for report rendering (optional)
# Uncomment to install Quarto in the image if needed
# RUN apt-get update && apt-get install -y wget && \
#     wget https://quarto.org/download/latest/quarto-linux-amd64.deb -O /tmp/quarto.deb && \
#     apt-get install -y /tmp/quarto.deb && \
#     rm /tmp/quarto.deb && rm -rf /var/lib/apt/lists/*

# Make the script executable
RUN chmod +x /app/run.sh

# Set the entrypoint to execute the script
ENTRYPOINT ["/app/run.sh"]
#ENTRYPOINT ["/bin/bash"]
