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
COPY Code/*.R /app/Code/

# Install packages
RUN Rscript Code/Install_packages.R

# Copy the R script to the container
COPY run.sh /app/run.sh

# Make the script executable
RUN chmod +x /app/run.sh

# Set the entrypoint to execute the script
ENTRYPOINT ["/app/run.sh"]
#ENTRYPOINT ["/bin/bash"]
