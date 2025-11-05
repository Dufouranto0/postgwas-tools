############################################################
# Dockerfile for postgwas-tools
############################################################

FROM python:3.10

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    libatlas-base-dev \
    wget \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy package source
COPY . .

# Install Python package
RUN pip install --upgrade pip \
    && pip install .

# Install PLINK 1.9
RUN wget -q https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20250819.zip \
    && unzip plink_linux_x86_64_20250819.zip \
    && rm plink_linux_x86_64_20250819.zip \
    && mv plink /app/plink \
    && chmod +x /app/plink

# Default command
CMD ["/bin/bash", "-i"]