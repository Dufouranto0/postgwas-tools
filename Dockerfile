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

# Install PLINK 1.9 safely into /app/plink
RUN wget -q https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20250819.zip \
    && mkdir -p /app/plink_bin \
    && unzip -o plink_linux_x86_64_20250819.zip -d /app/plink_bin \
    && rm plink_linux_x86_64_20250819.zip \
    && mv /app/plink_bin/plink /app/plink \
    && chmod +x /app/plink \
    && rm -rf /app/plink_bin

# Default command
CMD ["/bin/bash", "-i"]