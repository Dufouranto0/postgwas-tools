############################################################
# Dockerfile for postgwas-tools
############################################################

FROM python:3.10

RUN apt-get update && apt-get install -y \
    build-essential \
    libatlas-base-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY . .

RUN pip install --upgrade pip \
    && pip install .

CMD ["/bin/bash", "-i"]