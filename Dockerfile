# docker build -t easypqp:latest .

FROM python:3.11-slim

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1

# Minimal build/runtime dependencies. Add or remove system packages as needed
# if package compilation fails (e.g., pyopenms may require extra libs).
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
   build-essential \
   gcc \
   git \
   cmake \
   swig \
   pkg-config \
   libxml2-dev \
   zlib1g-dev \
   libbz2-dev \
   liblzma-dev \
   libcurl4-openssl-dev \
   libssl-dev \
   # Runtime libraries required by pyopenms
   libglib2.0-0 \
   libgomp1 \
 && rm -rf /var/lib/apt/lists/*

# Upgrade packaging tools
RUN pip install --no-cache-dir --upgrade pip setuptools wheel

# Copy project into the image
WORKDIR /tmp/easypqp
COPY . /tmp/easypqp

# Install EasyPQP with all optional features by default
RUN pip install --no-cache-dir ".[all]"

# Cleanup sources
WORKDIR /
RUN rm -rf /tmp/easypqp

# Default command prints help
CMD ["easypqp","--help"]
