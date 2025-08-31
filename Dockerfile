# ---- base ----
    FROM python:3.11-slim

    ENV PYTHONUNBUFFERED=1 PIP_NO_CACHE_DIR=1
    
    # System + R toolchain
    RUN apt-get update && apt-get install -y --no-install-recommends \
        curl tini \
        r-base r-base-dev \
        libxml2-dev libcurl4-openssl-dev libssl-dev \
        git make g++ pkg-config \
     && rm -rf /var/lib/apt/lists/*
    
    # Workdir
    WORKDIR /app
    
    # Python deps first (cache-friendly)
    COPY requirements.txt .
    RUN pip install --no-cache-dir -r requirements.txt
    
    # Copy the whole repo (contains your R package + Flask app + templates)
    COPY . .
    
    # Install your R package from this repo (since it has DESCRIPTION/NAMESPACE)
    RUN R -q -e "install.packages('remotes', repos='https://cloud.r-project.org')" \
     && R -q -e "remotes::install_local('.', upgrade='never', dependencies=TRUE)"
    
    # Runtime env
    ENV PORT=8080 \
        FLASK_SECRET_KEY=change-me \
        DEPTH_CSV=/mnt/data/depth.csv
    
    # Writable dirs your app uses
    RUN mkdir -p /mnt/data uploads
    
    EXPOSE 8080
    
    ENTRYPOINT ["/usr/bin/tini","--"]
    CMD ["gunicorn","-b","0.0.0.0:8080","microbiome_xai_chatbot:app","--workers","2","--threads","4","--timeout","120"]
    