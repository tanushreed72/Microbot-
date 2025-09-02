# ---- base ----
    FROM python:3.11-slim

    ENV PYTHONUNBUFFERED=1 PIP_NO_CACHE_DIR=1
    ENV MAKEFLAGS=-j4
    # Make reticulate use container Python
    ENV RETICULATE_PYTHON=/usr/local/bin/python3
    
    # System + R toolchain
    RUN apt-get update && apt-get install -y --no-install-recommends \
        curl tini \
        r-base r-base-dev \
        libxml2-dev libcurl4-openssl-dev libssl-dev \
        git make g++ pkg-config \
      && rm -rf /var/lib/apt/lists/*
    
    # Optional but speeds builds: preinstall igraph binary
    RUN apt-get update && apt-get install -y --no-install-recommends \
        r-cran-igraph \
      && rm -rf /var/lib/apt/lists/*
    
    # ---------- set workdir ----------
    WORKDIR /app
    
    # ---------- Python deps (cache-friendly) ----------
    COPY requirements.txt .
    RUN pip install --no-cache-dir -r requirements.txt
    
    # ---------- R package install (cache-friendly) ----------
    # Copy only the R package sources first so this layer caches unless R code changes
    # (Assuming your package is at repo root with DESCRIPTION/NAMESPACE and R/ folder)
    COPY DESCRIPTION NAMESPACE ./
    COPY R ./R
    # If you have inst/, man/, etc., include them too for a complete install:
    COPY inst ./inst
    COPY man ./man
    
    RUN R -q -e "options(warn=2); install.packages('remotes', repos='https://cloud.r-project.org')" \
     && R -q -e "options(warn=2); install.packages(c('plyr','pulsar','reticulate'), repos='https://cloud.r-project.org')" \
     && R -q -e "options(warn=2); remotes::install_local('.', upgrade='never', dependencies=TRUE)" \
     && R -q -e "library(InfIntE); cat('InfIntE installed: ', as.character(packageVersion('InfIntE')), '\n')"
    
    # ---------- App code ----------
    # Now copy the rest of the repo (Flask app, templates, static, etc.)
    COPY . .
    
    # Runtime env
    ENV PORT=8080 \
        FLASK_SECRET_KEY=change-me \
        DEPTH_CSV=/mnt/data/depth.csv
    
    # Writable dirs your app uses
    RUN mkdir -p /mnt/data uploads
    
    EXPOSE 8080
    
    ENTRYPOINT ["/usr/bin/tini","--"]
    CMD ["gunicorn","-b","0.0.0.0:8080","microbiome_xai_chatbot:app","--workers","2","--threads","4","--timeout","120"]
    