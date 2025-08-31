# ---- base ----
    FROM python:3.11-slim

    # Make logs unbuffered
    ENV PYTHONUNBUFFERED=1 \
        PIP_NO_CACHE_DIR=1
    
    # System deps (add more if wheels aren’t available on your platform)
    RUN apt-get update && apt-get install -y --no-install-recommends \
        curl tini && \
        rm -rf /var/lib/apt/lists/*
    
    # Create app dir
    WORKDIR /app
    
    # Copy and install deps first (better layer caching)
    COPY requirements.txt .
    RUN pip install -r requirements.txt
    
    # Copy your code
    COPY microbiome_xai_chatbot.py ./microbiome_xai_chatbot.py
    COPY templates ./templates
    
    # App runtime env (adjust as you like)
    ENV FLASK_SECRET_KEY=change-me \
        DEPTH_CSV=/mnt/data/depth.csv \
        PORT=8000
    
    # Create data dir your app expects/writes to
    RUN mkdir -p /mnt/data uploads
    
    # Expose port (informational)
    EXPOSE 8000
    
    # Use tini as init, run gunicorn in production
    ENTRYPOINT ["/usr/bin/tini", "--"]
    CMD ["gunicorn", "-b", "0.0.0.0:8000", "microbiome_xai_chatbot:app", "--workers", "2", "--threads", "4", "--timeout", "120"]
    