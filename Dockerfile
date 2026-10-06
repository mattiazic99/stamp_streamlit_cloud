# STAMP interface only: the atlas is bundled; statistical experiments run separately.
FROM python:3.11-slim@sha256:a2bc8c35469b6fe37735f7c4dae39049470b2ce068e73f799c02452de31d24c6
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    HOME=/home/stamp \
    MPLCONFIGDIR=/tmp/matplotlib \
    XDG_CACHE_HOME=/tmp/cache \
    PYTHONPATH=/app:/app/gui
WORKDIR /app
COPY requirements-interface.lock.txt ./
RUN pip install --no-cache-dir -r requirements-interface.lock.txt
COPY stamp/ ./stamp/
COPY gui/ ./gui/
COPY .streamlit/ ./.streamlit/
COPY output/v10_complete/normalized/ ./output/v10_complete/normalized/
COPY output/v10_complete/sets/ ./output/v10_complete/sets/
COPY output/v10_complete/jaccard/ ./output/v10_complete/jaccard/
COPY output/v8_complete/normalized/ ./output/v8_complete/normalized/
COPY output/v8_complete/sets/ ./output/v8_complete/sets/
COPY output/v8_complete/jaccard/ ./output/v8_complete/jaccard/
# Documentation downloads only: no statistical CLI, raw inputs or reference runs.
COPY reproducibility/README.md ./reproducibility/README.md
COPY requirements-reproduction.lock.txt ./
RUN groupadd --gid 1000 stamp && useradd --uid 1000 --gid stamp --create-home stamp \
    && chmod -R a-w /app/output
USER stamp
EXPOSE 8501
HEALTHCHECK --interval=30s --timeout=5s --start-period=20s --retries=3 \
    CMD python -c "import urllib.request,sys; sys.exit(0 if urllib.request.urlopen('http://localhost:8501/_stcore/health', timeout=4).status==200 else 1)"
CMD ["streamlit", "run", "gui/main.py", "--server.port=8501", "--server.address=0.0.0.0", "--server.headless=true", "--browser.gatherUsageStats=false"]
