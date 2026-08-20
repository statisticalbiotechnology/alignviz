# select starting image
FROM python:3.10-slim

# Create user name and home directory variables.
# The variables are later used as $USER and $HOME.
ENV USER=alignviz
ENV HOME=/home/$USER

# Add user to system
RUN useradd -m -u 1000 $USER

# Set working directory (this is where the code should go)
WORKDIR $HOME/app

# Update system and install dependencies.
# curl is needed by the HEALTHCHECK below.
RUN apt-get update && apt-get install --no-install-recommends -y \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Copy all files that the app needs.
# SciLifeLab Serve requires the main file to be called app.py and to sit in
# the working directory, so alignviz_app.py is renamed on the way in.
COPY alignviz/requirements.txt $HOME/app/requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

COPY alignviz/alignviz.py $HOME/app/alignviz.py
COPY alignviz/alignviz_app.py $HOME/app/app.py

USER $USER
EXPOSE 8501

HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health

ENTRYPOINT ["streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0", "--browser.gatherUsageStats=false"]
