# alignviz
An alignment vitalization app, for streamlit. 

## Running locally

```bash
pip install -r alignviz/requirements.txt
streamlit run alignviz/alignviz_app.py
```

## Deployment on SciLifeLab Serve

Every push to `main` builds a `linux/amd64` image and pushes it to
`ghcr.io/statisticalbiotechnology/alignviz` (see
[.github/workflows/docker-image-ghcr.yml](.github/workflows/docker-image-ghcr.yml)).
Each build gets a unique UTC timestamp tag, e.g. `20260820-160621`, alongside
`latest`; the exact tag to deploy is printed in the workflow run summary.

Serve needs to be able to pull the package anonymously. Check under repository
→ Packages → `alignviz` → Package settings that visibility is Public, and keep
it that way — Serve re-fetches the image at regular intervals.

To deploy or update the app, create/edit a Streamlit app at
[serve.scilifelab.se](https://serve.scilifelab.se) with port `8501` and the
image set to the timestamped tag. Serve caches by tag, so always point it at a
new timestamp rather than `latest` when updating.

Building and testing the image by hand:

```bash
docker build --platform linux/amd64 -t alignviz:test .
docker run --rm -it -p 8501:8501 alignviz:test   # then open http://localhost:8501
```
