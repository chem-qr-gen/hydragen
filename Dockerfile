# syntax=docker/dockerfile:1
FROM node:lts

ENV POETRY_VENV=/opt/poetry_venv
RUN apt-get update && apt-get install -y python3 python3-pip python3-venv
RUN yarn set version stable

# Install poetry separated from system interpreter
RUN python3 -m venv $POETRY_VENV \
    && $POETRY_VENV/bin/pip install -U pip \
    && $POETRY_VENV/bin/pip install poetry

# Add `poetry` to PATH
ENV PATH="${PATH}:${POETRY_VENV}/bin"

WORKDIR /app
COPY . .

RUN poetry config virtualenvs.in-project true
RUN poetry install

WORKDIR /app/chemquest_website/static
RUN yarn install
RUN yarn run build

WORKDIR /app
EXPOSE 5001
CMD poetry run gunicorn -w 4 -b 0.0.0.0:5001 chemquest_website:app