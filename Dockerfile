# syntax=docker/dockerfile:1
FROM nikolaik/python-nodejs:python3.11-nodejs18-slim

USER pn
WORKDIR /usr/pn/app
COPY . .
RUN poetry install
WORKDIR /usr/pn/app/chemquest_website/static
RUN yarn install
RUN yarn run build
WORKDIR /usr/pn/app
EXPOSE 5001
CMD poetry run gunicorn -w 4 -b 0.0.0.0:5001 chemquest_website:app
