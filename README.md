# Hydragen

<https://hydragen.fly.dev>

An automatic question generator for organic chemistry.

## Requirements

- Python 3.11+
- Node.js 18+
- [poetry](https://python-poetry.org/)
- [Yarn](https://yarnpkg.com/)

## Installation

For development and testing (port 5000):

```bash
$ poetry install
$ cd chemquest_website/static
$ yarn install
$ yarn run sass-watch & # compile sass into css, run in background
$ yarn run watch
# In a separate console window (in the main folder):
$ poetry run flask --app=chemquest_website --debug run
```

For a production environment (port 5001):

```bash
$ poetry install
$ cd chemquest_website/static
$ yarn install
$ yarn run build
# In a separate console window (in the main folder):
$ poetry run gunicorn -w 4 -b 0.0.0.0:5001 chemquest_website:app
```
Known Issues

msChart theme remains the same if user is using system default and changes browser theme without reloading the page 
