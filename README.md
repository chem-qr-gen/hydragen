# Chemquest

<https://chemquest.fly.dev>

An automatic question generator for organic chemistry.

## Requirements

- Python 3.6+
- Node.js
- [poetry](https://python-poetry.org/)

## Installation

```bash
$ poetry install
$ cd static
$ npm i
$ npm run build
$ cd ..
$ poetry run gunicorn -w 4 -b 0.0.0.0:5001 app:app
```

The website will be accessible at `localhost:5001` or `<YOUR_IP>:5001`.

## API Usage

The api for MS questions is accessible at <https://chemquest.fly.dev/ms_questions?id=1> or, if hosted locally, `<YOUR_IP>:5001/ms_questions?id=1`.

Parameters:

- `id`: Input the desired question id, or "random" for a redirect to a random question.