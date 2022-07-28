## Requirements

- Python 3.6+
- Node.js
- [poetry](https://python-poetry.org/)

## Installation

```bash
$ poetry install
$ cd static
$ npm i
$ cd ..
$ poetry run gunicorn -w 4 -b 0.0.0.0:5001 app:app
```

The website will be accessible at `localhost:5001` or `<YOUR_IP>:5001`.