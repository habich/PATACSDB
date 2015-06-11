#!/usr/bin/python
from PATACSDB import app

if __name__ == "__main__":
    app.debug=True
    app.run(host='0.0.0.0')
